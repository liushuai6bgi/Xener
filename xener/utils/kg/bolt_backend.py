from typing import Literal
import time
import pandas as pd
import scipy.sparse as sp

from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable, SessionExpired, TransientError

from ..logger import logger
from .base import KGBackend


# ---------------------------------------------------------------------------
# Retry / transient-error handling
# ---------------------------------------------------------------------------

_RETRYABLE_ERRORS = (ServiceUnavailable, SessionExpired, TransientError)
_MAX_RETRIES = 3
_RETRY_BASE_DELAY = 2.0  # seconds  →  2s, 4s, 8s exponential backoff


def _retry_on_transient(label, func):
    """Call *func()*, retrying on transient Neo4j errors with exponential backoff.

    *func* is called once per attempt (no arguments) and must return the query
    result.  Each attempt opens a **fresh** driver session so that a stuck /
    stale connection from a previous attempt never contaminates the retry.

    The three built-in retryable types cover the user's observed symptoms:

    * ``ServiceUnavailable`` – routing service not responding (the pingable
      host that nevertheless refuses Bolt work).
    * ``SessionExpired``   – server closed the session (e.g. connection-pool
      pressure).
    * ``TransientError``   – generic catch-all for "please retry" errors
      signalled by the server.
    """
    last_exc = None
    for attempt in range(_MAX_RETRIES):
        try:
            return func()
        except _RETRYABLE_ERRORS as exc:
            last_exc = exc
            if attempt < _MAX_RETRIES - 1:
                wait = _RETRY_BASE_DELAY * (2 ** attempt)
                logger.warning(
                    "KG transient error during %s, retry %s/%s in %.1fs: %s",
                    label, attempt + 1, _MAX_RETRIES, wait, exc,
                )
                time.sleep(wait)
    raise last_exc


# ---------------------------------------------------------------------------
# Backend
# ---------------------------------------------------------------------------

class KG_Neo4j_BoltBackend(KGBackend):
    """Neo4j Bolt protocol backend implementation."""

    def __init__(self, url, auth):
        """Initialize Neo4j database connection.

        Args:
            url: Neo4j database address.
            auth: Neo4j database username and password.

        Connection-pool tuning
        ----------------------
        The defaults shipped by the ``neo4j`` driver are geared toward
        multi-tenant application servers (``max_connection_pool_size=400``,
        no connection lifetime ceiling).  A batch-annotation pipeline is a
        **single** client issuing sequential read-only queries — a pool that
        large only pressures the server's routing layer without benefiting
        throughput.

        The values below bound the pool to 50 connections, cap individual
        connection lifetime at 1 hour (preventing the server-side idle-timeout
        disconnect that surfaces as ``ServiceUnavailable``), and enable TCP
        keep-alive so intermediate firewalls / NAT tables don't drop the
        connection silently.
        """
        logger.info("KG_Neo4j_BoltBackend connecting to %s", url)
        self.driver = GraphDatabase.driver(
            url,
            auth=auth,
            max_connection_pool_size=50,
            connection_acquisition_timeout=30,
            max_transaction_retry_time=30,
            keep_alive=True,
            max_connection_lifetime=3600,   # 1 hour — recycle before server idle-kill
        )

    # ------------------------------------------------------------------
    # Cypher builders (unchanged logic)
    # ------------------------------------------------------------------

    def _build_cypher_gene2celltype_path(
        self,
        organ: str | list = None,
        genes: list[str] = None,
        celltypes: list[str] = None,
    ) -> str:
        cypher = "MATCH path = (a:Gene)-[b:marker_of]->(c:Ontology) "

        WHERE = []
        if genes:
            WHERE.append(f" a.Name IN {genes} ")
        if celltypes:
            WHERE.append(f" c.Name IN {celltypes} ")
        if organ:
            organs = ["Unknown"]
            if isinstance(organ, str):
                organs.append(organ)
            if isinstance(organ, list):
                organs.extend(organ)
            WHERE.append(
                f" ANY(organ IN split(c.Organ, '|') WHERE organ IN {organs}) "
            )
        if len(WHERE) > 0:
            cypher += "WHERE " + " AND ".join(WHERE)

        cypher += "RETURN DISTINCT path "
        return cypher

    def _build_cypher_celltype2celltype_path(
        self, celltypes: list[str], max_step: int = 1
    ) -> str:
        relationships = ["is_a", "part_of", "intersection_of", "develops_from"]
        cypher = f"""
            MATCH path = (a:Ontology)-[b*1..{max_step}]->(c:Ontology)
            WHERE a.Name IN {celltypes}
            AND c.Name IN {celltypes}
            AND ALL(rel IN b WHERE rel.relation IN {relationships}) """

        cypher += "RETURN DISTINCT path"
        return cypher

    # ------------------------------------------------------------------
    # Query helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _execute(tx, cypher):
        """Run a Cypher statement inside a transaction function.

        Returns the raw ``neo4j.Result`` so callers can call
        ``.value()``, ``.values()``, ``.graph()``, etc.
        """
        return tx.run(cypher)

    # ------------------------------------------------------------------
    # KG queries  (each wrapped in transient-error retry)
    # ------------------------------------------------------------------

    def get_genecount_kg(self, celltype: str) -> int:
        cypher = (
            'MATCH (a:Gene)-->(:Ontology{Name:"' + celltype + '"})'
            + " RETURN COUNT(a)"
        )
        t0 = time.time()

        def _query():
            with self.driver.session() as session:
                result = session.execute_read(self._execute, cypher)
                values = result.value()
                if not values:
                    logger.warning(
                        "KG get_genecount_kg(celltype=%s) -> empty result, returning 0",
                        celltype,
                    )
                    return 0
                return values[0]

        data = _retry_on_transient("get_genecount_kg", _query)
        logger.info(
            "KG get_genecount_kg(celltype=%s) -> %s in %.3fs",
            celltype, data, time.time() - t0,
        )
        return data

    def get_celltypecount_kg(self, gene: str) -> int:
        cypher = (
            'MATCH (a:Gene{Name:"' + gene + '"})-->(:Ontology)'
            + " RETURN COUNT(a)"
        )
        t0 = time.time()

        def _query():
            with self.driver.session() as session:
                result = session.execute_read(self._execute, cypher)
                values = result.value()
                if not values:
                    logger.warning(
                        "KG get_celltypecount_kg(gene=%s) -> empty result, returning 0",
                        gene,
                    )
                    return 0
                return values[0]

        data = _retry_on_transient("get_celltypecount_kg", _query)
        logger.info(
            "KG get_celltypecount_kg(gene=%s) -> %s in %.3fs",
            gene, data, time.time() - t0,
        )
        return data

    def get_gene2celltype_kg(
        self,
        homolo_nodes: list[str] = None,
        organ: str = None,
        candidate_type: list[str] = None,
    ) -> tuple[list[str], list[str], sp.csr_matrix]:
        cypher = self._build_cypher_gene2celltype_path(
            organ, homolo_nodes, candidate_type
        )
        logger.info(
            "KG get_gene2celltype_kg request: %s homolos, organ=%s, %s candidate_types",
            len(homolo_nodes) if homolo_nodes else 0,
            organ,
            len(candidate_type) if candidate_type else 0,
        )

        def _tx_work(tx):
            """Transaction function: consumes graph results fully inside tx."""
            result = tx.run(cypher)
            src = []
            tgt = []
            edg = []
            for rel in result.graph().relationships:
                src.append(rel.nodes[0]["Name"])
                tgt.append(rel.nodes[-1]["Name"])
                edg.append((
                    rel.nodes[0]["Name"],
                    tgt[-1],
                    rel._properties["relation_confidence"],
                ))
            return src, tgt, edg

        t0 = time.time()

        def _query():
            with self.driver.session() as session:
                return session.execute_read(_tx_work)

        source_nodes, target_nodes, edges = _retry_on_transient(
            "get_gene2celltype_kg", _query,
        )
        cypher_elapsed = time.time() - t0

        # Build sparse adjacency matrix
        cellType_nodes = list(set(target_nodes))
        cellType2idx = {node: idx for idx, node in enumerate(cellType_nodes)}

        gene_nodes = list(set(source_nodes)) if homolo_nodes is None else homolo_nodes
        gene2idx = {node: idx for idx, node in enumerate(gene_nodes)}

        gene2celltype_matrix = sp.lil_matrix((len(gene_nodes), len(cellType_nodes)))
        for source_node, target_node, v in edges:
            gene2celltype_matrix[gene2idx[source_node], cellType2idx[target_node]] = v

        csr = gene2celltype_matrix.tocsr()
        logger.info(
            "KG get_gene2celltype_kg done: cypher=%.2fs, %s genes, %s celltypes, "
            "%s edges, matrix=%s, nnz=%s",
            cypher_elapsed, len(gene_nodes), len(cellType_nodes),
            len(edges), csr.shape, csr.nnz,
        )
        if csr.nnz == 0:
            logger.warning(
                "KG get_gene2celltype_kg returned an empty matrix (organ=%s, %s homolos).",
                organ,
                len(homolo_nodes) if homolo_nodes else 0,
            )
        return gene_nodes, cellType_nodes, csr

    def get_celltype2celltype_kg(
        self, nodes: list[str], symmetric: bool = False, max_step: int = 1
    ) -> tuple[sp.csr_matrix, list[str]]:
        cypher = self._build_cypher_celltype2celltype_path(nodes, max_step)
        logger.info(
            "KG get_celltype2celltype_kg request: %s nodes, symmetric=%s, max_step=%s",
            len(nodes), symmetric, max_step,
        )

        def _tx_work(tx):
            """Transaction function: consumes graph results fully inside tx."""
            result = tx.run(cypher)
            _rows = []
            _cols = []
            _new_nodes = []
            _new_nodes2idx = {}
            _idx = 0
            for rel in result.graph().relationships:
                source_node = rel.nodes[0]["Name"]
                if source_node not in _new_nodes2idx:
                    _new_nodes.append(source_node)
                    _new_nodes2idx[source_node] = _idx
                    _idx += 1
                _rows.append(_new_nodes2idx[source_node])

                target_node = rel.nodes[-1]["Name"]
                if target_node not in _new_nodes2idx:
                    _new_nodes.append(target_node)
                    _new_nodes2idx[target_node] = _idx
                    _idx += 1
                _cols.append(_new_nodes2idx[target_node])

                if symmetric:
                    _rows.append(_new_nodes2idx[target_node])
                    _cols.append(_new_nodes2idx[source_node])

            return _rows, _cols, _new_nodes, len(_new_nodes)

        t0 = time.time()

        def _query():
            with self.driver.session() as session:
                return session.execute_read(_tx_work)

        rows, cols, new_nodes, n_nodes = _retry_on_transient(
            "get_celltype2celltype_kg", _query,
        )
        cypher_elapsed = time.time() - t0

        matrix = sp.coo_matrix(
            ([1] * len(rows), (rows, cols)), shape=(n_nodes, n_nodes)
        )
        matrix = matrix.tocsr()
        matrix.data[matrix.data > 1] = 1  # deduplicate
        logger.info(
            "KG get_celltype2celltype_kg done: cypher=%.2fs, %s nodes, "
            "%s edges, matrix=%s, nnz=%s",
            cypher_elapsed, n_nodes, len(rows), matrix.shape, matrix.nnz,
        )
        return matrix, new_nodes

    def get_species_organ_cell(self) -> pd.DataFrame:
        cypher = (
            "MATCH (a:Gene)-[b]->(c:Ontology) "
            "RETURN DISTINCT a.Species,c.Organ,c.Name"
        )
        logger.info("KG get_species_organ_cell ...")

        def _tx_work(tx):
            """Transaction function: iterate all result rows."""
            records = []
            for species_name, organ_names, cell_name in tx.run(cypher).values():
                records.append((species_name, organ_names, cell_name))
            return records

        t0 = time.time()

        def _query():
            with self.driver.session() as session:
                return session.execute_read(_tx_work)

        records = _retry_on_transient("get_species_organ_cell", _query)

        species_organ_cell = []
        for species_name, organ_names, cell_name in records:
            species_name = species_name.strip().replace(" ", "_")
            organ_name_list = organ_names.split("|")
            for organ_name in organ_name_list:
                species_organ_cell.append([species_name, organ_name, cell_name])

        species_organ_cell = pd.DataFrame(
            species_organ_cell, columns=["species", "organ", "cell"]
        ).drop_duplicates()
        logger.info(
            "KG get_species_organ_cell done in %.2fs: %s rows, %s species, %s organs",
            time.time() - t0, len(species_organ_cell),
            species_organ_cell["species"].nunique(),
            species_organ_cell["organ"].nunique(),
        )
        return species_organ_cell

    def close(self):
        self.driver.close()
