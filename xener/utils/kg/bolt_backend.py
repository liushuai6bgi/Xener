from typing import Literal
import pandas as pd
import scipy.sparse as sp

from neo4j import GraphDatabase

from .base import KGBackend


class KG_Neo4j_BoltBackend(KGBackend):
    """Neo4j Bolt protocol backend implementation."""

    def __init__(self, url, auth):
        """
        Initialize Neo4j database connection.

        Args:
            url: Neo4j database address.
            auth: Neo4j database username and password.
        """
        self.driver = GraphDatabase.driver(url, auth=auth)
        
    def _build_cypher_gene2celltype_path(
        self,
        resolution: Literal["Cell", "Tissue"],
        organ: str | list = None,
        genes: list[str] = None,
        celltypes: list[str] = None,
    ) -> str:
        """
        Build Cypher query for gene-to-celltype paths in the knowledge graph.
        Can filter by cell type and gene lists.

        Args:
            resolution: Cell type resolution.
            organ: Organ name.
            genes: List of gene names.
            celltypes: List of cell type names.

        Returns:
            Cypher query string.
        """
        if resolution == "Cell":
            cypher = "MATCH path = (a:Gene)-[b:marker_of]->(c:Ontology) "
        elif resolution == "Tissue":
            cypher = (
                "MATCH path = (a:Gene)-[b:marker_of]->(c:Ontology)-[d:ontology_relation1..2]-(e:Ontology) "
            )
        else:
            raise ValueError("resolution must be Cell or Tissue")

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
        """
        Build Cypher query for cell type-to-cell type paths.

        Args:
            celltypes: List of cell type names.
            max_step: Maximum number of relationship steps.

        Returns:
            Cypher query string.
        """
        relationships = ["is_a", "part_of", "intersection_of", "develops_from"]
        cypher = f"""
            MATCH path = (a:Ontology)-[b*1..{max_step}]->(c:Ontology)
            WHERE a.Name IN {celltypes}
            AND c.Name IN {celltypes}
            AND ALL(rel IN b WHERE rel.relation IN {relationships}) """

        cypher += "RETURN DISTINCT path"
        return cypher

    def get_genecount_kg(self, celltype: str) -> int:
        """
        Count how many genes are connected to this cell type.

        Args:
            celltype: Cell type name.

        Returns:
            Number of connected genes.
        """
        cypher = (
            'MATCH (a:Gene)-->(:Ontology{Name:"' + celltype + '"})' + " RETURN COUNT(a)"
        )
        with self.driver.session() as session:
            results = session.run(cypher)
            data = results.value()[0]
            return data

    def get_celltypecount_kg(self, gene: str) -> int:
        """
        Count how many cell types are connected to this gene.

        Args:
            gene: Gene name.

        Returns:
            Number of connected cell types.
        """
        cypher = (
            'MATCH (a:Gene{Name:"' + gene + '"})-->(:Ontology)' + " RETURN COUNT(a)"
        )
        with self.driver.session() as session:
            results = session.run(cypher)
            data = results.value()[0]
            return data

    def get_gene2celltype_kg(
        self,
        homolo_nodes: list[str] = None,
        organ: str = None,
        resolution: Literal["Cell", "Tissue"] = "Cell",
        candidate_type: list[str] = None,
    ) -> tuple[list[str], list[str], sp.csr_matrix]:
        """
        Retrieve gene-to-celltype relationships.

        Args:
            homolo_nodes: List of homology nodes.
            organ: Organ name.
            resolution: Annotation granularity ('Cell' or 'Tissue').
            candidate_type: Candidate cell types.

        Returns:
            source: Deduplicated gene list.
            target: Deduplicated cell type list.
            matrix: Adjacency sparse matrix.
        """
        cypher = self._build_cypher_gene2celltype_path(
            resolution, organ, homolo_nodes, candidate_type
        )
        source_nodes = []
        target_nodes = []
        edges = []
        with self.driver.session() as session:
            results = session.run(cypher)
            for result in results.graph().relationships:
                source_nodes.append(result.nodes[0]["Name"])
                target_nodes.append(result.nodes[-1]["Name"])
                # 提取边的起点、终点以及权重
                edges.append(
                    (
                        result.nodes[0]["Name"],
                        target_nodes[-1],
                        result._properties["relation_confidence"],
                    )
                )
        # 构造邻接矩阵
        cellType_nodes = list(set(target_nodes))
        cellType2idx = {node: idx for idx, node in enumerate(cellType_nodes)}

        gene_nodes = list(set(source_nodes)) if homolo_nodes is None else homolo_nodes
        gene2idx = {node: idx for idx, node in enumerate(gene_nodes)}

        gene2celltype_matrix = sp.lil_matrix((len(gene_nodes), len(cellType_nodes)))
        for source_node, target_node, v in edges:
            gene2celltype_matrix[gene2idx[source_node], cellType2idx[target_node]] = v

        return gene_nodes, cellType_nodes, gene2celltype_matrix.tocsr()

    def get_celltype2celltype_kg(
        self, nodes: list[str], symmetric: bool = False, max_step: int = 1
    ) -> tuple[sp.csr_matrix, list[str]]:
        """
        Get the relationship matrix between cell types.

        Args:
            nodes: List of cell type nodes.
            symmetric: Whether to return a symmetric matrix.
            max_step: Maximum number of steps.

        Returns:
            celltype2celltype_matrix: Relationship matrix between cell types.
            new_nodes: New cell type node list. Cell types that appear in the matrix.
                       The count may differ from the input.
        """
        new_nodes = []
        rows, cols = [], []
        new_nodes2idx = dict()
        idx = 0

        cypher = self._build_cypher_celltype2celltype_path(nodes, max_step)
        with self.driver.session() as session:
            results = session.run(cypher)
            for result in results.graph().relationships:
                source_node = result.nodes[0]["Name"]
                if source_node not in new_nodes2idx.keys():
                    new_nodes.append(source_node)
                    new_nodes2idx[source_node] = idx
                    idx += 1
                rows.append(new_nodes2idx[source_node])

                target_node = result.nodes[-1]["Name"]
                if target_node not in new_nodes2idx.keys():
                    new_nodes.append(target_node)
                    new_nodes2idx[target_node] = idx
                    idx += 1
                cols.append(new_nodes2idx[target_node])

                if symmetric:
                    rows.append(new_nodes2idx[target_node])
                    cols.append(new_nodes2idx[source_node])

        matrix = sp.coo_matrix(
            ([1] * len(rows), (rows, cols)), shape=(len(new_nodes), len(new_nodes))
        )
        matrix = matrix.tocsr()
        # 去除重复边
        matrix.data[matrix.data > 1] = 1
        return matrix, new_nodes

    def get_species_organ_cell(self) -> pd.DataFrame:
        """
        Return the species-organ-celltype relationship table.

        Returns:
            DataFrame with columns: species, organ, cell.
        """
        species_organ_cell = []
        cypher = "MATCH (a:Gene)-[b]->(c:Ontology) RETURN DISTINCT c"
        with self.driver.session() as session:
            for value in session.run(cypher).values():
                cell_name = value[0]["Name"]
                organ_names = value[0]["Organ"].split("|")
                species_name = value[0]["Species_type"]
                for organ_name in organ_names:
                    species_organ_cell.append([species_name, organ_name, cell_name])
        species_organ_cell = pd.DataFrame(
            species_organ_cell, columns=["species", "organ", "cell"]
        ).drop_duplicates()
        return species_organ_cell

    def close(self):
        self.driver.close()
