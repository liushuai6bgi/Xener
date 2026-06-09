from typing import Literal
import time
import pandas as pd
import scipy.sparse as sp
import requests

from ..logger import logger
from .base import KGBackend


class KG_HttpBackend(KGBackend):
    """HTTP API backend implementation.

    Accesses functionality via REST API without relying on Cypher queries.
    """

    def __init__(self, url: str):
        """
        Initialize HTTP backend connection.

        Args:
            url: HTTP service root address (e.g., http://localhost:8080).
            auth: (username, password) tuple, default None.
        """
        self.url = url.rstrip("/")
        self._session = requests.Session()

        logger.info('KG_HttpBackend connecting to %s', self.url)
        t0 = time.time()
        self.species_organ_cell = self.get_species_organ_cell()
        logger.info('KG_HttpBackend species_organ_cell loaded in %.2fs (%s rows)',
                    time.time() - t0, len(self.species_organ_cell))

        # 统计所有可用的组织
        self.available_organ_dict = {
            organ.lower(): organ for organ in self.species_organ_cell["organ"].unique()
        }
        self.available_organ_set = set(self.available_organ_dict.keys())

    def _request(self, method: str, path: str, **kwargs) -> requests.Response:
        """Send an HTTP request."""
        url = f"{self.url}{path}"
        t0 = time.time()
        try:
            response = self._session.request(method, url, timeout=30, **kwargs)
        except requests.RequestException as e:
            logger.error('KG HTTP request failed: %s %s (elapsed=%.2fs): %s',
                         method, path, time.time() - t0, e)
            raise
        elapsed = time.time() - t0
        if not response.ok:
            detail = response.json().get("detail", response.text) if response.content else response.text
            logger.error('KG HTTP %s %s returned %s in %.2fs: %s',
                         method, path, response.status_code, elapsed, detail)
            raise requests.HTTPError(detail, response=response)
        logger.debug('KG HTTP %s %s %s in %.2fs (%s bytes)',
                     method, path, response.status_code, elapsed, len(response.content))
        return response



    def get_genecount_kg(self, celltype: str) -> int:
        """Return the number of genes connected to this cell type."""
        # GET /api/kg/genecount/{celltype}
        resp = self._request("GET", f"/api/kg/genecount/{celltype}")
        count = resp.json()["count"]
        logger.info('KG get_genecount_kg(celltype=%s) -> %s', celltype, count)
        return count

    def get_celltypecount_kg(self, gene: str) -> int:
        """Return the number of cell types connected to this gene."""
        # GET /api/kg/celltypecount/{gene}
        resp = self._request("GET", f"/api/kg/celltypecount/{gene}")
        count = resp.json()["count"]
        logger.info('KG get_celltypecount_kg(gene=%s) -> %s', gene, count)
        return count

    def get_gene2celltype_kg(
        self,
        homolo_nodes: list[str] = None,
        organ: str = None,
        candidate_type: list[str] = None,
    ) -> tuple[list[str], list[str], sp.csr_matrix]:
        """Retrieve gene-to-celltype relationships."""
        # POST /api/kg/gene2celltype
        params = {}
        if homolo_nodes:
            params["homolo_nodes"] = homolo_nodes
        if organ:
            params["organ"] = organ
        if candidate_type:
            params["candidate_type"] = candidate_type
        logger.info('KG get_gene2celltype_kg request: %s homolos, organ=%s, %s candidate_types',
                    len(homolo_nodes) if homolo_nodes else 0, organ,
                    len(candidate_type) if candidate_type else 0)
        t0 = time.time()
        resp = self._request("POST", "/api/kg/gene2celltype", json=params)
        data = resp.json()
        gene_nodes = data["genes"]
        cellType_nodes = data["celltypes"]
        matrix_data = data["matrix"]
        matrix = sp.csr_matrix(
            (matrix_data["data"], matrix_data["indices"], matrix_data["indptr"]),
            shape=tuple(matrix_data["shape"])
        )
        nnz = matrix.nnz
        logger.info('KG get_gene2celltype_kg done in %.2fs: %s genes, %s celltypes, matrix %s, nnz=%s',
                    time.time() - t0, len(gene_nodes), len(cellType_nodes), matrix.shape, nnz)
        if nnz == 0:
            logger.warning('KG get_gene2celltype_kg returned an empty matrix (organ=%s, %s homolos).',
                           organ, len(homolo_nodes) if homolo_nodes else 0)
        return gene_nodes, cellType_nodes, matrix

    def get_celltype2celltype_kg(
        self, nodes: list[str], symmetric: bool = False, max_step: int = 1
    ) -> tuple[sp.csr_matrix, list[str]]:
        """Get the relationship matrix between cell types."""
        # POST /api/kg/celltype2celltype
        params = {"nodes": nodes, "symmetric": symmetric, "max_step": max_step}
        logger.info('KG get_celltype2celltype_kg request: %s nodes, symmetric=%s, max_step=%s',
                    len(nodes), symmetric, max_step)
        t0 = time.time()
        resp = self._request("POST", "/api/kg/celltype2celltype", json=params)
        data = resp.json()
        new_nodes = data["nodes"]
        matrix_data = data["matrix"]
        matrix = sp.csr_matrix(
            (matrix_data["data"], matrix_data["indices"], matrix_data["indptr"]),
            shape=tuple(matrix_data["shape"])
        )
        logger.info('KG get_celltype2celltype_kg done in %.2fs: %s nodes, matrix %s, nnz=%s',
                    time.time() - t0, len(new_nodes), matrix.shape, matrix.nnz)
        return matrix, new_nodes

    def get_species_organ_cell(self) -> pd.DataFrame:
        """Return the species-organ-celltype relationship table."""
        # GET /api/kg/species-organ-cell
        logger.info('KG get_species_organ_cell ...')
        t0 = time.time()
        resp = self._request("GET", "/api/kg/species-organ-cell")
        data = resp.json()['data']
        df = pd.DataFrame(data, columns=["species", "organ", "cell"]).drop_duplicates()
        logger.info('KG get_species_organ_cell done in %.2fs: %s rows, %s species, %s organs',
                    time.time() - t0, len(df), df['species'].nunique(), df['organ'].nunique())
        return df

    def close(self):
        self._session.close()
