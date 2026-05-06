from typing import Literal
import pandas as pd
import scipy.sparse as sp
import requests

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

        self.species_organ_cell = self.get_species_organ_cell()

        # 统计所有可用的组织
        self.available_organ_dict = {
            organ.lower(): organ for organ in self.species_organ_cell["organ"].unique()
        }
        self.available_organ_set = set(self.available_organ_dict.keys())

    def _request(self, method: str, path: str, **kwargs) -> requests.Response:
        """Send an HTTP request."""
        url = f"{self.url}{path}"
        response = self._session.request(method, url, timeout=30, **kwargs)
        response.raise_for_status()
        return response



    def get_genecount_kg(self, celltype: str) -> int:
        """Return the number of genes connected to this cell type."""
        # GET /api/kg/genecount/{celltype}
        resp = self._request("GET", f"/api/kg/genecount/{celltype}")
        return resp.json()["count"]

    def get_celltypecount_kg(self, gene: str) -> int:
        """Return the number of cell types connected to this gene."""
        # GET /api/kg/celltypecount/{gene}
        resp = self._request("GET", f"/api/kg/celltypecount/{gene}")
        return resp.json()["count"]

    def get_gene2celltype_kg(
        self,
        homolo_nodes: list[str] = None,
        organ: str = None,
        resolution: Literal["Cell", "Tissue"] = "Cell",
        candidate_type: list[str] = None,
    ) -> tuple[list[str], list[str], sp.csr_matrix]:
        """Retrieve gene-to-celltype relationships."""
        # GET /api/kg/gene2celltype
        params = {"resolution": resolution}
        if homolo_nodes:
            params["homolo_nodes"] = ",".join(homolo_nodes)
        if organ:
            params["organ"] = organ
        if candidate_type:
            params["candidate_type"] = ",".join(candidate_type)
        resp = self._request("GET", "/api/kg/gene2celltype", params=params)
        data = resp.json()
        gene_nodes = data["genes"]
        cellType_nodes = data["celltypes"]
        matrix_data = data["matrix"]
        matrix = sp.csr_matrix(
            (matrix_data["data"], matrix_data["indices"], matrix_data["indptr"]),
            shape=tuple(matrix_data["shape"])
        )
        return gene_nodes, cellType_nodes, matrix

    def get_celltype2celltype_kg(
        self, nodes: list[str], symmetric: bool = False, max_step: int = 1
    ) -> tuple[sp.csr_matrix, list[str]]:
        """Get the relationship matrix between cell types."""
        # POST /api/kg/celltype2celltype
        params = {"nodes": ",".join(nodes), "symmetric": symmetric, "max_step": max_step}
        resp = self._request("GET", "/api/kg/celltype2celltype", params=params)
        data = resp.json()
        new_nodes = data["nodes"]
        matrix_data = data["matrix"]
        matrix = sp.csr_matrix(
            (matrix_data["data"], matrix_data["indices"], matrix_data["indptr"]),
            shape=tuple(matrix_data["shape"])
        )
        return matrix, new_nodes

    def get_species_organ_cell(self) -> pd.DataFrame:
        """Return the species-organ-celltype relationship table."""
        # GET /api/kg/species-organ-cell
        resp = self._request("GET", "/api/kg/species-organ-cell")
        data = resp.json()
        return pd.DataFrame(data, columns=["species", "organ", "cell"]).drop_duplicates()

    def close(self):
        self._session.close()
