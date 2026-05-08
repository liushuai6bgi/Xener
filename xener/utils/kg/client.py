import yaml
from typing import Literal

from .base import KGBackend
from .bolt_backend import KG_Neo4j_BoltBackend
from .http_backend import KG_HttpBackend


class KGClient(KGBackend):
    """Knowledge graph client factory class.

    Automatically selects Bolt or HTTP backend implementation based on configuration.
    """

    def __init__(self, url:str, usr=None, pwd=None):
        """Initialize KGClient."""
        backend_type = "http" if url.startswith("http") else "bolt"
        if backend_type == "bolt":
            backend = KG_Neo4j_BoltBackend(url=url, auth=(usr, pwd))
        elif backend_type == "http":
            backend = KG_HttpBackend(url=url)
        else:
            raise ValueError(f"Unknown backend type: {backend_type}")
        self._backend = backend
        self.species_organ_cell = self._backend.get_species_organ_cell()

        # 统计所有可用的组织
        self.available_organ_dict = {
            organ.lower(): organ for organ in self.species_organ_cell["organ"].unique()
        }
        self.available_organ_set = set(self.available_organ_dict.keys())
    def check_organ(self, organ: str | list) -> str | list:
        """Check and normalize organ input, auto-handling capitalization."""
        if organ is None:
            return organ
        if isinstance(organ, str):
            organ_low = organ.lower()
            if organ_low in self.available_organ_set:
                organ = self.available_organ_dict[organ_low]
            else:
                organ = None
        elif isinstance(organ, list):
            organ = [self.check_organ(x) for x in organ]
        return organ
    @staticmethod
    def init_from_yaml(cls, yaml_file: str) -> "KGClient":
        """
        Initialize from a YAML configuration file.

        Config file format:
            KG_url: bolt://localhost:7687  # or http://localhost:7474
            KG_usr: username
            KG_pwd: password
        """
        with open(yaml_file, "r", encoding="utf8") as f:
            config = yaml.load(f, Loader=yaml.FullLoader)

        url = config["KG_url"]
        usr = config.get("KG_usr")
        pwd = config.get("KG_pwd")
        return cls(url, usr, pwd)

    @staticmethod
    def from_bolt(url: str, auth: tuple) -> "KGClient":
        """Create client from a Bolt connection."""
        return KGClient(backend=KG_Neo4j_BoltBackend(url=url, auth=auth))

    @staticmethod
    def from_http(url: str, auth: tuple = None) -> "KGClient":
        """Create client from an HTTP connection."""
        return KGClient(backend=KG_HttpBackend(url=url, auth=auth))

    def get_genecount_kg(self, celltype: str) -> int:
        return self._backend.get_genecount_kg(celltype)

    def get_celltypecount_kg(self, gene: str) -> int:
        return self._backend.get_celltypecount_kg(gene)

    def get_gene2celltype_kg(
        self,
        homolo_nodes: list[str] = None,
        organ: str = None,
        candidate_type: list[str] = None,
    ) -> tuple[list[str], list[str], "sp.csr_matrix"]:
        return self._backend.get_gene2celltype_kg(
            homolo_nodes, organ, candidate_type
        )

    def get_celltype2celltype_kg(
        self, nodes: list[str], symmetric: bool = False, max_step: int = 1
    ) -> tuple["sp.csr_matrix", list[str]]:
        return self._backend.get_celltype2celltype_kg(nodes, symmetric, max_step)

    def get_species_organ_cell(self) -> "pd.DataFrame":
        return self._backend.get_species_organ_cell()

    def close(self):
        self._backend.close()