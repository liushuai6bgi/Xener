from abc import ABC, abstractmethod
from typing import Literal
import pandas as pd
import scipy.sparse as sp


class KGBackend(ABC):
    """Abstract base class for knowledge graph backends, defining a unified interface."""

    @abstractmethod
    def get_genecount_kg(self, celltype: str) -> int:
        """Return the number of genes connected to this cell type."""
        pass

    @abstractmethod
    def get_celltypecount_kg(self, gene: str) -> int:
        """Return the number of cell types connected to this gene."""
        pass

    @abstractmethod
    def get_gene2celltype_kg(
        self,
        homolo_nodes: list[str] = None,
        organ: str = None,
        resolution: Literal["Cell", "Tissue"] = "Cell",
        candidate_type: list[str] = None,
    ) -> tuple[list[str], list[str], sp.csr_matrix]:
        """Retrieve gene-to-celltype relationships."""
        pass

    @abstractmethod
    def get_celltype2celltype_kg(
        self, nodes: list[str], symmetric: bool = False, max_step: int = 1
    ) -> tuple[sp.csr_matrix, list[str]]:
        """Get the relationship matrix between cell types."""
        pass

    @abstractmethod
    def get_species_organ_cell(self) -> pd.DataFrame:
        """Return the species-organ-celltype relationship table."""
        pass

    def close(self):
        """Close the connection; subclasses may override."""
        pass
