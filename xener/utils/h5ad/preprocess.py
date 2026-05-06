import os
import sys

import scanpy as sc

from ..logger import logger

def read_h5ad(path:str) -> sc.AnnData:
    """
    Read an h5ad file.

    Args:
        path: Path to the h5ad file.

    Returns:
        AnnData object.
    """
    adata = sc.read_h5ad(path)
    if hasattr(adata, 'raw'):
        del adata.raw
    return adata

def write_h5ad(adata:sc.AnnData, path:str):
    """
    Write an AnnData object to an h5ad file.

    Args:
        adata: AnnData object to write.
        path: Output file path.
    """
    if hasattr(adata, 'raw'):
        del adata.raw
    adata.write_h5ad(path)

def quality_control(adata:sc.AnnData) -> sc.AnnData:
    """
    Perform quality control on single-cell data.

    Args:
        adata: AnnData object containing single-cell data.

    Returns:
        Quality-controlled AnnData object.
    """
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=200)
    adata.X[adata.X < 0] = 1e-10
    sc.pp.normalize_total(adata)
    if 'log1p' not in adata.uns or adata.X.max() > 1000:
        logger.info('Applying log1p transformation')
        sc.pp.log1p(adata)
    return adata

if __name__ == '__main__':
    adata = read_h5ad(sys.argv[1])
    logger.info(f'Input: {adata}')
    dir = os.path.dirname(sys.argv[1])
    adata = quality_control(adata)
    logger.info(f'Output: {adata}')
    write_h5ad(adata, os.path.join(dir, 'adata_p.h5ad'))
