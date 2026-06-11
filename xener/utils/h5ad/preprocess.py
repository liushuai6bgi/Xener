import os
import sys
import gc

import scanpy as sc
import anndata as ad
from scipy.sparse import issparse

from ..logger import logger

def read_h5ad(path:str) -> sc.AnnData:
    """
    Read an h5ad file and make sure adata.X is sparse.

    Args:
        path: Path to the h5ad file.

    Returns:
        AnnData object.
    """
    adata = sc.read_h5ad(path)
    raw_available = hasattr(adata, 'raw') and adata.raw is not None
    
    # the old version of scanpy, rebuild the raw
    if raw_available and '_index' in adata.raw.var.columns:
        raw_var = adata.raw.var.set_index('_index').copy()
        raw_var.index = raw_var.index.astype(str)
        raw_var.index.name = None
        adata.raw = ad.AnnData(X=adata.raw.X, var=raw_var)
        logger.info("it's saved by old version of scanpy, rebuild the raw.")

    # check gene names
    if raw_available:
        try:
            adata.raw.var_names.astype(float)# if it's gene names, it will raise a ValueError
            raise TypeError('adata.raw.var_names isn\'t gene names!')
        except ValueError:
            pass
    if not raw_available:
        try:
            adata.var_names.astype(float)# if it's gene names, it will raise a ValueError
            raise TypeError('adata.var_names isn\'t gene names, cann\'t find available gene names!')
        except ValueError:
            pass

    return adata

def write_h5ad(adata:sc.AnnData, path:str):
    """
    Write an AnnData object to an h5ad file.

    Args:
        adata: AnnData object to write.
        path: Output file path.
    """
    # if hasattr(adata, 'raw') and adata.raw is not None:
    #     logger.warning('adata.raw will be deleted in write_h5ad!')
    #     del adata.raw
    adata.write_h5ad(path, compression='gzip')

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
    # raw counts
    if adata.X.max() > 1e5:
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.filter_cells(adata, min_genes=200)
        logger.info('normalize_total adata')
        sc.pp.normalize_total(adata)
    # need log1p
    if adata.X.max() > 10:
        logger.info('log1p adata')
        sc.pp.log1p(adata)
    # update raw
    if hasattr(adata, 'raw') and adata.raw is not None:
        raw_adata = adata.raw.to_adata()
        if raw_adata.X.max() > 1e5:
            logger.info('normalize_total raw_adata')
            sc.pp.normalize_total(raw_adata)
        if raw_adata.X.max() > 10:
            if 'log1p' in raw_adata.uns:
                del raw_adata.uns['log1p']
            logger.info('log1p raw_adata')
            sc.pp.log1p(raw_adata)
        adata.raw = raw_adata
    elif issparse(adata.X):
        # logger.info('set adata.raw to adata')
        # adata.raw = adata.copy()
        logger.info('skip setting adata.raw; use sparse adata.X directly to avoid large raw copies')
    return adata

def process(adata:sc.AnnData, force_HVG:bool=False, n_top_genes:int=2000):
    threshold_genes = 4000
    adata_genes = adata.shape[1]
    adata = quality_control(adata)
    logger.info(f'gene counts that [adata:threshold]=[{adata_genes}:{threshold_genes}].')
    if force_HVG or adata_genes > threshold_genes:
        logger.info(f'highly_variable_genes[{n_top_genes}].')
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,       # 指定选取2000个高变基因
            subset=True            # 自动筛选并保留高变基因
        )
    else:
        logger.info('skip highly_variable_genes.')

    if 'X_pca' not in adata.obsm.keys():
        if 'highly_variable' in adata.var.columns:
            adata.var['highly_variable'] = adata.var['highly_variable'].astype(bool)
        logger.info('pca')
        sc.pp.pca(adata)
    if 'connectivities' not in adata.obsp.keys():
        logger.info('neighbors')
        sc.pp.neighbors(adata)
    if 'leiden' not in adata.obs.keys():
        # 使用参数flavor=igraph会导致代码报奇怪的错：ValueError: high is out of bounds for int32
        logger.info('leiden')
        sc.tl.leiden(adata, flavor="leidenalg", n_iterations=-1)
    if 'X_umap' not in adata.obsm.keys():
        logger.info('umap')
        sc.tl.umap(adata)
    return adata

if __name__ == '__main__':
    adata = read_h5ad(sys.argv[1])
    logger.info(f'Input: {adata}')
    dir = os.path.dirname(sys.argv[1])
    adata = quality_control(adata)
    logger.info(f'Output: {adata}')
    write_h5ad(adata, os.path.join(dir, 'adata_p.h5ad'))
