import os
import time
import psutil
from pathlib import Path
from typing import Literal

import yaml
import numpy as np
import scanpy as sc
import pandas as pd
import networkx as nx
import scipy.sparse as sp
import scanpy.external as sce

from .utils import logger
from .utils import (
    deliver_weight_on_graph_max, deliver_weight_on_graph_sum, build_graph_from_adjust_matrix, 
    tempered_softmax_contributions)#, check_celltype)
from .utils.blast import makeblastdb, get_blastdb, blastp
from .utils.seq import extract_fasta_by_name, translate2protein
from .utils.kg import KGClient
from .utils.h5ad import read_h5ad, write_h5ad, quality_control
from .config import _XENER_DATA_DIR, _ensure_data

class Xener:
    @staticmethod
    def read_h5ad(path:str):
        return read_h5ad(path)
    @staticmethod
    def write_h5ad(adata:sc.AnnData, path:str):
        return write_h5ad(adata, path)
    
    @classmethod
    def init_from_yaml(cls, yaml_file:str, **kwargs):
        f =  open(yaml_file, 'r', encoding='utf-8')
        config = yaml.load(f, Loader=yaml.FullLoader)
        f.close()
        config.update(kwargs)

        kg_kwargs = {
            'url': config['KG_url'],
            'usr': config.get('KG_usr', None),
            'pwd': config.get('KG_pwd', None),
        }
        return cls(kg_kwargs, blastdb_path=config['blastdb_path'], 
                   blastp_result_path=config.get('blastp_result_path', None))

    def __init__(self, kg_kwargs:dict=None,
                 blastdb_path:str=None, blastp_result_path=None):
        '''
        Initialize the Xener object.
        Parameters:
            kg_kwargs (dict[str,str]): Neo4j connection information
            blastdb_path (str): Path to the existing BLAST database
            blastp_result_path (str): Path to BLASTP result files
        '''
        if blastdb_path is None:
            self.blastdb_path = _XENER_DATA_DIR / Path('blastdb/prot')
            if not self.blastdb_path.exists():
                _ensure_data(DATA_KEY='blastdb.zip')
                
        else:
            self.blastdb_path = Path(blastdb_path)
            assert self.blastdb_path.exists(), f"blastdb_path does not exist! Please verify the provided parameters"
        self.blastdb = get_blastdb(self.blastdb_path)
        self.blastp_result_path = blastp_result_path
        if kg_kwargs is None:
            kg_kwargs = {'url': 'https://xenor.dcs.cloud'}
        self.KG = KGClient(**kg_kwargs)
        logger.info('Xener initialized!')

    def run_from_yaml(self, yaml_file:str, Resource_occupancy_record:dict=None, save=True, configs:dict={}) -> dict:
        '''
        Encapsulates the entire workflow for testing.
        Parameters:
            yaml_file (str): Path to the YAML configuration file
            Resource_occupancy_record (dict): Record of resource usage. If not None, should contain keys: steps, cpu_time, and memory, with corresponding values as lists.
            save (bool): Whether to save intermediate results (used for subsequent analysis)
            configs (dict): Parameters in the configuration file will be updated, but the file content remains unchanged.
        Returns:
            dict: cluster2celltype mapping
            dict: cluster2max_initweight_celltype mapping
        '''
        process = psutil.Process(os.getpid())
        def recorder(phase:str):
            Resource_occupancy_record['steps'].append(phase)
            Resource_occupancy_record['cpu_time'].append(time.process_time())
            Resource_occupancy_record['memory'].append(process.memory_info().rss / 1024 / 1024)
            
        if Resource_occupancy_record:
            recorder('init')

        _ensure_data(DATA_KEY='default_run_conf.yaml')
        with open(_XENER_DATA_DIR / 'default_run_conf.yaml', 'r', encoding='utf-8') as f:
            config = yaml.load(f, Loader=yaml.FullLoader)
        with open(yaml_file, 'r', encoding='utf-8') as f:
            config_from_yaml = yaml.load(f, Loader=yaml.FullLoader)
        config.update(config_from_yaml)
        config.update(configs)

        if 'outdir' not in config:
            raise 'outdir not found in yaml file!'
        os.makedirs(config['outdir'], exist_ok=True)
        self.outdir = Path(config['outdir'])

        if isinstance(config.get('marker_gene', None), str) and os.path.exists(config['marker_gene']):
            marker_gene = pd.read_csv(config['marker_gene'], sep=None, engine='python', header=0)
        elif config.get('marker_gene', None) is None or config['marker_gene']:
            logger.info('generating marker_gene ...')
            adata = read_h5ad(config['non_model_h5ad'])
            logger.info('adata: %s', adata)
            if Resource_occupancy_record:
                recorder('load h5ad')

            marker_gene = self.get_markers(adata, config.get('cluster_key', None))
            if save:
                marker_gene_path = self.outdir / 'marker_gene.zip'
                marker_gene.to_csv(marker_gene_path, index=False)
                config['marker_gene'] = str(marker_gene_path)
                logger.info('marker_gene saved to %s', marker_gene_path)
            logger.info('marker_gene.shape: %s', marker_gene.shape)
            if Resource_occupancy_record:
                recorder('generate marker_gene')
        
        if isinstance(config.get('marker_weight', None), str) and os.path.exists(config['marker_weight']):
            marker_weight = pd.read_csv(config['marker_weight'], sep=None, engine='python', header=0)
        elif config.get('marker_weight', None) is None or config['marker_weight']:
            logger.info('generating marker_weight ...')
            marker_weight = self.get_gene_weight(marker_gene, config.get('marker_weight_method', 'prod'))
            if save:
                marker_weight_path = self.outdir / 'marker_weight.zip'
                marker_weight.to_csv(marker_weight_path, index=False)
                config['marker_weight'] = str(marker_weight_path)
                logger.info('marker_weight saved to %s', marker_weight_path)
            logger.info('marker_weight.shape: %s', marker_weight.shape)
            if Resource_occupancy_record:
                recorder('generate marker_weight')

        if isinstance(config.get('gene_homolo_weight', None), str) and os.path.exists(config['gene_homolo_weight']):
            gene_homolo_weight = pd.read_csv(config['gene_homolo_weight'], sep=None, engine='python', header=0)
        elif config.get('gene_homolo_weight', None) is None or config['gene_homolo_weight']:
            logger.info('generating gene_homolo_weight ...')
            if not config.get('is_protein_fasta', True):
                config['non_model_fasta'] = translate2protein(config['non_model_fasta'])
            as_homolo_weight_key = config.get('as_homolo_weight_key', 'pident')
            gene_homolo_weight = self.mapping(marker_weight, config['non_model_fasta'], config['model_species'], self.outdir, as_homolo_weight_key, config.get('gene_name_prefix', None), config.get('kv_blastp_args', None))
            if save:
                gene_homolo_weight_path = self.outdir / 'gene_homolo_weight.zip'
                gene_homolo_weight.to_csv(gene_homolo_weight_path, index=False)
                config['gene_homolo_weight'] = str(gene_homolo_weight_path)
                logger.info('gene_homolo_weight saved to %s', gene_homolo_weight_path)
            logger.info('gene_homolo_weight.shape: %s', gene_homolo_weight.shape)
            if Resource_occupancy_record:
                recorder('generate gene_homolo_weight')

        regenerate = True
        if isinstance(config.get('topk_markers', None), str) and os.path.exists(config['topk_markers']):
            topk_markers = pd.read_csv(config['topk_markers'], sep=None, engine='python', header=0)
            gene_count = topk_markers[['group', 'gene']].drop_duplicates()
            regenerate = gene_count.shape[0] != len(gene_count['group'].unique()) * config['top_num']
        if config.get('topk_markers', None) is None or config['topk_markers'] or regenerate:
            logger.info('generating topk_markers ...')
            topk_markers = self.get_topk_gene(gene_homolo_weight, config['top_num'], config.get('multihomolo', True))
            if save:
                topk_markers_path = self.outdir / 'topk_markers.zip'
                topk_markers.to_csv(topk_markers_path, index=False)
                config['topk_markers'] = str(topk_markers_path)
                logger.info('topk_markers saved to %s', topk_markers_path)
            logger.info('topk_markers.shape: %s', topk_markers.shape)
            if Resource_occupancy_record:
                recorder('generate topk_markers')

        annotation_info_path = self.outdir / 'annotation'
        cluster2celltype, cluster2max_initweight_celltype, celltype_weight, cluster_celltype_ann, homolo2celltype = self.cell_annotation(
            topk_markers, annotation_info_path, config.get('organ', None), 
            config.get('threshold', None), config.get('candidate_annotation', None), 
            mode=config.get('mode', 'path'), decay_factor=config.get('decay_factor', 0.7))
        if Resource_occupancy_record:
            recorder('cell_annotation')
        # Cell type information
        if save:
            celltype_weight_path = self.outdir / 'celltype_weight.zip'
            celltype_weight.to_csv(celltype_weight_path, index=False)
            homolo2celltype.to_csv(self.outdir / 'homolo2celltype.zip', index=False)
            logger.info('celltype_weight saved to %s', celltype_weight_path)
        logger.info('celltype_weight.shape: %s', celltype_weight.shape)
        if cluster_celltype_ann.shape[0] > 0:
            cluster_celltype_ann_path = self.outdir / 'cluster_celltype_ann.zip'
            cluster_celltype_ann.to_csv(cluster_celltype_ann_path)
        if Resource_occupancy_record:
            recorder('save celltype_weight')

        with open(self.outdir / 'config.yaml', 'w') as f:
            yaml.dump(config, f, default_flow_style=False)

        return cluster2celltype, cluster2max_initweight_celltype
    
    def get_markers(self, adata:sc.AnnData, cluster_key:str=None,
                    preprocess:bool=False, batch_key:str=None) -> pd.DataFrame:
        '''
        Preprocess and analyze single-cell data, returning marker genes.
        Parameters:
            adata (sc.AnnData): Object containing single-cell data, including gene expression matrix and associated metadata
            preprocess (bool, optional): Whether to perform quality control preprocessing. Defaults to False
            batch_key (str, optional): Batch key. Defaults to None; if provided, batch effect correction will be applied
        Returns:
            pd.DataFrame: Result of scanpy.get.rank_genes_groups_df
        '''
        if preprocess:
            adata = quality_control(adata)
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=2000,       # Select 2000 highly variable genes
                subset=True            # Automatically filter and retain highly variable genes
            )
            sc.pp.pca(adata)
        if batch_key:
            # bbknn can run PCA by default, and its output as an optimized version of neighbors is consistent with sc.pp.neighbors
            # bbknn does not modify expression data; since we don't directly use expression data afterward, using bbknn should be fine
            sce.pp.bbknn(adata, batch_key=batch_key)
        if adata.obs[cluster_key].dtype != 'str':
            adata.obs[cluster_key] = adata.obs[cluster_key].astype('str')
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, key_added = "rank_genes_groups", pts=True)
        # Get gene ranking table = number of clusters * number of genes
        markers = sc.get.rank_genes_groups_df(adata, group=None, key='rank_genes_groups')
        return markers.drop(columns=['scores'])
    
    def get_gene_weight(self, markers:pd.DataFrame, method:Literal['prod','sum']='prod'):
        # Determine data source: scanpy/Seurat
        col_maps = {
            'scanpy':{'pct_nz_group':'pct1',    'pct_nz_reference':'pct2',  'pvals_adj':'pv', 'logfoldchanges':'logfc'},
            'Seurat':{'pct.1':'pct1',           'pct.2':'pct2',             'p_val_adj':'pv', 'avg_log2FC':'logfc'}
        }
        marker_type = None
        if 'names' in markers.columns:# scanpy
            marker_type = 'scanpy'
        elif 'gene' in markers.columns:# Seurat
            marker_type = 'Seurat'
        else:
            raise ValueError('Invalid marker type')
        markers.rename(columns=col_maps[marker_type], inplace=True)

        # Get positive values
        # def to_norm(ser:pd.Series):
        #     return (ser - ser.min()) / (ser.max() - ser.min())

        # Check for infinite values in logfc
        if any(np.isinf(markers['logfc'])):
            replace_inf = markers[np.isfinite(markers['logfc'])]['logfc'].abs().max() * 10
            logger.info('replace_inf logfc: %s', replace_inf)
            markers.loc[markers['logfc'] == np.inf, 'logfc'] = replace_inf
            markers.loc[markers['logfc'] == -np.inf, 'logfc'] = -replace_inf
        # Check for zero values in pvals_adj
        if any(markers['pv']==0):
            replace_zero = markers[markers['pv'] != 0]['pv'].min() * 0.5
            logger.info('replace_zero pvals_adj: %s', replace_zero)
            markers.loc[markers['pv'] == 0, 'pv'] = replace_zero

        markers['logfc'] = np.clip(markers['logfc'], 0, markers['logfc'].max())
        if any(np.isinf(markers['logfc'])):
            logger.warning('logfc has inf values, please check the data')
        markers['pts_delta'] = markers['pct1'] - markers['pct2']
        markers['pts_delta'] = 1 + np.clip(markers['pts_delta'], 0, markers['pts_delta'].max())
        if any(np.isinf(markers['pts_delta'])):
            logger.warning('pts_delta has inf values, please check the data')
        markers['-log10_pvals_adj'] = markers['pv'].map(lambda x: -np.log10(x))
        if any(np.isinf(markers['-log10_pvals_adj'])):
            logger.warning('-log10_pvals_adj has inf values, please check the data')
        if method == 'prod':
            markers['weight'] = markers['logfc'] * markers['pts_delta'] * markers['-log10_pvals_adj']
        elif method == 'sum':
            markers['weight'] = markers['logfc'] + markers['pts_delta'] + markers['-log10_pvals_adj']
        markers.loc[markers['weight'] < 0, 'weight'] = 0

        if marker_type == 'scanpy':
            markers.rename(columns={'names': 'gene'}, inplace=True)
        elif marker_type == 'Seurat':
            markers.rename(columns={'cluster': 'group'}, inplace=True)
        markers = markers[['gene', 'weight', 'group']]
        if any(np.isinf(markers['weight'])):
            logger.warning('weight has inf values, please check the data')
        return markers

    def get_topk_gene(self, markers:pd.DataFrame, k:int=10, multihomolo:bool=True) -> pd.DataFrame:
        # Compatibility with old version
        if 'weight' not in markers.columns:
            if 'energy' in markers.columns:
                markers['weight'] = markers['energy']
            else:
                raise ValueError('No weight column found in markers')
        # Start getting top genes
        ## Step 1: Group by group and gene, get unique weight for each gene
        grouped = markers.groupby(['group', 'gene'], as_index=False)['weight'].first()
        ## Step 2: Sort by weight in descending order within each group
        grouped['rank'] = grouped.groupby('group')['weight'].rank(method='first', ascending=False)
        ## Step 3: Take top_num genes for each group
        top_gene = grouped[grouped['rank'] <= k].drop(columns=['rank'])
        ## Step 4: Keep gene to homolo mapping
        result = pd.merge(top_gene, markers[['gene', 'homolo', 'homolo_weight']].drop_duplicates(), on='gene', how='left')
        result.sort_values(by=['group', 'weight'], inplace=True, ascending=[True, False])

        # Single homolo
        if not multihomolo:
            result = result.loc[result.groupby(['group', 'gene'])['homolo_weight'].idxmax()]
        
        for group_id in result['group'].unique():
            group_data = result[result['group'] == group_id]
            gene_count = len(group_data['gene'].unique())
            homolo_count = len(group_data['homolo'].unique())
            logger.info('group %s has %s gene and %s homolo', group_id, gene_count, homolo_count)
        result['homolo_weight'] = result['homolo_weight'] / result.groupby(['group', 'gene'])['homolo_weight'].transform('sum')
        return result

    def mapping(self, markers:pd.DataFrame, non_model_fasta:str, 
            model_species:list[str], outdir:str, as_homolo_weight_key='pident', gene_name_prefix:str=None,
            num_threads:int=None, pident=60, evalue=0.05, bitscore=200, **kv_blastp_args) -> pd.DataFrame:
        '''
        Run BLAST to find homologous genes and integrate them into markers.
        Parameters:
            markers (pd.DataFrame): Relationship between each cluster and each gene
            non_model_fasta (str): Path to the non-model species FASTA file
            model_species (list[str]): List of model species
        Returns:
            pd.DataFrame: Relationship between marker genes and their homologous genes, with weight column named homolo_weight
        '''
        # Compatibility with old version
        if 'weight' not in markers.columns:
            if 'energy' in markers.columns:
                markers['weight'] = markers['energy']
            else:
                raise ValueError('No weight column found in markers')
        origin_group_num = len(markers['group'].unique())
        # Start finding homologs
        marker_list = markers['gene'].unique().tolist()
        logger.info('mapping %s genes from %s', len(marker_list), model_species)
        seq_file = extract_fasta_by_name(outdir, marker_list, non_model_fasta, gene_name_prefix)
        blastp_result = []
        for species in model_species:
            non_model_name = os.path.basename(non_model_fasta).split('.')[0]
            blast_result_name = f'blastp4{non_model_name}2{species}.zip'
            if self.blastp_result_path:
                blast_result_path = os.path.join(self.blastp_result_path, blast_result_name)
                if os.path.exists(blast_result_path):
                    # If full gene alignment results exist, read directly
                    data = pd.read_csv(blast_result_path, compression='zip')
                    data = data[data['qseqid'].isin(marker_list)]
                    blastp_result.append(data)
            else:
                logger.info('run blastp for %s', species)
                # Align only the extracted genes
                blastp_result.append(
                    blastp(seq_file, self.blastdb[species], outdir / f'blastp_{species}.zip', num_threads, **kv_blastp_args)
                )
        os.remove(seq_file)
        blast_result = pd.concat(blastp_result)

        # Filter alignment results
        condition = True
        if pident is not None:
            condition &= blast_result['pident'] > pident # Percentage of identical matches
        if evalue is not None:
            condition &= blast_result['evalue'] < evalue # Expect value
        if bitscore is not None:
            condition &= blast_result['bitscore'] > bitscore # Alignment score
        blast_result = blast_result[condition]

        blast_result = blast_result[['qseqid', 'sseqid', as_homolo_weight_key]]
        blast_result.columns = ['gene', 'homolo', 'homolo_weight']

        # Merge homologs
        merged_data = pd.merge(markers, blast_result, on='gene', how='inner')
        merged_data = merged_data[['group', 'gene', 'homolo', 'weight', 'homolo_weight']]
        # Normalize homolo_weight
        # merged_data['homolo_weight'] = merged_data['homolo_weight'] / merged_data.groupby('gene')['homolo_weight'].transform('sum')
        if origin_group_num!=len(merged_data['group'].unique()):
            logger.warning('%s groups merged to %s groups in mapping', origin_group_num, len(merged_data["group"].unique()))
        blast_result_with_weight = merged_data.sort_values(by=['group', 'weight'], ascending=[True, False])
        return blast_result_with_weight

    def cell_annotation(self, blast_result, 
            outdir:Path, organ=None, threshold:int=None, 
            resolution:Literal['Cell', 'Tissue']='Cell', mode:Literal['node','path']='node', decay_factor:float=0.7) -> tuple[dict[str, str], pd.DataFrame, pd.DataFrame]:
        '''
        Annotate each cell cluster and save celltype weight information for each cluster.
        Parameters:
            blast_result (pandas.DataFrame): BLAST result. Key `group` contains cluster information. Key `gene` contains gene names. Key `homolo` contains BLAST homologous genes. Key `weight` contains gene weights. Key `pct` contains confidence between homologous gene and gene.
            resolution (str): Annotation granularity
            organ (str, optional): Organ name. Defaults to None
            threshold (int, optional): Threshold for filtering cell types. Defaults to None
            candidate_annotation (list[str], optional): Candidate cell types. Defaults to None
        Returns:
            cluster2celltype (dict[str, str]): Mapping from cluster to celltype
            cluster2max_initweight_celltype (dict[str, str]): Mapping from cluster to celltype with maximum initial weight
            celltype_weight (pandas.DataFrame): Cell type information
            cluster_celltype_ann (pandas.DataFrame): Mapping from predicted type to candidate type, indexed by cluster
            homolo2celltype (pandas.DataFrame): Relationship between homologous genes and cell types, columns: 'homolo', 'celltype', 'homolo2celltype'
        '''
        outdir = Path(outdir)
        os.makedirs(outdir, exist_ok=True)
        # Validate organ
        if organ is not None:
            checked_organ = self.KG.check_organ(organ)
            if checked_organ is None:
                logger.warning('organ %s is unavailable, set to None.', organ)
            organ = checked_organ
            
        # Lists to save all celltype weights for each cluster, these three lists should have equal length
        cluster_list, celltype_list, weight_list, init_weight_list, genecount_graph_list, genecount_KG_list, homolo2celltype_list =\
        [], [], [], [], [], [], []

        cluster_celltype_ann = pd.DataFrame(index=blast_result['group'].unique(),columns=['celltype', 'ann_celltype'])
        cluster2celltype = dict()
        cluster2max_initweight_celltype = dict()
        for group in blast_result['group'].unique():
            logger.info('processing cluster %s', group)
            cluster2celltype[group], cluster2max_initweight_celltype[group], celltypes, weights, init_weights,\
                  genecount_graph, genecount_KG, homolo2celltype, gene2celltype_matrix = self.cell_annotation_cluster_singletype(
                blast_result[blast_result['group'] == group], resolution, 
                f'cluster_{group}', outdir, organ, candidate_annotation=None, 
                threshold=threshold, mode=mode, decay_factor=decay_factor)

            cluster_list += [group] * len(celltypes)
            celltype_list.extend(celltypes)
            weight_list.extend(weights)
            init_weight_list.extend(init_weights)
            genecount_graph_list.extend(genecount_graph)
            genecount_KG_list.extend(genecount_KG)
            homolo2celltype_list.extend(homolo2celltype)
        celltype_weight = pd.DataFrame(
            {'cluster': cluster_list, 
             'celltype': celltype_list, 
             'genecount_graph': genecount_graph_list,
             'genecount_KG': genecount_KG_list,
             'weight': weight_list,
             'init_weight':init_weight_list}
             )
        homolo2celltype = pd.DataFrame(homolo2celltype_list, columns=['homolo', 'celltype', 'homolo2celltype']).drop_duplicates()
        return cluster2celltype, cluster2max_initweight_celltype, celltype_weight, cluster_celltype_ann, homolo2celltype

    def cell_annotation_cluster_singletype(self, blast_result:pd.DataFrame, 
            resolution:Literal['Cell', 'Tissue'], output_prefix, outdir:Path, 
            organ=None, candidate_annotation:list[str]=None, 
            threshold:int=None,return_genecount_graph=True, return_genecount_KG=True,
            mode:Literal['node','path']='node', decay_factor:float=0.9
        ) -> tuple[str, list[str], list[float]]:
        '''
        Annotate a single cluster with cell types.
        Parameters:
            blast_result (pd.DataFrame): BLAST result for a single cluster. Key `group` contains cluster information. Key `gene` contains gene names. Key `homolo` contains BLAST homologous genes. Key `weight` contains gene weights. Key `pct` contains confidence between homologous gene and gene.
            resolution (str): Annotation granularity
            organ (str, optional): Organ name. Defaults to None
            candidate_annotation (list[str], optional): Candidate cell types to constrain knowledge graph query results. Defaults to None
            threshold (int, optional): Threshold for filtering candidate types using normal distribution of cell type weights
        Returns:
            celltype (str): Cell type; when mode is set to 'path', the return value may contain multiple cell types separated by '>' and difficult-to-distinguish types connected by '|'
            max_initweight_celltype (str): Same as celltype when mode is 'node'; when mode is 'path', this is the cell type with the maximum initial weight
            celltype_list (list[str]): All possible cell types
            weight_list (list[float]): Cell type weights, sorted in descending order; these are aggregated results
            initial_weight_list (list[float]): Cell type weights before aggregation
            genecount_graph (list[int]): Number of genes for each cell type in the annotation graph
            genecount_KG (list[int]): Number of genes for each cell type in the knowledge graph
            homolo2celltype_list (list[tuple[str,str,float]]): Records relationship between homologous genes and cell types, in order: homolo_node, celltype_node, weight
            gene2celltype_matrix (sp.csr_matrix): Sparse matrix representing gene-to-celltype relationships
        '''
        outdir = Path(outdir)
        os.makedirs(outdir, exist_ok=True)
        celltype = None
        output_prefix = output_prefix.replace(':', '')

        # Query knowledge graph based on homologs
        homolo_nodes = blast_result['homolo'].unique().tolist()
        _, celltype_nodes, homolo2celltype_matrix = self.KG.get_gene2celltype_kg(
            homolo_nodes, organ, resolution, candidate_type=candidate_annotation)# Do not limit when querying knowledge graph

        # Decay the relationship between homologous genes and cell types; the more cell types connected to a homologous gene, the less distinguishable they are, and the more decay is applied
        sum_result = homolo2celltype_matrix.sum(axis=1)
        # Check for homolos not appearing in the knowledge graph
        index = np.where(sum_result == 0)
        for i in index[0].tolist():  # Ignore zero-weight edges
            logger.warning('hmolo[%s] of organ[%s] not in kg', homolo_nodes[i], organ)
        logger.info('total %s %% homolos of organ[%s] not in kg', round(len(index[0])*100/len(homolo_nodes), 1), organ)
        max_result = homolo2celltype_matrix.max(axis=1)
        divisor = homolo2celltype_matrix.shape[1] * max_result.tocsr()
        # Using 1.1 to avoid decay_coefficient being set to 0 when there is only one candidate type
        decay_coefficient = 1.1 - sp.csr_matrix(sum_result) / divisor
        homolo2celltype_matrix = homolo2celltype_matrix.multiply(decay_coefficient).tocsr()

        # Record the relationship between homologous genes and cell types
        homolo2celltype_list = []
        coo_matrix = homolo2celltype_matrix.tocoo()
        for row, col, weight in zip(coo_matrix.row, coo_matrix.col, coo_matrix.data):
            if weight == 0:  # Ignore zero-weight edges
                continue
            homolo_node = homolo_nodes[row]
            celltype_node = celltype_nodes[col]
            homolo2celltype_list.append([homolo_node, celltype_node, weight])

        sub_g_homolo2celltype = build_graph_from_adjust_matrix(
            homolo2celltype_matrix,
            ['homolo_'+i for i in homolo_nodes],
            celltype_nodes,
            'relation_confidence',
            source_attr=[{'name':i,'weight':0,'type':'homolo'} for i in homolo_nodes],
            target_attr=[{'name':i,'weight':0,'type':'celltype'} for i in celltype_nodes],
            mode='n'
        )

        # Build relationship between genes and homologous genes
        gene_nodes = blast_result['gene'].unique().tolist()
        gene2homolo_matrix = sp.lil_matrix((len(gene_nodes), len(homolo_nodes)))
        for _, row in blast_result.iterrows():
            gene_idx = gene_nodes.index(row['gene'])
            homolo_idx = homolo_nodes.index(row['homolo'])
            gene2homolo_matrix[gene_idx, homolo_idx] = row['homolo_weight']
        # Get gene weights
        gene_weight = []
        for gene_node in gene_nodes:
            gene_weight.append(
                blast_result[blast_result['gene'] == gene_node]['weight'].iloc[0]
            )
        logger.info('%s: min=%s, max=%s', output_prefix, min(gene_weight), max(gene_weight))
        sub_g_gene2homolo = build_graph_from_adjust_matrix(
            gene2homolo_matrix.tocsr(), 
            ['gene_'+i for i in gene_nodes], 
            ['homolo_'+i for i in homolo_nodes], 
            'relation_confidence',
            source_attr=[{'name': gene,'weight':weight,'type':'gene'} for gene, weight in zip(gene_nodes, gene_weight)],
            target_attr=[{'name': i,'weight':0,'type':'homolo'} for i in homolo_nodes],
            mode='n'
        )
        # Merge the two graphs
        sub_g_gene2celltype = nx.compose(sub_g_homolo2celltype, sub_g_gene2homolo)

        sub_g_gene2celltype = deliver_weight_on_graph_sum(sub_g_gene2celltype, 'weight', 'relation_confidence')
        # Extract gene to celltype relationships
        gene2celltype_matrix = sp.lil_matrix((len(gene_nodes), len(celltype_nodes)))
        for celltype in celltype_nodes:
            if 'recorder' not in sub_g_gene2celltype.nodes[celltype].keys():
                # Did not participate in weight propagation
                logger.info('%s: Did not participate in weight propagation %s', output_prefix, celltype)
                continue
            for key, value in sub_g_gene2celltype.nodes[celltype]['recorder'].items():
                if sub_g_gene2celltype.nodes[key]['type'] == 'gene':
                    gene_idx = gene_nodes.index(sub_g_gene2celltype.nodes[key]['name'])
                    celltype_idx = celltype_nodes.index(celltype)
                    gene2celltype_matrix[gene_idx, celltype_idx] = value
        gene2celltype_matrix = gene2celltype_matrix.tocsr()

        # Save cell type information
        weight_distribution = {sub_g_gene2celltype.nodes[node_name]['name']: sub_g_gene2celltype.nodes[node_name]
                               for node_name in sub_g_gene2celltype.nodes if sub_g_gene2celltype.nodes[node_name]['type'] == 'celltype'}
        sorted_nodes = sorted(weight_distribution.items(), key=lambda item: item[1]['weight'], reverse=True)

        # Return values
        celltypes = []
        weight = []
        for node, attrs in sorted_nodes:
            celltypes.append(attrs['name'])
            weight.append(attrs['weight'])

        # Filter
        weight_array = np.array(weight)
        # Apply normal distribution to weight_array
        if weight_array.shape[0] == 0:
            logger.warning('weight_array has no element')
        elif weight_array.std() == 0: # Standard deviation is 0
            logger.warning('%s: weight_array.std() == 0', output_prefix)
        else:
            weight_array = (weight_array - weight_array.mean()) / weight_array.std()

        if threshold is not None:
            filtered_nodes = [celltypes[idx] for idx in np.where(weight_array > threshold)[0].tolist()]
        else:
            filtered_nodes = celltypes

        if mode == 'node':
            if len(filtered_nodes) == 0:# No candidate type
                celltype = 'unknown'
            elif len(filtered_nodes) == 1:# Only one candidate type, return directly
                celltype = filtered_nodes[0]
            elif weight_array[0] > 3 and weight_array[1] < 3: # Dominance is too obvious, return directly
                celltype = filtered_nodes[0]
            else:# Dominance is not obvious, perform aggregation
                celltype2celltype_matrix, new_celltypes = self.KG.get_celltype2celltype_kg(filtered_nodes)

                # Set edges corresponding to unselected nodes to negative
                for idx, node in enumerate(new_celltypes): # Iterate through new celltypes
                    if node not in filtered_nodes: # Unselected nodes
                        celltype2celltype_matrix[idx, :] = -abs(celltype2celltype_matrix[idx, :])
                        celltype2celltype_matrix[:, idx] = -abs(celltype2celltype_matrix[:, idx])

                node_attr=[{'weight':0} for _ in new_celltypes]
                sub_g_celltype2celltype = build_graph_from_adjust_matrix(
                    celltype2celltype_matrix,
                    new_celltypes,
                    new_celltypes,
                    'relation_confidence',
                    node_attr, node_attr,mode='n')
                if sub_g_celltype2celltype.number_of_edges() == 0: # These celltypes have no connection, return the result with maximum weight
                    celltype = filtered_nodes[0]
                else: # Has connection
                    # Add weights
                    for node in sub_g_celltype2celltype.nodes:
                        if node in weight_distribution:
                            attrs = weight_distribution[node]
                            sub_g_celltype2celltype.nodes[attrs['name']]['weight'] = attrs["weight"]
                            sub_g_celltype2celltype.nodes[attrs['name']]['initial_weight'] = attrs["weight"] # This node carries initial weight
                        else:
                            sub_g_celltype2celltype.nodes[node]['weight'] = 0
                            sub_g_celltype2celltype.nodes[node]['initial_weight'] = 0
                    sub_g_celltype2celltype = deliver_weight_on_graph_max(sub_g_celltype2celltype, 'weight', 'relation_confidence', alpha=decay_factor)
                    # Update weight information
                    weight_distribution = {node: sub_g_celltype2celltype.nodes[node] for node in sub_g_celltype2celltype.nodes}
                    sorted_nodes = sorted(weight_distribution.items(), key=lambda item: item[1]['weight'], reverse=True)
                    celltype = sorted_nodes[0][0]
            max_initweight_celltype = celltype
        elif mode == 'path':
            # Direct aggregation
            celltype2celltype_matrix, new_celltypes = self.KG.get_celltype2celltype_kg(celltypes)
            # Supplement celltypes that did not appear in edges
            new_celltypes.extend(list(set(celltypes)-set(new_celltypes)))

            # Set edges corresponding to unselected nodes to negative
            for idx, node in enumerate(new_celltypes): # Iterate through new celltypes
                if idx >= celltype2celltype_matrix.shape[0]:
                    break
                if node not in filtered_nodes: # Unselected nodes
                    celltype2celltype_matrix[idx, :] = -abs(celltype2celltype_matrix[idx, :])
                    celltype2celltype_matrix[:, idx] = -abs(celltype2celltype_matrix[:, idx])

            node_attr=[{'weight':0} for _ in new_celltypes]
            sub_g_celltype2celltype = build_graph_from_adjust_matrix(
                celltype2celltype_matrix,
                new_celltypes, new_celltypes,
                'relation_confidence',
                node_attr,node_attr, mode='n')
            # Add weights
            for node in sub_g_celltype2celltype.nodes:
                if node in weight_distribution:
                    attrs = weight_distribution[node]
                    sub_g_celltype2celltype.nodes[attrs['name']]['weight'] = attrs["weight"]
                    sub_g_celltype2celltype.nodes[attrs['name']]['initial_weight'] = attrs["weight"]
                else:
                    sub_g_celltype2celltype.nodes[node]['weight'] = 0
                    sub_g_celltype2celltype.nodes[node]['initial_weight'] = 0
            sub_g_celltype2celltype = deliver_weight_on_graph_max(sub_g_celltype2celltype, 'weight', 'relation_confidence')
            # Update weight information
            weight_distribution = {node: sub_g_celltype2celltype.nodes[node] for node in sub_g_celltype2celltype.nodes}
            sorted_nodes = sorted(weight_distribution.items(), key=lambda item: item[1]['weight'], reverse=True)
            celltype = [sorted_nodes[0][0]]
            max_initweight_celltype = [weight[celltypes.index(celltype[-1])], sorted_nodes[0][0]]
            # Find significant path
            while sub_g_celltype2celltype.in_degree(celltype[-1]) > 0:
                # Get all parent nodes
                parents = [node for node in sub_g_celltype2celltype.predecessors(celltype[-1])]
                # Get weights of all parent nodes
                parent_weights = [weight_distribution[parent]['weight'] for parent in parents]
                weight_array = tempered_softmax_contributions(parent_weights)
                # parent_index = np.argmax(weight_array)
                # celltype.append(parents[parent_index])
                indices = np.where(weight_array > 0.15)[0].tolist()
                if len(indices) == 1: # Only one candidate parent node
                    celltype.append(parents[indices[0]])
                    if weight[celltypes.index(celltype[-1])] > max_initweight_celltype[0]: # Parent node's initial weight is greater than max initial weight
                        max_initweight_celltype = [weight[celltypes.index(celltype[-1])], celltype[-1]]
                else: # Multiple candidate parent nodes
                    multiple_parents = [parents[idx] for idx in indices]
                    # Add all candidate parent nodes
                    celltype.append('|'.join(multiple_parents))
                    init_weights = [weight[celltypes.index(parent)] for parent in multiple_parents]
                    maxindex = np.argmax(init_weights)
                    if init_weights[maxindex] > max_initweight_celltype[0]: # Parent node's initial weight is greater than max initial weight
                        max_initweight_celltype = [init_weights[maxindex], multiple_parents[maxindex]]
                    break
            celltype = '>'.join(celltype)
            max_initweight_celltype = max_initweight_celltype[1]


        celltype_list = [node[0] for node in sorted_nodes]
        weight_list = [node[1]['weight'] for node in sorted_nodes]
        initial_weight_list = [weight[celltypes.index(node[0])] for node in sorted_nodes]

        genecount_KG = []
        if return_genecount_KG:
            for i in celltype_list:
                genecount_KG.append(self.KG.get_genecount_kg(i))

        genecount_graph = []
        if return_genecount_graph:
            graph_count = homolo2celltype_matrix.getnnz(axis=0)
            for i in celltype_list:
                if i not in celltype_nodes:
                    genecount_graph.append(0)
                else:
                    genecount_graph.append(graph_count[celltype_nodes.index(i)])
        
        graph_gene2celltypes = nx.compose(sub_g_gene2celltype, sub_g_celltype2celltype)
        # graph_gene2celltypes = clean_graph_attributes(graph_gene2celltypes)
        # Convert recorder to native type
        signal4recorder = False
        for node in graph_gene2celltypes.nodes:
            if 'recorder' in graph_gene2celltypes.nodes[node]:
                signal4recorder = True
                for key in graph_gene2celltypes.nodes[node]['recorder'].keys():
                    graph_gene2celltypes.nodes[node]['recorder'][key] = float(graph_gene2celltypes.nodes[node]['recorder'][key])
        if not signal4recorder: # No recorder
            logger.warning('No recorder in graph_gene2celltype')
        nx.write_gexf(graph_gene2celltypes, outdir / f'{output_prefix}_gene2celltype.xml')
        return celltype, max_initweight_celltype, celltype_list, weight_list, initial_weight_list, genecount_graph, genecount_KG, homolo2celltype_list, gene2celltype_matrix

    def refine(self, adata:sc.AnnData, group_gene_homolo_weight:pd.DataFrame, 
               celltypes_weight:pd.DataFrame, key_added:str, 
               cluster_key:str='leiden', topk=5, organ:str=None):
        '''
        Perform finer annotation for each cluster.
        Parameters:
            adata (sc.AnnData): Expression data
            group_gene_homolo_weight (pd.DataFrame): Contains cluster (group), gene, and homologous gene (homolo) information
            celltypes_weight (pd.DataFrame): Each cluster's celltypes and their corresponding weights
            key_added (str): New celltype column name
        '''
        data_need_refine = []
        # Convert weights to normal distribution
        celltypes_weight.sort_values(by=['cluster', 'init_weight'], ascending=[True, False], inplace=True)
        for cluster in celltypes_weight['cluster'].unique():
            init_weight = celltypes_weight[celltypes_weight['cluster']==cluster]
            if init_weight.shape[0] < 2: continue
            diff = np.diff(init_weight['init_weight'])
            max_change_idx = np.argmax(abs(diff))
            index = min(max_change_idx, topk)
            if index > 0:
                data_need_refine.append([
                    cluster, init_weight['celltype'][:index+1].tolist()
                ])

        # Annotate each cluster
        for key, celltypes in data_need_refine:
            self.refine_single_cluster(adata, group_gene_homolo_weight, cluster_key, key, 
                                        celltypes, key_added, organ=organ)

    @staticmethod
    def filt_marker_by_moranI(adata:sc.AnnData, moranI_threshold:float=0.5) -> list[str]:
        '''
        Filter var_names in adata using Moran's I.
        Parameters:
            adata: AnnData object
            moranI_threshold: Moran's I threshold
        Returns:
            marker_list: Filtered marker list
        '''
        if hasattr(adata.X, 'toarray'):
            adata.X = adata.X.toarray()
        marker_df = pd.DataFrame(adata.X)
        marker_df.columns = adata.var_names.tolist()
        marker_df_T = marker_df.T
        marker_df_T['moransI'] = sc.metrics.morans_i(adata, vals=marker_df_T)
        qualified_marker = marker_df_T.loc[(marker_df_T['moransI'] >= moranI_threshold)].index
        marker_df = marker_df[qualified_marker]
        return marker_df.columns.tolist()

    def refine_single_cluster(self, adata:sc.AnnData, group_gene_homolo_weight:pd.DataFrame, 
            cluster_key:str, cluster_id:str, candidate_celltype:list[str], 
            key_added:str, organ:str=None, moranI_threshold=0.5, 
            split_method:Literal['bindiv','argmax']='bindiv', 
            markergene_method:Literal['diff','all']='diff', celltype_geneCount_gene=None):
        '''
        Annotate a single cluster.
        Parameters:
            adata: AnnData object
            markers: Contains cluster (group), gene, homologous gene (homolo), and weight information
            cluster_key: Column name in adata.obs where clustering results are stored
            cluster_id: Cluster ID
            candidate_celltype: Candidate cell types
            key_added: New column name
            moranI_threshold: Filtering threshold, [-1, 1]
        Returns:
            list[tuple[str, int, list[str]]]: Cell type, differential gene count, differential gene list
        Side effects:
            Adds a new column to adata.obs named key_added
            Determines gene composition for each cell type
            Determines the order of cell type partitioning
        '''
        logger.info('refine_single_cluster: %s, %s', cluster_id, cluster_key)
        if 'pca' not in adata.uns.keys():
            logger.info('run sc.pp.pca(adata)')
            sc.pp.pca(adata)
        if 'connectivities' not in adata.obsp.keys():
            logger.info('run sc.pp.neighbors(adata)')
            sc.pp.neighbors(adata)

        # Validate organ
        organ = self.KG.check_organ(organ)
        if organ is None:
            logger.warning("organ is unavailable, set to None.")

        if celltype_geneCount_gene is None:
            # Get mapping from homologous gene to marker
            homolo2gene = {}
            group_gene_homolo_weight.apply(lambda x: homolo2gene.update({x['homolo']: x['gene']}), axis=1)
            # Look up marker genes in knowledge graph for each cell type
            source, target, matrix = self.KG.get_gene2celltype_kg(
                homolo_nodes=list(homolo2gene.keys()),
                organ=organ, 
                candidate_type=candidate_celltype)
            celltype2markers = dict()
            X_coo = matrix.tocoo()
            for i, j, v in zip(X_coo.row, X_coo.col, X_coo.data):
                if v == 0: continue
                if target[j] in celltype2markers.keys():
                    celltype2markers[target[j]].append(homolo2gene[source[i]])
                else:
                    celltype2markers[target[j]] = [homolo2gene[source[i]]]
            # Remove duplicates
            if len(celltype2markers) == 0:
                raise Exception('no celltype2markers!')
            celltype2markers = {k:list(set(v)) for k, v in celltype2markers.items()}
            # Calculate Moran's I to filter available markers, determine whether to continue based on remaining markers for each cell type
            all_markers = set()
            for markers in celltype2markers.values():
                all_markers.update(markers)

            adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)
            adata_sub = adata[adata.obs[cluster_key] == str(cluster_id), list(all_markers)].copy()
            marker_available = all_markers
            if moranI_threshold>-1 and moranI_threshold<1:
                marker_available = Xener.filt_marker_by_moranI(adata_sub, moranI_threshold)
                if len(marker_available) == 0:
                    marker_available = all_markers
                    logger.info('no marker filtered by moranI, use all markers')
            else:
                logger.info('skip moranI filter')
            # all_markers = set(marker_available)
            # Filter markers
            for key in list(celltype2markers.keys()):
                celltype2markers[key] = [marker for marker in celltype2markers[key] if marker in marker_available]

            '''
            Determine cell type selection order
            '''
            # Convert to sets
            celltype2markers_set = {key: set(celltype2markers[key]) for key in celltype2markers.keys()}
            celltype_geneCount_gene = []
            all_gene_set = set()
            for key, markers_set in celltype2markers_set.items():
                if len(markers_set) == 0:
                    logger.info('%s deleted from geneCount gene', key)
                    continue
                all_gene_set |= markers_set
                celltype_geneCount_gene.append((key, len(markers_set), list(markers_set)))
        else:
            all_gene_set = set()
            for key, _, markers_set in celltype_geneCount_gene:
                all_gene_set |= set(markers_set)
            adata_sub = adata[adata.obs[cluster_key] == str(cluster_id), list(all_gene_set)].copy()
            logger.info('use celltype_geneCount_gene')
        
        celltype_diffgeneCount_gene = []
        all_difference_set = set()
        for key, _, markers_set in celltype_geneCount_gene:
            difference_set = set(markers_set)
            for key2, _, markers_set2 in celltype_geneCount_gene:
                if key != key2:
                    difference_set -= set(markers_set2)
            if len(difference_set) > 0:
                # Use difference set length
                all_difference_set |= difference_set
                celltype_diffgeneCount_gene.append((key, len(difference_set), list(difference_set)))
            else:
                logger.info('%s deleted from diff gene', key)
        celltype_diffgeneCount_gene = sorted(celltype_diffgeneCount_gene, key=lambda item: item[1], reverse=True)

        adata_sub.obs[key_added] = 'waitting'

        if markergene_method == 'diff':
            celltype_geneCount_gene4reine = celltype_diffgeneCount_gene.copy()
            gene_set4refine = all_difference_set.copy()
        elif markergene_method == 'all':
            celltype_geneCount_gene4reine = celltype_geneCount_gene.copy()
            gene_set4refine = all_gene_set.copy()

        # Start refinement
        if split_method == 'bindiv':
            adata_sub.obs['need_refine'] = True# Not yet refined marker
            adata_sub.obs['score'] = 0
            for idx, (goal_type, goal_gene_num, goal_gene) in enumerate(celltype_geneCount_gene4reine):
                to_refine = adata_sub[adata_sub.obs['need_refine'],:].copy()
                gene_set4refine -= set(goal_gene)
                # to_refine.obs[] = to_refine[:,goal_gene].X.mean(axis=1)
                if len(gene_set4refine) == 0:
                    to_refine.obs['score'] = 1
                else:
                    # Use weights
                    score =  to_refine[:,goal_gene].X.mean(axis=1) \
                        - to_refine[:,list(gene_set4refine)].X.mean(axis=1)
                    score = to_refine.obsp['connectivities'] @ score / to_refine.obsp['connectivities'].sum(axis=1)
                    to_refine.obs['score'] = score

                mean_exp = adata[:,goal_gene].X.mean(axis=1)
                adata.obs[f'{idx}_{goal_type}_EXP'] = mean_exp
                if len(gene_set4refine) > 0:
                    mean_exp = adata[:,list(gene_set4refine)].X.mean(axis=1)
                    adata.obs[f'{idx}_!{goal_type}_EXP'] = mean_exp
                    adata.obs[f'{idx}_diff_EXP'] = adata.obs[f'{idx}_{goal_type}_EXP'] - adata.obs[f'{idx}_!{goal_type}_EXP']
                # Annotate on adata_sub
                sub_type_count = to_refine[to_refine.obs['score']>0].shape[0]
                adata_sub.obs.loc[
                    to_refine[to_refine.obs['score']>0].obs_names, key_added] = goal_type
                adata_sub.obs.loc[
                    to_refine[to_refine.obs['score']>0].obs_names, 'need_refine'] = False
                logger.info('refine_single_cluster, cluster_id[%s], goal_type[%s], goal_gene_num[%s], sub_type_count[%s]', cluster_id, goal_type, goal_gene_num, sub_type_count)

        elif split_method == 'argmax':
            exps=[]
            for goal_type, goal_gene_num, goal_gene in celltype_geneCount_gene4reine:
                # to_refine = adata_sub[adata_sub.obs['need_refine'],:].copy()
                mean_exp = adata[:,goal_gene].X.mean(axis=1)
                adata.obs[f'{goal_type}_EXP'] = mean_exp
                exps.append(adata_sub[:,goal_gene].X.mean(axis=1))
            exps = np.array(exps)
            # Remove extra dimensions
            exps = np.squeeze(exps)
            exps = adata_sub.obsp['connectivities'] @ exps.T / adata_sub.obsp['connectivities'].sum(axis=1)
            type_id = np.argmax(exps, axis=1)
            adata_sub.obs[key_added] = [celltype_geneCount_gene4reine[i][0] for i in type_id]

        logger.info('refine_single_cluster %s total: %s', cluster_id, adata_sub.obs[key_added].unique().tolist())
        adata.obs.loc[adata_sub.obs_names, key_added] = adata_sub.obs[key_added]

        return celltype_geneCount_gene, celltype_diffgeneCount_gene