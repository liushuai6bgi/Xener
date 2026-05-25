import os
from typing import Literal

import numpy as np
import scipy.sparse as sp
import networkx as nx
from .logger import logger

def name2path(path:str, postfix:list[str]) -> dict[str, str]:
    """
    Collect all file names and their paths under a directory, returning a dictionary.

    Args:
        path: Directory path to traverse.
        postfix: List of file suffixes to filter target files.

    Returns:
        Dictionary mapping file names to file paths.
    """
    name2path = {}
    for subdir, _, files in os.walk(path):
        for file in files:
            if file.split('.')[-1] in postfix:
                file_name = file.split('.')[0]
                if file_name in name2path.keys():
                    raise Exception(f"重复的名称:{file_name}")
                name2path[file_name] = os.path.join(subdir, file)
    return name2path

def deliver_weight_on_graph_sum(graph:nx.DiGraph, node_key:str, edge_key:str, weight_recorder_key:str='recorder', neg_edge:bool=False, alpha:float=0.9) -> nx.DiGraph:
    """
    Propagate node weights on a directed weighted NetworkX graph by summing all child node weights (duplicate nodes retain the maximum).

    Args:
        graph: Directed graph; must be acyclic.
        node_key: Key to store node weights.
        edge_key: Key to store edge weights.
        weight_recorder_key: Key to track the composition of node weights during propagation, avoiding duplicate weight addition for the same child node.
        neg_edge: Whether to propagate negative edge weights to child nodes.
        alpha: Weight decay factor.

    Returns:
        The directed graph after weight propagation.
    """
    assert nx.is_directed_acyclic_graph(graph)
    in_degrees = dict(graph.in_degree)
    start_nodes = [node for node, degree in in_degrees.items() if degree == 0]
    while len(start_nodes) > 0:
        for start_node in start_nodes:
            del in_degrees[start_node]
            if weight_recorder_key not in graph.nodes[start_node].keys():
                start_weight_recorder = {}
            else:
                start_weight_recorder = graph.nodes[start_node][weight_recorder_key]

            for next_node in graph.successors(start_node):
                start_node_weight = graph.nodes[start_node][node_key]
                if weight_recorder_key not in graph.nodes[next_node].keys():
                    next_weight_recorder = {}
                else:
                    next_weight_recorder = graph.nodes[next_node][weight_recorder_key]
                relation_confidence = float(graph[start_node][next_node][edge_key])
                if not neg_edge and relation_confidence < 0:
                    continue
                graph.nodes[next_node][node_key] += graph.nodes[start_node][node_key] * relation_confidence * alpha
                for node, weight in start_weight_recorder.items():
                    start_node_weight -= weight
                    if node not in next_weight_recorder.keys():
                        next_weight_recorder[node] = weight
                        continue
                    new_weight = float(weight) * relation_confidence * alpha
                    old_weight = next_weight_recorder[node]
                    next_weight_recorder[node] = max(new_weight, old_weight)
                    graph.nodes[next_node][node_key] -= min(new_weight, old_weight)
                in_degrees[next_node] -= 1
                next_weight_recorder[start_node] = start_node_weight * relation_confidence * alpha
                graph.nodes[next_node][weight_recorder_key] = next_weight_recorder
        start_nodes = [node for node, degree in in_degrees.items() if degree == 0]
    return graph

def deliver_weight_on_graph_max(graph:nx.DiGraph, node_key:str, edge_key:str, weight_recorder_key:str='recorder', neg_edge:bool=False, alpha:float=0.9) -> nx.DiGraph:
        '''
        Propagate node weights on a directed weighted NetworkX graph, retaining only the highest-weight path. The weight_recorder_key stores direct child nodes.

        Args:
            graph (networkx.DiGraph): A directed graph; must be acyclic.
            node_key (str): Key to store node weights.
            edge_key (str): Key to store edge weights.
            weight_recorder_key (str): Key to track the composition of node weights during propagation, avoiding duplicate weight addition for the same child node.
            neg_edge (bool): Whether to propagate negative edge weights to child nodes.
            alpha (float): Weight decay factor.
        '''
        # Check for cycles
        assert nx.is_directed_acyclic_graph(graph)
        # Get nodes with zero in-degree, starting points for weight propagation
        in_degrees = dict(graph.in_degree) # Track in-degrees; use this dict to find starting points
        start_nodes = [node for node, degree in in_degrees.items() if degree == 0]
        while len(start_nodes) > 0:
            for start_node in start_nodes:
                # Current start node has finished propagation; remove it
                del in_degrees[start_node]
                # Accept the highest-weight child node path
                if weight_recorder_key in graph.nodes[start_node].keys() and len(graph.nodes[start_node][weight_recorder_key]) > 0:
                    graph.nodes[start_node][node_key] += max(graph.nodes[start_node][weight_recorder_key].values())
                # Propagate weight to parent node
                for next_node in graph.successors(start_node):
                    relation_confidence = float(graph[start_node][next_node][edge_key])
                    if not neg_edge and relation_confidence < 0:
                        continue
                    if weight_recorder_key not in graph.nodes[next_node].keys():
                        graph.nodes[next_node][weight_recorder_key] = {}
                    graph.nodes[next_node][weight_recorder_key][start_node] = graph.nodes[start_node][node_key] * relation_confidence * alpha
                    # Propagation complete; decrement in-degree
                    in_degrees[next_node] -= 1
            # Get new starting points: nodes with zero in-degree
            start_nodes = [node for node, degree in in_degrees.items() if degree == 0]
        return graph

def build_graph_from_adjust_matrix(matrix:sp.csr_matrix, source_names:list[str],
                                   target_names:list[str], edge_weight_key:str, source_attr:list[dict]=None,
                                   target_attr:list[dict]=None, mode:Literal['e','n']='e'):
    '''
    Build a graph from an adjacency matrix by adding edges. The resulting graph may not contain all nodes.

    Args:
        matrix (sp.csr_matrix): Adjacency matrix storing edge weights.
        source_names (list[str]): Source node names; count should be >= number of non-zero rows in matrix.
        target_names (list[str]): Target node names; count should be >= number of non-zero columns in matrix.
        edge_weight_key (str): Edge weight attribute name.
        source_attr (list[dict]): Source node attributes, defaults to empty.
        target_attr (list[float]|float): Target node attributes, defaults to empty.
        mode (Literal['e','n']): Node addition mode; defaults to 'e' (edge-only). If 'n', all nodes are added.
    Returns:
        nx.DiGraph: Directed weighted graph.
    '''
    graph = nx.DiGraph()
    # Add nodes

    X_coo = matrix.tocoo()
    for col, raw, v in zip(X_coo.row, X_coo.col, X_coo.data):
        if v == 0:  # Skip zero-weight edges
            continue
        if not graph.has_node(source_names[col]):
            attr = source_attr[col] if isinstance(source_attr,list) else {}
            graph.add_node(source_names[col], **attr)
        if not graph.has_node(target_names[raw]):
            attr = target_attr[raw] if isinstance(target_attr,list) else {}
            graph.add_node(target_names[raw], **attr)

        graph.add_edge(source_names[col], target_names[raw], **{edge_weight_key:float(v)})
    if mode == 'n':
        for node in source_names:
            if not graph.has_node(node):  # Add missing nodes
                col = source_names.index(node)
                attr = source_attr[col] if isinstance(source_attr,list) else {}
                graph.add_node(node, **attr)
        for node in target_names:
            if not graph.has_node(node):  # Add missing nodes
                raw = target_names.index(node)
                attr = target_attr[raw] if isinstance(target_attr,list) else {}
                graph.add_node(node, **attr)
    return graph

def tempered_softmax_contributions(weights, temperature=1.0):
    """
    Softmax function with a temperature parameter.
    temperature < 1: Enhances the contribution of the maximum value.
    temperature > 1: Smooths contribution differences.
    temperature = 1: Standard softmax.
    """
    weights = np.array(weights)
    weights_shifted = weights - np.max(weights)
    exp_weights = np.exp(weights_shifted / temperature)
    return exp_weights / np.sum(exp_weights)

def keep_max(sp_m:sp.csr_matrix, axis:int) -> sp.csr_matrix:
    """
    Keep only the maximum value(s) along the specified axis.
    
    Parameters
    ----------
    axis : int
        - axis=1 : Keep only the maximum value in each row.
        - axis=0 : Keep only the maximum value in each column.
        - else : Keep only the maximum value globally.
    
    Shape is preserved. Returns a new sparse matrix.
    """
    C = sp_m.tocoo()

    def _flatten_max(max_result):
        if hasattr(max_result, 'toarray'):
            max_result = max_result.toarray()
        return np.asarray(max_result).flatten()

    if axis == 1:
        max_vals = _flatten_max(sp_m.max(axis=1))
        keep = C.data == max_vals[C.row]
    elif axis == 0:
        max_vals = _flatten_max(sp_m.max(axis=0))
        keep = C.data == max_vals[C.col]
    else:
        global_max = sp_m.max()
        keep = C.data == global_max
    
    result = sp.csr_matrix(
        (C.data[keep], (C.row[keep], C.col[keep])),
        shape=sp_m.shape
    )
    return result

def remove_zeros(sp_m: sp.spmatrix, axis: int | None = None) -> tuple[sp.csr_matrix, list[int], list[int]]:
    """
    Remove zero rows and/or columns from a sparse matrix.
    
    Parameters
    ----------
    axis : int or None
        - 0 : remove zero rows only
        - 1 : remove zero columns only
        - None : remove both zero rows and zero columns
                 (rows and cols are detected on the original matrix, then removed together)
    
    Returns
    -------
    If axis=0: (clean_matrix, removed_rows)
    If axis=1: (clean_matrix, removed_cols)
    If axis=None: (clean_matrix, removed_rows, removed_cols)
    """
    sp_m = sp_m.tocsr()
    removed_rows, removed_cols = [], []
    
    if axis == 0:
        row_nnz = np.diff(sp_m.indptr)
        keep_rows = np.where(row_nnz > 0)[0]
        removed_rows = np.where(row_nnz == 0)[0].tolist()
        clean = sp_m[keep_rows, :]
    
    elif axis == 1:
        csc = sp_m.tocsc()
        col_nnz = np.diff(csc.indptr)
        keep_cols = np.where(col_nnz > 0)[0]
        removed_cols = np.where(col_nnz == 0)[0].tolist()
        clean = csc[:, keep_cols]
    
    else:
        row_nnz = np.diff(sp_m.indptr)
        keep_rows = np.where(row_nnz > 0)[0]
        removed_rows = np.where(row_nnz == 0)[0].tolist()
        
        csc = sp_m.tocsc()
        col_nnz = np.diff(csc.indptr)
        keep_cols = np.where(col_nnz > 0)[0]
        removed_cols = np.where(col_nnz == 0)[0].tolist()
        
        clean = sp_m[keep_rows, :][:, keep_cols]
    return clean, removed_rows, removed_cols
