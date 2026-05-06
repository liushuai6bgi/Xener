import os
import json
from typing import Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp
import networkx as nx

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

def filter_top_k_col_per_row(adj_matrix:sp.csr_matrix, k:int, col:list) -> tuple[sp.csr_matrix, list]:
    """
    Filter the top-k highest-weight non-zero columns per row.

    Args:
        adj_matrix: Sparse adjacency matrix.
        k: Number of top-weight columns to retain per row.
        col: List of column names.

    Returns:
        Filtered sparse matrix and corresponding column name list.
    """
    adj_matrix = adj_matrix.tocsr()
    num_row, num_col = adj_matrix.shape
    result_matrix = sp.lil_matrix((num_row, num_col))

    for i in range(num_row):
        gene_indices = adj_matrix[i, :].nonzero()[1]
        weights = adj_matrix[i, gene_indices].toarray().flatten()
        top_k_indices = np.argsort(weights)[-k:]
        for idx in top_k_indices:
            result_matrix[i, gene_indices[idx]] = adj_matrix[i, gene_indices[idx]]
    result_matrix = result_matrix.tocsr()
    col_sums = result_matrix.sum(axis=0)
    non_empty_cols = np.where(col_sums > 0)[1]
    return result_matrix[:, non_empty_cols], [ col[idx] for idx in non_empty_cols.tolist() ]

def joint_matrix(matrix_NxN:sp.csr_matrix, matrix_NxM:sp.csr_matrix, matrix_MxM:sp.csr_matrix):
    """
    Joint embedding: combine three matrices into a single sparse matrix.

    Args:
        matrix_NxN: Top-left matrix (N x N).
        matrix_NxM: Top-right matrix (N x M).
        matrix_MxM: Bottom-right matrix (M x M).

    Returns:
        Combined sparse matrix ((N+M) x (N+M)).
    """
    return sp.vstack([
                sp.hstack([matrix_NxN, matrix_NxM]),
                sp.hstack([matrix_NxM.T, matrix_MxM])
            ]).tocsr()

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