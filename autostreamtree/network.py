"""
River network functions of autoStreamTree: reading a network, finding paths
between sites along the river, extracting the subgraph they use, and the
pairwise river distances and incidence matrix. They need only the network
libraries (no VCF reading), so other packages (e.g. ResistNet) can use
them; autostreamtree.functions re-exports them.
"""
import itertools
import math
import pickle
import sys
import traceback

import momepy
import networkx as nx
import numpy as np
import pandas as pd
import pyogrio
from networkx import NodeNotFound


def prune_graph(G, edge_list, reachid_col):
    """
    Prunes a graph to only retain edges whose 'reachid_col' matches values in
    'edge_list'.

    Args:
        G: The input NetworkX Graph.
        edge_list: A list of values to filter edges.
        reachid_col: The edge attribute to be checked against 'edge_list'.

    Returns:
        A pruned NetworkX Graph.
    """
    for u, v in G.edges():
        _ = G.get_edge_data(u, v)
        # print(f"Edge ({u}, {v}): {data}")
    # Get edges to be retained using get_edge_data method
    edges_to_keep = [(u, v) for u, v in G.edges() if
                     G.get_edge_data(u, v).get(reachid_col) in edge_list]

    # Check if there are no edges to keep
    if not edges_to_keep:
        raise ValueError(
            "There are no edges to retain based on the given edge",
            "list and attribute."
        )

    # Create a new graph with only the edges to be retained
    pruned_G = G.edge_subgraph(edges_to_keep).copy()

    # Remove isolated nodes
    isolated_nodes = list(nx.isolates(pruned_G))
    pruned_G.remove_nodes_from(isolated_nodes)

    return pruned_G


def _reach_ids(edge_list, reachid_col):
    """
    The reach IDs of an edge list (a DataFrame with reachid_col, or any
    iterable of IDs) as a list, without missing values; None if no edge
    list.
    """
    if edge_list is None:
        return None
    if isinstance(edge_list, pd.DataFrame):
        if reachid_col not in edge_list.columns:
            raise ValueError(f"{reachid_col} not found in the edge list")
        return edge_list[reachid_col].dropna().tolist()
    return pd.Series(list(edge_list)).dropna().tolist()


def _where_in(col, ids):
    """An OGR SQL filter selecting the features whose col is in ids."""
    if all(isinstance(x, (int, float, np.integer, np.floating))
           and not isinstance(x, bool) for x in ids):
        values = ",".join(str(x) for x in ids)
    else:
        values = ",".join("'" + str(x).replace("'", "''") + "'"
                          for x in ids)
    return f"{col} IN ({values})"


def read_network(network, shapefile, reachid_col="HYRIV_ID",
                 edge_list=None):
    """
    Reads a network from a saved file or builds a network from a shapefile.

    Args:
        network: Path to the saved network file (pickle format). If provided,
                  the function will read the network from this file.
        shapefile: Path to the shapefile (or other vector file) to build the
                   network from. This is used if the `network` parameter is
                   not provided.
        reachid_col: Attribute giving each reach's ID.
        edge_list: Optional reaches to keep (a DataFrame with reachid_col, or
                   an iterable of IDs). From a shapefile, only these reaches
                   are read; a saved network is pruned to them (see
                   prune_graph()).

    Returns:
        A NetworkX Graph object representing the network.
    """
    ids = _reach_ids(edge_list, reachid_col)

    # Check if a saved network file is provided
    if network:
        print("Reading network from saved file: ", network)
        # Read the network from the saved file and convert it to an undirected
        # graph
        with open(network, 'rb') as f:
            G = nx.Graph(pickle.load(f)).to_undirected()
        if ids is not None:
            G = prune_graph(G, ids, reachid_col)
        return G

    # If no saved network file is provided, build the network from the
    # shapefile
    print("Building network from shapefile:", shapefile)
    print("WARNING: This can take a while with very large files!")
    if ids is None:
        rivers = pyogrio.read_dataframe(shapefile)
    elif not ids:
        raise ValueError("The edge list has no reach IDs")
    else:
        # Read only the listed reaches; if the filter fails or matches
        # nothing (e.g. IDs stored as text), read all and filter here
        rivers = None
        try:
            rivers = pyogrio.read_dataframe(
                shapefile, where=_where_in(reachid_col, ids))
        except Exception as e:
            print(f"WARNING: filtering reaches while reading failed ({e}); "
                  "reading all of them instead")
        if rivers is None or rivers.empty:
            rivers = pyogrio.read_dataframe(shapefile)
            if reachid_col not in rivers.columns:
                raise ValueError(
                    f"{reachid_col} not found in the shapefile attributes")
            keep = pd.Series(ids)
            try:
                keep = keep.astype(rivers[reachid_col].dtype)
            except (TypeError, ValueError):
                rivers[reachid_col] = rivers[reachid_col].astype(str)
                keep = keep.astype(str)
            rivers = rivers[rivers[reachid_col].isin(keep)]
        if rivers.empty:
            raise ValueError(
                "None of the reaches in the edge list are in the shapefile")
    # Convert the GeoDataFrame to a NetworkX Graph object; multi-part
    # reaches are split into their parts (momepy needs single lines)
    G = momepy.gdf_to_nx(rivers.explode(index_parts=False),
                         approach="primal", directed=False, multigraph=False)
    return G


def river_length(len_col):
    """
    Edge weight for shortest paths along the river: each edge's length
    (len_col). On a network without loops (e.g. HydroRIVERS) the path
    between two nodes is unique; with loops (braided channels, canals), the
    shortest path by river length is used. (The weight was the inverse of
    the length, which picked the route through the longest edges.)

    Args:
        len_col (str): Attribute giving each edge's length.

    Returns:
        callable: The weight function, for networkx.
    """
    def weight(u, v, attributes):
        return float(attributes[len_col])
    return weight


def get_stream_mats(points, graph, len_col):
    """
    Computes pairwise stream distances and 0/1 incidence matrix for StreamTree
    calculations.

    Args:
        points (dict): Dictionary of point indices and their corresponding node
                       IDs in the graph.
        graph (networkx.Graph): NetworkX graph object representing the stream
                                network.
        len_col (str): Attribute name for the length of the edges in the graph.

    Returns:
        tuple: Pair of numpy.ndarray representing the pairwise stream distance
               matrix and incidence matrix.
    """
    dist = np.zeros((len(points), len(points)))
    inc = np.zeros((nCr(len(points), 2), len(graph.edges())), dtype=int)

    # Establish as NaN
    dist[:] = np.nan

    index = 0
    for ia, ib in itertools.combinations(range(0, len(points)), 2):
        path = nx.bidirectional_dijkstra(graph, points.values()[ia],
                                         points.values()[ib],
                                         weight=river_length(len_col))
        if path:
            dist[ia, ib] = float(sum(path_edge_attributes(graph, path[1],
                                                          len_col)))
            dist[ib, ia] = dist[ia, ib]

        # Incidence matrix: assign 1 if edge is in the path, 0 otherwise
        for ie, edge in enumerate(graph.edges()):
            if find_pair(path[1], edge[0], edge[1]):
                inc[index, ie] = 1
            else:
                inc[index, ie] = 0
        index += 1

    np.fill_diagonal(dist, 0.0)
    return dist, inc


def find_pair(lst, x, y):
    """
    Check if two elements are consecutive in a list, irrespective of their
    order.

    Args:
        lst (list): The list to search for the pair.
        x (Any): The first element of the pair.
        y (Any): The second element of the pair.

    Returns:
        bool: True if the elements are consecutive in the list, False otherwise
    """
    if x not in lst or y not in lst:
        return False
    elif abs(lst.index(x) - lst.index(y)) == 1:
        return True
    else:
        return False


def nCr(n, k):
    """
    Calculate the number of combinations, n choose k.

    Args:
        n (int): The number of elements.
        k (int): The number of elements to choose.

    Returns:
        int: The number of possible combinations.
    """
    f = math.factorial
    return f(n) // f(k) // f(n - k)


def path_edge_attributes(graph, path, attribute):
    """
    Get the attribute values for edges in a given path.

    Args:
        graph (NetworkX Graph): The graph containing the edges.
        path (list): The list of nodes forming the path.
        attribute (str): The edge attribute to get the values for.

    Returns:
        list: A list of attribute values for the edges in the path.
    """
    return [graph[u][v][attribute] for (u, v) in zip(path, path[1:])]


def path_subgraph(graph, nodes, method, id_col, len_col):
    """
    Find and extract paths between points from a graph.

    Args:
        graph (NetworkX Graph): The input graph.
        nodes (dict): A dictionary of nodes to extract paths between.
        method (callable): The method to build the subgraph.
        id_col (str): The column name for node ID.
        len_col (str): The column name for edge length.

    Returns:
        NetworkX Graph: A subgraph containing the extracted paths.
    """
    k = nx.Graph()

    p1 = list(nodes.values())[0]
    for p2 in list(nodes.values())[1:]:
        try:
            # Find the shortest path between the two points
            path = nx.bidirectional_dijkstra(graph, p1, p2,
                                             weight=river_length(len_col))

            # Traverse the nodes in the path to build a minimal set of edges
            method(k, graph, nodes.values(), id_col, len_col, path[1])

            if p1 not in k:
                k.add_node(p1)
            if p2 not in k:
                k.add_node(p2)

        except NodeNotFound as e:
            print("Node not found:", e)
        except Exception as e:
            traceback.print_exc()
            print("Something unexpected happened:", e)
            sys.exit(1)
    return k


def extract_full_subgraph(subgraph, graph, nodelist, id_col, len_col, path):
    """
    Extracts the full subgraph from the given nodes.

    Args:
        subgraph (NetworkX Graph): The subgraph to be modified.
        graph (NetworkX Graph): The input graph.
        nodelist (list): The list of nodes.
        id_col (str): The column name for node ID.
        len_col (str): The column name for edge length.
        path (list): The path between nodes.
    """
    # Iterate through the nodes in the path
    for first, second in zip(path, path[1:]):
        # Add nodes to the subgraph if they are not already present
        if first not in subgraph:
            subgraph.add_node(first)
        if second not in subgraph:
            subgraph.add_node(second)

        # Add the edge between the nodes with the corresponding edge data
        dat = graph.get_edge_data(first, second)
        subgraph.add_edge(first, second, **dat)


def extract_minimal_subgraph(subgraph, graph, nodelist, id_col, len_col, path):
    """
    Extracts a simplified subgraph from paths, keeping only terminal and
    junction nodes.

    Args:
        subgraph (NetworkX Graph): The subgraph to be modified.
        graph (NetworkX Graph): The input graph.
        nodelist (list): The list of nodes.
        id_col (str): The column name for edge ID.
        len_col (str): The column name for edge length.
        path (list): The path between nodes.
    """
    curr_edge = {id_col: list(), len_col: 0.0}
    curr_start = None

    # Iterate through each pair of nodes in the path
    for first, second in zip(path, path[1:]):
        if not curr_start:
            curr_start = first
            if first in nodelist or len(graph[first]) > 2:
                subgraph.add_node(first)

        # Add path attributes to current edge
        dat = graph.get_edge_data(first, second)
        curr_edge[id_col].extend([dat[id_col]] if not
                                 isinstance(dat[id_col], list) else
                                 dat[id_col])
        curr_edge[len_col] = float(curr_edge[len_col]) + float(dat[len_col])

        # If the second node is a STOP node (in nodelist or is a junction)
        if second in nodelist or len(graph[second]) > 2:
            # Add node to subgraph
            subgraph.add_node(second)
            # Link current attribute data
            subgraph.add_edge(curr_start, second, **curr_edge)
            # Empty edge attributes and set current second to curr_start
            curr_edge = {id_col: list(), len_col: 0}
            curr_start = second
        else:
            # Otherwise, continue building the current edge
            continue


def snap_to_node(graph, pos):
    """
    Finds the closest node to the given [x, y] coordinates in the graph.

    Args:
        graph (NetworkX Graph): The input graph.
        pos (tuple): A tuple of [x, y] coordinates.

    Returns:
        tuple: The closest node to the input coordinates.
    """
    nodes = np.array(graph.nodes())
    node_pos = np.argmin(np.sum((nodes - pos) ** 2, axis=1))
    return tuple(nodes[node_pos])
