#!/usr/bin/env python
import sys 
import os

import argparse
import networkx as nx
import numpy as np
import geopandas as gpd
import momepy
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt

# Add the parent directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

import autostreamtree.functions as ast

def main():

    # read arguments 
    parser = argparse.ArgumentParser(description='Plot a stream tree with fittedD attribute as dendrogram')
    parser.add_argument('--shp', type=str, required=True, help='Network as shapefile')
    parser.add_argument("--points", type=str, required=True, help="Sample coordinates")
    parser.add_argument("--pop", type=str, required=False, default=None, help="Population map file (optional)")
    parser.add_argument('--edge_id', type=str, default='EDGE_ID', help='Edge ID attribute.')
    parser.add_argument('--dist', type=str, default='fittedD', help='Dist attribute.')
    parser.add_argument('--out', type=str, default='out', help='Output prefix')
    args = parser.parse_args()
    setattr(args, "geopop", False)
    setattr(args, "clusterpop", False)
    setattr(args, "run", "STREAMTREE")

    params = args

    shp = gpd.read_file(params.shp)
    # Convert the GeoDataFrame to a NetworkX Graph object
    G = momepy.gdf_to_nx(shp, approach="primal", directed=False, multigraph=False)
    #print_edge_attributes(G)

    # read point coordinate data
    points = pd.read_csv(params.points, sep="\t", header=0)
    (point_coords, pop_coords, popmap) = ast.process_samples(params, points, G)


    K = ast.path_subgraph(G, pop_coords, ast.extract_minimal_existing, params.edge_id, params.dist)

    plot_network(K, params.edge_id, point_coords, params.out)

    # convert to newick
    if params.pop is not None:
        node_names = {value: key for key, value in pop_coords.items()}
    else:
        node_names = {value: key for key, value in point_coords.items()}

    root = list(K.nodes)[0]
    newick_str = to_newick(root, None, K, params.dist, node_names)
    print(newick_str)

# def plot_newick(newick_string):
#     # Create a toytree object
#     tre = toytree.tree(newick_string)
    
#     # Draw the tree
#     tre.draw(tip_labels_align=True)

def to_newick(node, parent_node, graph, dist_attr, node_names=None):
    # Get the neighbors of the current node
    neighbors = list(graph.neighbors(node))

    # If the current node has a parent node, remove it from the neighbors list
    if parent_node in neighbors:
        neighbors.remove(parent_node)

    # If there are no more neighbors, this is a leaf node
    if not neighbors:
        # Use name from dict if available
        node_name = node_names.get(node, node) if node_names else node
        return str(node_name) + ":" + str(graph.edges[node, parent_node][dist_attr])

    # If there are neighbors, this is an internal node
    # We need to build a subtree for each neighbor
    subtrees = [to_newick(neighbor, node, graph, dist_attr, node_names) for neighbor in neighbors]

    # Join the subtrees with commas and enclose them in parentheses
    return "(" + ",".join(subtrees) + "):" + str(graph.edges[node, parent_node][dist_attr] if parent_node is not None else 0.0)


def plot_network(G, edge_id_attr, points, oname):
    # Grab real coordinates as node positions for plotting
    pos = {n: n for n in G.nodes}

    # Make a color map to color sample points and junctions differently
    color_map = ["blue" if node in points.values() else "black" for node in G]

    nx.draw_networkx(G, pos, with_labels=False, node_color=color_map, node_size=10)
    #nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=6)

    network_plot = str(oname) + ".subGraph.pdf"
    plt.savefig(network_plot)

def print_edge_attributes(network):
    for u, v, attrs in network.edges(data=True):
        print(f"Edge ({u}, {v}): {attrs}")

if __name__ == "__main__":
    main()