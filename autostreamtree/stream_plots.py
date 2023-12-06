import os
import sys

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from typing import List

def plot_edges_to_streams(network: str, res: pd.DataFrame, oname: str) -> None:
    """
    Plot a stream network with resistance values colored by a gradient.

    Args:
        network (str): The filename of a NetworkX graph object in .gpickle format.
        res (pandas.DataFrame): A DataFrame with an "EDGE_ID" column and a "Resistance" column.
        oname (str): The filename to save the output plot as.
    """
    G = nx.Graph(nx.read_gpickle(network).to_undirected())
    geoDF = momepy.nx_to_Gdf(G, lines=True)
    geoDF = geoDF.merge(res, on="EDGE_ID")
    sns.set(style="ticks")
    geoDF.plot(column="Resistance", cmap="RdYlGn_r", legend=True)
    plt.title("Stream network colored by resistance")
    plt.savefig(str(oname) + ".streamsByResistance.pdf")
    plt.clf()
    plt.close()
