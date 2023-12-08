import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from shapely.geometry import MultiPoint
from sortedcontainers import SortedDict

from typing import Tuple, Dict, Optional


# credit to Geoff Boeing at
# https://geoffboeing.com/2014/08/clustering-to-reduce-spatial-data-set-size/
# for the tutorial
# here coords should be a dictionary with key=ID and value=tuple of lat,long
# epsilon=max distance (in km) between points in a cluster
def dbscan_cluster(coords: Dict[str, Tuple[float, float]], epsilon: float,
                   min_samples: int, out: Optional[str] = None) -> SortedDict:
    """
    Cluster geographic coordinates using DBSCAN algorithm.

    Parameters:
    coords (Dict[str, Tuple[float, float]]): A dictionary of coordinate pairs
                                             with the key as the location
                                             identifier and values as tuples of
                                             (latitude, longitude).
    epsilon (float): Maximum distance between two samples for them to be
                     considered as in the same neighborhood.
    min_samples (int): The number of samples (or total weight) in a
                       neighborhood for a point to be considered as a core
                       point.
    out (Optional[str]): Output file path to save cluster labels.

    Returns:
    SortedDict: A sorted dictionary of cluster labels with keys as cluster name
                and values as lists of location identifiers in that cluster.

    Raises:
    ValueError: If any coordinate value is missing or not a float.

    Examples:
    >>> coords = {'A': (33.4484, -112.0740), 'B': (37.7749,
                  -122.4194), 'C': (40.7128, -74.0060)}
    >>> dbscan_cluster(coords, 300, 2)
    SortedDict({'DB_-1': ['C'], 'DB_0': ['A'], 'DB_1': ['B']})
    """
    # Check if coordinates have correct values
    for k, v in coords.items():
        if not isinstance(v, tuple) or len(v) != 2:
            raise ValueError(
                f"Coordinate value for '{k}' is missing or",
                "not a tuple of size 2."
            )
        if not all(isinstance(coord, float) for coord in v):
            raise ValueError(
                f"Coordinate value for '{k}' has non-float values."
            )

    # Convert coordinates to a numpy array for clustering
    points = np.array(list(coords.values()))
    kms_per_radian = 6371.0088
    eps = epsilon / kms_per_radian

    # Run DBSCAN clustering
    db = DBSCAN(
        eps=eps, min_samples=min_samples, algorithm='ball_tree',
        metric='haversine').fit(np.radians(points))

    # Save cluster labels
    cluster_labels = db.labels_

    # Build SortedDict to return
    popmap = SortedDict()
    for i, k in enumerate(coords.keys()):
        pop = "DB_" + str(cluster_labels[i])
        if pop not in popmap:
            names = [k]
            popmap[pop] = names
        else:
            popmap[pop].append(k)

    # Save cluster labels to file if output path provided
    if out:
        with open(out, 'w') as f:
            for i, k in enumerate(coords.keys()):
                f.write(f"{k}\t{cluster_labels[i]}\n")

    return popmap


def coords_to_matrix(coords):
    """
    Convert a dictionary of geographic coordinates to a NumPy array.

    Args:
    - coords (dict): A dictionary containing the coordinates in the format
                     {id:(long, lat)}.

    Returns:
    - A NumPy array of coordinates in the format [[long, lat], ...].

    """
    if not coords:
        return np.array([], dtype=float).reshape(0, 2)

    return pd.DataFrame(
        [[coords[k][0], coords[k][1]] for k in coords],
        columns=["long", "lat"]).to_numpy()


def coords_to_dataframe(coords):
    """
    Convert a dictionary of geographic coordinates to a Pandas DataFrame.

    Args:
    - coords (dict): A dictionary containing the coordinates in the format
                     {id:(long, lat)}.

    Returns:
    - A Pandas DataFrame of coordinates with columns "long" and "lat".

    """
    df = pd.DataFrame(
        [[coords[k][0], coords[k][1]] for k in coords],
        index=list(coords.keys()),
        columns=["long", "lat"]
    )
    if df.empty:
        df = df.astype({'long': 'float64', 'lat': 'float64'})
    return df


# function to find the centroid of a set of points
# requires a SortedDict of coordinates and a SortedDict giving population IDs
# """Coords:
#     key     value
#     SampleName    Tuple(Lat, Long)
#     popmap:
#     PopulationName    list(SampleName,...)
# """
def get_cluster_centroid(coords, popmap, out=None):
    centroids = SortedDict()
    ofh = None
    if out:
        ofh = out + ".clusterCentroids.txt"
    log = ""
    for i, pop in enumerate(list(popmap.keys())):
        cluster = get_pop_coords_matrix(coords, popmap[pop])
        if len(cluster) < 1:
            raise ValueError(f"No coordinates in cluster: {pop}")

        # add cluster to logfile (if provided)
        log = log + "Population=" + str(i) + "\n"
        log = log + str(cluster) + "\n"

        # get centroid point
        centroid = (
            MultiPoint(cluster).centroid.x, MultiPoint(cluster).centroid.y
        )
        log = log + "Centroid=" + str(centroid) + "\n"
        centroids[pop] = centroid

    # write optional logfile
    if out:
        f = open(ofh, "w")
        f.write(log)
        f.close()
    return centroids


def get_pop_coords_matrix(d, subset):
    if not isinstance(d, dict):
        raise TypeError("Expected 'd' to be a dictionary")
    if not isinstance(subset, list):
        raise TypeError("Expected 'subset' to be a list")

    ret = []
    for k in subset:
        if k in d:
            if isinstance(d[k], tuple) and len(d[k]) == 2:
                ret.append([d[k][0], d[k][1]])
            else:
                raise ValueError(f"Invalid coordinate format for key: {k}")
        else:
            raise ValueError(f"Key not found in dictionary: {k}")

    return np.array(ret, dtype=float).reshape(-1, 2)


# function plots clustered coordinates given a SortedDict of coords and a
# population map
def plot_clustered_points(point_coords, popmap, out, centroids=None):
    # set output file name
    ofh = out+".clusteredPoints.pdf"
    sns.set(style="ticks")

    # get 1D popmap
    pmap = flatten_popmap(popmap)

    df = pd.DataFrame(
        [[ind, pmap[ind], point_coords[ind][0], point_coords[ind][1]]
         for ind in point_coords], columns=["sample", "pop", "long", "lat"])

    ax = sns.scatterplot(x="long", y="lat", hue="pop", palette="Set2", data=df)

    # plot centroid positions if available
    if centroids:
        cdf = pd.DataFrame(
            [[p, centroids[p][0], centroids[p][1]] for p in centroids],
            columns=["pop", "long", "lat"])
        sns.scatterplot(
            x="long", y="lat", hue="pop", palette="Set2", data=cdf,
            legend=False, marker="X", ax=ax
        )

    plt.savefig(ofh)
    plt.clf()


def plot_histogram(dat, out):
    """
    Plots a histogram of snap distances and saves to a PDF file.

    Args:
        dat (list): List of snap distances.
        out (str): Output file path without extension.

    Returns:
        None
    """
    of = f"{out}.snapDistances.pdf"
    sns.set(style="ticks")

    x = pd.Series(dat, name="Snap distance (km)")
    _ = sns.displot(x, kde=True, rug=True)
    plt.savefig(of)
    plt.clf()


def flatten_popmap(popmap):
    """
    Flattens the popmap dictionary from a form of key=pop; value=list(inds) to
    key=ind; value=pop.

    Args:
        popmap (dict): Dictionary in key=pop; value=list(inds) format.

    Returns:
        dict: Dictionary in key=ind; value=pop format.
    """
    new_popmap = {}
    for k, v in popmap.items():
        for i in v:
            new_popmap[i] = k
    return new_popmap
