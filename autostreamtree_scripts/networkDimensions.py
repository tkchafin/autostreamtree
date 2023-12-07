#!/usr/bin/env python
import argparse
import networkx as nx
from pyproj import Proj
from pyproj import Transformer
from scipy.spatial import ConvexHull
import numpy as np


def main():

    # read arguments 
    parser = argparse.ArgumentParser(description='Calculate total lineage length and bounding box area of a network.')
    parser.add_argument('network_file', type=str, help='Path to the pickled Networkx graph file.')
    args = parser.parse_args()

    # Load the pickled Networkx graph
    with open(args.network_file, 'rb') as file:
        network = nx.read_gpickle(file)

    total_length = calculate_total_length(network)
    print(f"Total linear length: {total_length} km")

    bounding_box_area = calculate_bounding_box_area(network)
    print(f"Bounding box area: {bounding_box_area} sq. km")

    convex_hull_area = calculate_convex_hull_area(network)
    print(f"Convex hull area: {convex_hull_area} sq. km")

def calculate_convex_hull_area(network):
    in_proj = Proj(proj='latlong', datum='WGS84')

    longitudes, latitudes = zip(*network.nodes)

    center_lat = (min(latitudes) + max(latitudes)) / 2
    center_lon = (min(longitudes) + max(longitudes)) / 2
    out_proj = Proj(proj='laea', lat_0=center_lat, lon_0=center_lon, datum='WGS84')

    transformer = Transformer.from_proj(in_proj, out_proj)

    coords_projected = [transformer.transform(lon, lat) for lon, lat in zip(longitudes, latitudes)]
    coords_projected_np = np.array(coords_projected)

    hull = ConvexHull(coords_projected_np)
    area = hull.volume / (1000 * 1000)  # Convert from sq. meters to sq. km

    return area

def calculate_total_length(network):
    total_length = 0
    for edge in network.edges(data=True):
        total_length += edge[2]['LENGTH_KM']
    return total_length

def calculate_bounding_box_area(network):
    in_proj = Proj(proj='latlong', datum='WGS84')

    longitudes, latitudes = zip(*network.nodes)
    min_lat, max_lat = min(latitudes), max(latitudes)
    min_lon, max_lon = min(longitudes), max(longitudes)

    center_lat = (min_lat + max_lat) / 2
    center_lon = (min_lon + max_lon) / 2
    out_proj = Proj(proj='laea', lat_0=center_lat, lon_0=center_lon, datum='WGS84')

    transformer = Transformer.from_proj(in_proj, out_proj)

    min_x, min_y = transformer.transform(min_lon, min_lat)
    max_x, max_y = transformer.transform(max_lon, max_lat)

    width = max_x - min_x
    height = max_y - min_y

    return width * height / (1000 * 1000)  # Convert from sq. meters to sq. km

if __name__ == "__main__":
    main()
