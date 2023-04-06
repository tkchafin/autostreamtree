
import sys
import os
import geopy
import itertools
import traceback
import math
import getopt
import scipy
import momepy
import pandas as pd
import geopandas as gpd
import numpy as np
import networkx as nx
import seaborn as sns
from pysam import VariantFile
from scipy import stats
from os import listdir
from os.path import isfile, join
from sortedcontainers import SortedDict
from sklearn.linear_model import LinearRegression
from shapely.geometry import LineString, point, Point
from networkx import NodeNotFound
from networkx import NetworkXNoPath
import matplotlib.pyplot as plt
from geopy.distance import geodesic
import pickle
from math import radians, degrees, sin, cos, asin, acos, sqrt

import autostreamtree.cluster_pops as clust
import autostreamtree.sequence as seq
import autostreamtree.genetic_distances as gendist
import autostreamtree.report_refs as ref
import autostreamtree.aggregators as agg
import autostreamtree.Mantel as Mantel

from typing import List, Tuple, Dict, Any, Union, Optional

def read_vcf(vcf, concat="none", popmap=None):
    """
    Reads a VCF file and returns a dictionary of sample genotypes.

    ARgs:
        vcf: Path to the input VCF file.
        concat: Specifies the concatenation mode for genotypes. Options are "all", "loc", and "none".
                   "all": Concatenate genotypes of all loci for each sample.
                   "loc": Concatenate genotypes within the same chromosome for each sample.
                   "none": Do not concatenate genotypes.
        popmap: Optional dictionary that maps populations to a list of samples. If provided,
                   only samples in the popmap will be retained in the output dictionary.
    Returns:
        A dictionary with sample names as keys and lists of genotypes as values.
    """

    bcf_in = VariantFile(vcf)

    # get all samples in the VCF
    vcf_samples = list(bcf_in.header.samples)

    # set up data dict
    dat=dict()
    samples = list((bcf_in.header.samples))
    for s in samples:
        if concat == "all":
            dat[s] = list()
            dat[s].append(["",""])
        else:
            dat[s] = list()

    # if popmap, make list of samples to drop that aren't in a pop
    if popmap:
        keep = list()
        for pop in popmap:
            keep.extend(popmap[pop])
        keep = [s for s in keep if s in vcf_samples]
        bcf_in.subset_samples(keep)

    chrom="FIRST"
    for record in bcf_in.fetch():
        for i, sample in enumerate(record.samples):
            if concat=="all":
                loc = seq.decode(record.samples[i]['GT'], record.ref, record.alts, as_list=True)
                dat[sample][-1][0]=dat[sample][-1][0]+loc[0]
                dat[sample][-1][1]=dat[sample][-1][1]+loc[1]
            elif concat=="loc":
                if record.chrom != chrom:
                    dat[sample].append(["",""])
                loc = seq.decode(record.samples[i]['GT'], record.ref, record.alts, as_list=True)
                dat[sample][-1][0]=dat[sample][-1][0]+loc[0]
                dat[sample][-1][1]=dat[sample][-1][1]+loc[1]
            else:
                loc = seq.decode(record.samples[i]['GT'], record.ref, record.alts)
                dat[sample].append(loc)
        chrom=record.chrom
    if concat != "none":
        for sample in dat:
            dat[sample] = ["/".join(x) for x in dat[sample]]
    for sample in list(dat.keys()):
        if len(dat[sample]) < 1:
            del dat[sample]
        elif len(dat[sample]) == 1 and dat[sample][0][0] == "":
            del dat[sample]
    return(dat)

def read_network(network, shapefile):
    """
    Reads a network from a saved file or builds a network from a shapefile.

    Args:
        network: Path to the saved network file (pickle format). If provided, the function will read
                    the network from this file.
        shapefile: Path to the shapefile to build the network from. This is used if the `network`
                      parameter is not provided.
    Returns:
        A NetworkX Graph object representing the network.
    """

    # Check if a saved network file is provided
    if network:
        print("Reading network from saved file: ", network)
        # Read the network from the saved file and convert it to an undirected graph
        G = nx.Graph(nx.read_gpickle(network).to_undirected())
    else:
        # If no saved network file is provided, build the network from the shapefile
        print("Building network from shapefile:", shapefile)
        print("WARNING: This can take a while with very large files!")
        # Read the shapefile
        rivers = gpd.read_file(shapefile)
        # Convert the GeoDataFrame to a NetworkX Graph object
        G = momepy.gdf_to_nx(rivers, approach="primal", directed=False, multigraph=False)

    return G

def parse_subgraph_from_points(params, point_coords, pop_coords, G):
    """
    Extracts a subgraph from a given graph based on input points.

        params: A custom object containing various input parameters.
        point_coords: A list of point coordinates to use for extracting the subgraph.
        pop_coords: A list of population coordinates to use for extracting the subgraph.
        G: A NetworkX Graph object representing the input graph.

    Returns:
        A NetworkX Graph object representing the extracted subgraph.
    """

    # Choose the appropriate set of points based on input parameters
    if params.pop or params.geopop or params.clusterpop:
        points = pop_coords
    else:
        points = point_coords

    # Output points to a table
    p = get_point_table(points)
    p.to_csv((str(params.out) + ".pointCoords.txt"), sep="\t", index=False)
    del p

    # First pass extracts a subgraph from the master shapefile graph
    print("\nExtracting full subgraph...")
    ktemp = path_subgraph(G, points, extract_full_subgraph, params.reachid_col, params.length_col)
    del G

    # Second pass to simplify subgraph and collapse redundant nodes
    print("\nMerging redundant paths...\n")
    K = path_subgraph(ktemp, points, extract_minimal_subgraph, params.reachid_col, params.length_col)

    # Grab real coordinates as node positions for plotting
    pos = {n: n for n in K.nodes}

    # Make a color map to color sample points and junctions differently
    color_map = ["blue" if node in points.values() else "black" for node in K]

    # Draw networkx
    nx.draw_networkx(K, pos, with_labels=False, node_color=color_map, node_size=50)

    # Get LENGTH_KM attributes for labelling edges
    edge_labels = nx.get_edge_attributes(K, params.length_col)
    for e in edge_labels:
        edge_labels[e] = "{:.2f}".format(edge_labels[e])

    nx.draw_networkx_edge_labels(K, pos, edge_labels=edge_labels, font_size=6)

    # Save minimized network to file (unless we already read from one)
    if not params.network or params.overwrite:
        net_out = str(params.out) + ".network"
        nx.write_gpickle(K, net_out, pickle.HIGHEST_PROTOCOL)
        net_full_out = str(params.out) + ".full.network"
        nx.write_gpickle(ktemp, net_full_out, pickle.HIGHEST_PROTOCOL)
    else:
        print("NOTE: Not over-writing existing network. To change this, use --overwrite")

    network_plot = str(params.out) + ".subGraph.pdf"
    plt.savefig(network_plot)

    del ktemp
    return K


def report_genmats(params, gen, pop_gen, point_coords, pop_coords):
    """
    Prints genetic distance matrices and writes them to files.

    ARgs:
        params: A custom object containing various input parameters.
        gen: A NumPy array representing the individual genetic distance matrix.
        pop_gen: A NumPy array representing the population genetic distance matrix.
        point_coords: A dictionary containing individual point coordinates.
        pop_coords: A dictionary containing population point coordinates.
    """

    # If the individual genetic distance matrix is not None, print and write it to file
    if gen is not None:
        print("Genetic distances:")
        np.set_printoptions(precision=3)
        print(gen, "\n")

        # Write individual genetic distances to file
        ind_genDF = pd.DataFrame(gen, columns=list(point_coords.keys()), index=list(point_coords.keys()))
        ind_genDF.to_csv((str(params.out) + ".indGenDistMat.txt"), sep="\t", index=True)

    # If the population genetic distance matrix is not None, print and write it to file
    if pop_gen is not None:
        print("Population genetic distances:")
        np.set_printoptions(precision=3)
        print(pop_gen, "\n")

        # Write population genetic distances to file
        pop_genDF = pd.DataFrame(pop_gen, columns=list(pop_coords.keys()), index=list(pop_coords.keys()))
        pop_genDF.to_csv((str(params.out) + ".popGenDistMat.txt"), sep="\t", index=True)
        del pop_genDF


def get_loc_data(seqs):
    """
    Generator function that yields a dictionary of individual loci data for each locus in the input sequences.

    Args:
        seqs: A dictionary containing sequences as values and individual identifiers as keys.
    
    Returns:
        A generator that yields a dictionary with individual identifiers as keys and a list containing the
             corresponding locus as the value.
    """

    # Iterate through the loci in the sequences
    for loc in range(0, len(seqs[list(seqs.keys())[0]])):
        d = dict()

        # Iterate through the individuals in the sequences
        for ind in seqs.keys():
            # Add the locus data for the individual to the dictionary
            d[ind] = [seqs[ind][loc]]

        # Yield the dictionary containing the locus data for the current locus
        yield(d)


def report_genmats_list(params, genlist, popgenlist, point_coords, pop_coords):
    """
    Writes individual and population genetic distance matrices to files for each locus in genlist and popgenlist.

    Args:
        params: A namespace object containing parameters, including the output directory.
        genlist: A list of individual genetic distance matrices for each locus.
        popgenlist: A list of population genetic distance matrices for each locus.
        point_coords: A dictionary containing individual point coordinates.
        pop_coords: A dictionary containing population point coordinates.
    """

    # Create an output directory for the genetic distance matrices
    dir = str(params.out) + "_locmats"
    os.makedirs(dir, exist_ok=True)

    i = 0
    # Iterate through the individual genetic distance matrices
    for gen in genlist:
        if gen is not None:
            # Write individual genetic distances to a file
            ind_genDF = pd.DataFrame(gen, columns=list(point_coords.keys()), index=list(point_coords.keys()))
            ind_genDF.to_csv((str(dir) + "/loc_" + str(i) + ".indGenDistMat.txt"), sep="\t", index=True)
            del ind_genDF
            i += 1

    j = 0
    # Iterate through the population genetic distance matrices
    for pop_gen in popgenlist:
        if pop_gen is not None:
            # Write population genetic distances to a file
            pop_genDF = pd.DataFrame(pop_gen, columns=list(pop_coords.keys()), index=list(pop_coords.keys()))
            pop_genDF.to_csv((str(dir) + "/loc_" + str(j) + ".popGenDistMat.txt"), sep="\t", index=True)
            del pop_genDF
            j += 1

def block_print():
    """
    Disables standard output by redirecting it to a null device, effectively blocking any print statements.
    """
    sys.stdout = open(os.devnull, 'w')


def enable_print():
    """
    Restores standard output to its original state, allowing print statements to be displayed again.
    """
    sys.stdout = sys.__stdout__


def parse_input_genmat(params, inmat, point_coords, popmap):
    """
    Parses an input genetic distance matrix and verifies if it matches the user input
    parameters. Aggregates individual distances if required by the user input.
    
    Args:
        params: Input parameters provided by the user.
        inmat (pd.DataFrame): The input genetic distance matrix.
        point_coords (dict): Dictionary containing point coordinates.
        popmap (dict): Dictionary containing the population map.

    Returns:
        tuple: A tuple containing the genetic distance matrix (gen) and
        population genetic distance matrix (pop_gen).
    """
    gen = None
    pop_gen = None
    if params.coercemat:
        inmat[inmat < 0.0] = 0.0
    if set(list(inmat.columns.values)) != set(list(inmat.index.values)):
        print(inmat.columns.values)
        print(inmat.index.values)
        print("Oh no! Input matrix columns and/ or rows don't appear to be labelled. Please provide an input matrix with column and row names!")
        sys.exit(1)
    else:
        agg=False
        #first check if it fits whatever the user input was (i.e. --pop)
        if params.pop:
            if len(inmat.columns) != len(popmap.keys()):
                print("Found",str(len(inmat.columns)), "columns in provided matrix. This doesn't match number of populations from popmap.")
                if (len(inmat.columns)) != len(point_coords):
                    print("Doesn't match number of individuals either! Please check your matrix format.")
                    sys.exit(1)
                else:
                    print("Assuming input matrix has individual distances... Aggregating using the following method (--pop_agg):", str(params.pop_agg))
                    agg=True
            else:
                #re-order using pop orders
                inmat = inmat.reindex(list(popmap.keys()))
                inmat = inmat[list(popmap.keys())]
                pop_gen = inmat.to_numpy()
                del(inmat)
        elif params.geopop or params.clusterpop:
            if (len(inmat.columns)) != len(point_coords):
                print("Found",str(len(inmat.columns)), "columns in provided matrix. This doesn't match number of individuals.")
                print("When using --geopop or --clusterpop, the provided matrix must represent individual-level distances.")
                sys.exit(1)
            else:
                #re-order using pop orders
                inmat = inmat.reindex(list(point_coords.keys()))
                inmat = inmat[list(point_coords.keys())]
                gen = inmat.to_numpy()
                agg = True
                del(inmat)
        else:
            if (len(inmat.columns)) != len(point_coords):
                print("Found",str(len(inmat.columns)), "columns in provided matrix. This doesn't match number of individuals.")
                sys.exit(1)
            else:
                #re-order using pop orders
                inmat = inmat.reindex(list(point_coords.keys()))
                inmat = inmat[list(point_coords.keys())]
                gen = inmat.to_numpy()
                del(inmat)
        #if --geopop or --clusterpop, it should be an ind matrix
        #if so, need to aggregate according to --pop_agg
        #print(pop_gen)
        if agg:
            print("Aggregating user-provided individual-level distance matrix using:",params.pop_agg)
            pop_gen = gendist.get_pop_genmat("PDIST", gen, popmap, point_coords, seqs, pop_agg=params.pop_agg, loc_agg=params.loc_agg, ploidy=params.ploidy, global_het=params.global_het)

    return (gen, pop_gen)

def read_popmap(popmap):
    """
    Reads a population map file and returns a dictionary with individuals as keys
    and populations as values.

    Args:
        popmap (str): Path to the population map file.

    Returns:
        dict: A dictionary with individuals as keys and populations as values.
    """
    popdict = dict()
    with open(popmap, "r") as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            cols = line.split()
            ind = cols[0]
            pop = cols[1]
            popdict[ind] = pop
    return(popdict)

def process_samples(params, points, G):
    """
    Processes input sample data by snapping points to a graph, calculating coordinates,
    and processing populations if required. This function is separated to declutter __main__.
    
    Args:
        params: Input parameters provided by the user.
        points (pd.DataFrame): DataFrame containing sample points.
        G (networkx.Graph): Graph object representing the road network.

    Returns:
        tuple: A tuple containing point coordinates, population coordinates, and the population map.
    """
    popmap = SortedDict()
    point_coords = SortedDict()
    pop_coords = SortedDict()
    snapDists = dict()
    verb = True
    first = True

    if params.pop:
        popmap_temp = read_popmap(params.pop)
        mask = points[points.columns[0]].isin(popmap_temp)
        points = points[mask]
        
    for idx, row in points.iterrows():
        name = None
        data = None
        row[1] = float(row[1])
        row[2] = float(row[2])
        if params.run == "GENDIST":
            name = row[0]
            data = tuple([row[2], row[1]])
        else:
            if not params.pop and not params.clusterpop:
                #print(tuple([row[3], row[2]]))
                #--geopop and individual-level snap coordinates to nodes here
                node = snap_to_node(G, tuple([row[2], row[1]]))
                snapDists[row[0]] = great_circle(node[0], node[1], row[2], row[1])
            else:
                #if pop or clusterpop, extract centroid later
                node = tuple([row[2], row[1]])
            #print(node)
            data = node
            name = row[0]
            #point_labels[node]=str(row[0])
        point_coords[name] = data
        #seq_data = parseLoci(params, list(row[4:]), verbose=verb)
        #print(seq_data)

        # Process population-level analyses
        if params.geopop:
            if point_coords[name] not in popmap:
                l = [name]
                popmap[point_coords[name]] = l
            else:
                popmap[point_coords[name]].append(row[0])
        elif params.pop:
            if popmap_temp[row[0]] not in popmap:
                l = [name]
                popmap[popmap_temp[row[0]]] = l
            else:
                popmap[popmap_temp[row[0]]].append(name)
    #print("Found",numLoci,"loci.\n")
    #points["node"]=point_coords

    print("Read",str(len(point_coords.keys())),"individuals.")
    #print(list(point_coords.keys()))
    print()

    if params.pop or params.geopop:
        print("Read",str(len(popmap.keys())),"populations.")
        #print(list(popmap.keys()))
        print()

    """
    For population-level analyses, generate population maps and centroids here
    according to user-input options: --pop, --geopop, --clusterpop
    """
    #get population centroid
    if params.pop or params.geopop or params.clusterpop:
        if params.clusterpop:
            #create population clusters using DBSCAN
            print("Running DBSCAN clustering with min_samples=",params.min_samples,"and epsilon=",params.epsilon)
            popmap=clust.dbscan_cluster(point_coords, params.epsilon, params.min_samples)
            num_clusters=len(popmap.keys())
            print("Found",str(num_clusters),"clusters!")
            print(popmap)
            print("\n")

            #calculate centroids for clusters
            pop_temp=clust.get_cluster_centroid(point_coords, popmap, params.out)

        elif params.pop or params.geopop:
            #popmap generated earlier when parsing input file!
            #still need to calculate centroids:
            print("Calculating population centroids...")
            pop_temp=clust.get_cluster_centroid(point_coords, popmap, params.out)
            #note in the case of --geopop the centroid is the joint snapped-to location

        #now, snap pop_coords to nodes
        pop_coords = SortedDict()
        for p in pop_temp:
            node = snap_to_node(G, pop_temp[p])
            #print(node[0], node[1], pop_temp[p][0], pop_temp[p][1])
            snapDists[p] = great_circle(node[0], node[1], pop_temp[p][0], pop_temp[p][1])
            pop_coords[p]=node
        #write popmap to file
        flat = clust.flatten_popmap(popmap)
        temp = pd.DataFrame(list(flat.items()), columns=['IND_ID', 'POP_ID'])
        temp.to_csv((str(params.out) + ".popmap.txt"), sep="\t", index=False)
        del flat
        del temp

        #plot grouped samples
        #TODO: for --geopop maybe plot original coordinates with "snap" as centroid here??
        clust.plot_clustered_points(point_coords, popmap, params.out, pop_coords)

    # Plot histogram of snap distances and write snap distances to a file
    clust.plot_histogram(list(snapDists.values()), params.out)
    dtemp = pd.DataFrame(list(snapDists.items()), columns=['name', 'km'])
    dtout = str(params.out) + ".snapDistances.txt"
    dtemp.to_csv(dtout, sep="\t", index=False)
    del dtemp
    del dtout
    del snapDists

    #return everything
    return(point_coords, pop_coords, popmap)


def get_gendist_mats(params, point_coords, popmap, seqs):
    """Returns population genetic distance matrices.

    Args:
        params (object): An object that contains parameters for the analysis.
        point_coords (list): A list of coordinates for the sampled points.
        popmap (dict): A dictionary that maps each sample to its corresponding population.
        seqs (list): A list of DNA sequences for the sampled points.

    Returns:
        tuple: A tuple containing two matrices. The first is a pairwise distance matrix for all
        samples. The second is a matrix of pairwise genetic distances between populations.

    Raises:
        SystemExit: If distance metric is not possible without population data.

    """
    gen = None  # Initialize variable to hold pairwise distance matrix
    pop_gen = None  # Initialize variable to hold population genetic distance matrix

    if params.dist in ["PDIST", "TN84", "TN93", "K2P", "JC69"]:
        # Calculate pairwise distance matrix using selected method
        gen = gendist.get_genmat(params.dist, point_coords, seqs, ploidy=params.ploidy, het=params.het, loc_agg=params.loc_agg)

        if params.pop or params.geopop or params.clusterpop:
            print("Aggregating pairwise population genetic distances from individual distances using:",params.pop_agg)

    else:
        # If distance metric requires population data, but none is provided
        if not params.pop and not params.geopop:
            print("ERROR: Distance metric",params.dist,"not possible without population data.")
            sys.exit(1)

    # Calculate population genetic distance matrix
    if params.pop or params.geopop or params.clusterpop:
        pop_gen = gendist.get_pop_genmat(params.dist, gen, popmap, point_coords, seqs, pop_agg=params.pop_agg, loc_agg=params.loc_agg, ploidy=params.ploidy, global_het=params.global_het)

    # Return pairwise distance matrix and population genetic distance matrix
    return(gen, pop_gen)


def get_point_table(points):
    """Returns a pandas DataFrame from a dictionary of points.

    Args:
        points (dict): A dictionary of points with their latitude and longitude values.

    Returns:
        pandas.DataFrame: A dataframe with columns for 'sample', 'lat', and 'long'.

    """
    temp = []  # Initialize list to hold temporary values
    for p in points:
        # Append values for each point to temporary list
        temp.append([p, points[p][1], points[p][0]])
    p = pd.DataFrame(temp, columns=['sample', 'lat', 'long'])  # Create dataframe from temporary list
    return p


def r2(x, y):
    """Returns the Pearson correlation coefficient squared between two arrays.

    Args:
        x (array): An array of values.
        y (array): An array of values.

    Returns:
        float: The squared Pearson correlation coefficient between the two arrays.

    """
    return (stats.pearsonr(x, y)[0] ** 2)  # Calculate Pearson correlation coefficient squared and return it


def get_fitted_d(points, genmat, inc, r):
    """Calculates predicted genetic distances based on fitted streamtree distances.

    Args:
        points (dict): A dictionary of points with their latitude and longitude values.
        genmat (ndarray): A pairwise genetic distance matrix.
        inc (ndarray): An incidence matrix representing the presence or absence of streams for each point.
        r (ndarray): A fitted streamtree distance matrix.

    Returns:
        pandas.DataFrame: A dataframe with columns for 'from', 'to', 'observed_D', 'predicted_D', and 'abs_diff'.

    """
    rows = []  # Initialize list to hold temporary values
    names = list(points.keys())  # Get names of points

    # Iterate over all pairwise combinations of points
    inc_row=0 #tracks what ROW we are in the incidence matrix (streams are columns!)
    for ia, ib in itertools.combinations(range(0,len(points)),2):
        obs=genmat[ia,ib]
        inc_streams=inc[inc_row,]
        pred_dist=np.sum(r[inc_streams==1])
        inc_row+=1
        rows.append([names[ia], names[ib], obs, pred_dist, np.abs(obs-pred_dist)])
     # Create dataframe from temporary list and return it
    D=pd.DataFrame(rows, columns=['from','to','observed_D', 'predicted_D', 'abs_diff'])
    return(D)

def plot_gen_by_geo(gen, sdist, out, log=False):
    """
    Plots genetic distance against geographic distance to visualize isolation by distance.

    Args:
        gen (numpy.ndarray): Genetic distance matrix.
        sdist (numpy.ndarray): Spatial distance matrix.
        out (str): Output file prefix for the generated plot.
        log (bool, optional): If True, the geographic distance axis will be log-transformed. Defaults to False.
    """

    genetic_distance = get_lower_tri(gen)
    geographic_distance = get_lower_tri(sdist)

    # Create a DataFrame from the geographic_distance and genetic_distance arrays
    data = pd.DataFrame({'Geographic Distance': geographic_distance, 'Genetic Distance': genetic_distance})

    if not log:
        sns.jointplot(data=data, x='Geographic Distance', y='Genetic Distance', kind="reg")
        plt.savefig(str(out) + ".isolationByDistance.pdf")
    else:
        geographic_distance = replace_zeroes(geographic_distance)
        log_geo = np.log(geographic_distance)
        data['Log Geographic Distance'] = log_geo

        sns.jointplot(data=data, x='Log Geographic Distance', y='Genetic Distance', kind="reg")
        plt.savefig(str(out) + ".isolationByDistance.pdf")

    del geographic_distance
    del genetic_distance


def great_circle(lon1, lat1, lon2, lat2, thresh=0.0000001):
    """
    Calculates the great circle distance between two points on a sphere, using their longitudes and latitudes.

    Args:
        lon1 (float): Longitude of the first point.
        lat1 (float): Latitude of the first point.
        lon2 (float): Longitude of the second point.
        lat2 (float): Latitude of the second point.
        thresh (float, optional): Threshold for determining if the points are the same. Defaults to 0.0000001.

    Returns:
        float: Great circle distance in kilometers.
    """
    if (abs(lon1 - lon2)) < thresh and (abs(lat1 - lat2)) < thresh:
        return 0.0
    else:
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
        return 6371 * (
            acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
        )

def get_lower_tri(mat):
    """
    Extracts the lower triangular elements from a square matrix.

    Args:
        mat (numpy.ndarray): Input square matrix.

    Returns:
        numpy.ndarray: 1D array containing the lower triangular elements of the input matrix.
    """
    n = mat.shape[0]
    i = np.tril_indices(n, -1)
    return mat[i]

def replace_zeroes(data):
    """
    Replaces zeroes in the input array with the smallest non-zero value.

    Args:
        data (numpy.ndarray): Input array.

    Returns:
        numpy.ndarray: Array with zeroes replaced by the smallest non-zero value.
    """
    min_nonzero = np.min(data[np.nonzero(data)])
    data[data == 0] = min_nonzero
    return data


#function computes Mantel test using various transformations
def test_ibd(gen, geo, out, perms, log=False):
    #get flattened lower triangle of each matrix
    gen = get_lower_tri(gen)
    geo = get_lower_tri(geo)

    if log==True:
        geo=replace_zeroes(geo)
        geo=np.log(geo)

    #non-log pearson
    res=list(Mantel.test(geo, gen, perms=int(perms), method='pearson'))
    rows=list()
    rows.append(['genXgeo','pearson', str(perms), res[0], res[1], res[2]])

    #non-log spearman
    res=list(Mantel.test(geo, gen, perms=int(perms), method='spearman'))
    rows.append(['genXgeo','spearman', str(perms), res[0], res[1], res[2]])

    #print(rows)
    ibd=pd.DataFrame(rows,  columns=['test', 'method', 'perms', 'r', 'p' ,'z'])
    print("Mantel test results:")
    print(ibd)
    ibd.to_csv((str(out) + ".isolationByDistance.txt"), sep="\t", index=False)
    print()

def output_fitted_d(pred, out):
    pred.to_csv((str(out)+".obsVersusFittedD.txt"), sep="\t", index=False)
    sns.jointplot(x="observed_D", y="predicted_D", data=pred, kind="reg")
    plt.savefig((str(out)+".obsVersusFittedD.pdf"))
    del pred
    #plt.show()

def fit_least_squares_distances(D, X, iterative, out, weight=None):
    """
    Computes least-squares branch lengths from a vector of genetic distances D and incidence matrix X.
    When iterative=True, negative distances are constrained to 0 and then recomputed.

    Args:
        D (numpy.ndarray): Vector of genetic distances.
        X (numpy.ndarray): Incidence matrix.
        iterative (bool): Whether to use an iterative approach to constrain negative distances.
        out (str): Output file prefix.
        weight (numpy.ndarray, optional): Weight matrix for weighted least squares optimization. Defaults to None.

    Returns:
        numpy.ndarray: Least-squared optimized distances.
    """
    num_segments = (np.size(X,1))
    #print(num_segments)
    ls = np.zeros(num_segments)
    d = vectorize_mat(D)

    #calculate weights matrix and write to file
    W=generate_weights_matrix(d, weight)
    print("Weights matrix:")
    print(W)
    #ofh=out+".weightsMatrix.txt"
    #np.savetxt(ofh, W, delimiter="\t")

    #weighted least-squares optimization
    ls = np.matmul(np.linalg.inv(np.matmul(np.matmul(X.transpose(),W),X)), np.matmul(np.matmul(X.transpose(), W),d))

    if iterative:
        ls_old=ls
        if(np.count_nonzero(ls<0.0) > 0):
            print("\nLS-optimized distances contain negative values: Using iterative approach to re-calculate...")
        constrains = list() #save indices of all constrained values

        #if negative distances, use iterative procedure to re-calculate
        while (np.count_nonzero(ls<0.0) > 0):
            bad_ind = np.argmin(ls)
            constrains.append(bad_ind)
            #constrain to 0 by removing from incidence matrix
            X = np.delete(X, bad_ind, 1)
            #re-compute values
            ls = np.matmul(np.linalg.inv(np.matmul(np.matmul(X.transpose(),W),X)), np.matmul(np.matmul(X.transpose(), W),d))
        for i in reversed(constrains):
            ls=np.insert(ls, i, 0.0)
        #print(ls)

        #write original and constrained results to log file
        ofh=out+".leastSquaresConstrained.txt"
        df=pd.DataFrame({'LS.original':ls_old, 'LS.constrained':ls})
        df.to_csv(ofh, sep="\t", index=False)

        return(ls)
    else:
        return(ls)

def generate_weights_matrix(d, weight):
    """
    Generates a weights matrix for the least-squares method, where weights are on the diagonals.

    Args:
        d (numpy.ndarray): Vector of genetic distances.
        weight (str): Weighting method to use, options: 'CSE67', 'BEYER74', 'FM67'.

    Returns:
        numpy.ndarray: Weights matrix.
    """
    W = np.zeros((len(d), len(d)), dtype=float)
    row, col = np.diag_indices(W.shape[0])

    if weight.upper() == "CSE67":
        W[row, col] = np.ones(len(d))
    elif weight.upper() == "BEYER74":
        if np.count_nonzero(d == 0) > 0:
            print("WARNING: Divide-by-zero in weighted least-squares (weight=1/D).")
        W[row, col] = np.divide(1.0, d, out=np.zeros_like(d), where=d != 0)
    elif weight.upper() == "FM67":
        if np.count_nonzero(d == 0) > 0:
            print("WARNING: Divide-by-zero in weighted least-squares (weight=1/D^2).")
        W[row, col] = np.divide(1.0, np.square(d), out=np.zeros_like(d), where=d != 0)
    else:
        print(f"ERROR: Weight option {weight} not recognized. Using ordinary least-squares instead.")
        W[row, col] = np.ones(len(d))

    return W


def vectorize_mat(mat):
    """
    Converts a pairwise matrix to a 1D vector.

    Args:
        mat (numpy.ndarray): Pairwise matrix.

    Returns:
        numpy.ndarray: 1D vector of matrix elements.
    """
    size = nCr(np.size(mat, 0), 2)
    vec = np.zeros(size)
    index = 0

    for ia, ib in itertools.combinations(range(0, np.size(mat, 0)), 2):
        vec[index] = mat[ia, ib]
        index += 1

    return vec

def get_stream_mats(points, graph, len_col):
    """
    Computes pairwise stream distances and 0/1 incidence matrix for StreamTree calculations.

    Args:
        points (dict): Dictionary of point indices and their corresponding node IDs in the graph.
        graph (networkx.Graph): NetworkX graph object representing the stream network.
        len_col (str): Attribute name for the length of the edges in the graph.

    Returns:
        tuple: Pair of numpy.ndarray representing the pairwise stream distance matrix and incidence matrix.
    """
    dist = np.zeros((len(points), len(points)))
    inc = np.zeros((nCr(len(points), 2), len(graph.edges())), dtype=int)

    # Establish as NaN
    dist[:] = np.nan

    def dijkstra_weight(left, right, attributes):
        """
        Calculates weights for Dijkstra's shortest path algorithm by inverting the edge length with a small constant to avoid division by zero.
        """
        epsilon = 1e-9
        return 1 / (attributes[len_col] + epsilon)

    index = 0
    for ia, ib in itertools.combinations(range(0, len(points)), 2):
        path = nx.bidirectional_dijkstra(graph, points.values()[ia], points.values()[ib], weight=dijkstra_weight)
        if path:
            dist[ia, ib] = float(sum(path_edge_attributes(graph, path[1], len_col)))
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
    Check if two elements are consecutive in a list, irrespective of their order.

    Args:
        lst (list): The list to search for the pair.
        x (Any): The first element of the pair.
        y (Any): The second element of the pair.

    Returns:
        bool: True if the elements are consecutive in the list, False otherwise.
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

    def dijkstra_weight(left, right, attributes):
        """
        Calculates weights for Dijkstra's shortest path algorithm by inverting the edge length with a small constant to avoid division by zero.
        """
        epsilon = 1e-9
        return 1 / (attributes[len_col] + epsilon)

    p1 = list(nodes.values())[0]
    for p2 in list(nodes.values())[1:]:
        try:
            # Find the shortest path between the two points
            path = nx.bidirectional_dijkstra(graph, p1, p2, weight=dijkstra_weight)

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
    Extracts a simplified subgraph from paths, keeping only terminal and junction nodes.

    Args:
        subgraph (NetworkX Graph): The subgraph to be modified.
        graph (NetworkX Graph): The input graph.
        nodelist (list): The list of nodes.
        id_col (str): The column name for node ID.
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
        curr_edge[id_col].extend([dat[id_col]] if not isinstance(dat[id_col], list) else dat[id_col])
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


