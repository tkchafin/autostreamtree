import sys
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt

import autostreamtree.functions as ast
from autostreamtree.params import parseArgs
import autostreamtree.cluster_pops as clust
import autostreamtree.report_refs as ref
import autostreamtree.aggregators as agg

def main():

    params = parseArgs()

    #########################################################
    # Step 1: Reading in datasets
    #########################################################

    print("Running autostream tree", params.run, "workflow\n")
    G = ast.read_network(params.network, params.shapefile)

    # read point coordinate data
    points = pd.read_csv(params.geodb, sep="\t", header=0)
    (point_coords, pop_coords, popmap) = ast.process_samples(params, points, G)

    # read genetic data from VCF
    if params.run != "STREAMDIST":
        seqs = ast.read_vcf(params.vcf, concat=params.concat, popmap=popmap)

    #########################################################
    # Step 2: Calculating genetic distance matrix(-ces)
    # Optional: Calculate separately for each locus
    #########################################################

    gen=None
    pop_gen=None
    # calculate genetic distance matrix
    if params.run != "STREAMDIST" and params.run != "RUNLOCI":
        if not params.genmat:
            print("\nCalculating genetic distances...")
            (gen, pop_gen) = ast.get_gendist_mats(params, point_coords, popmap, seqs)
        else:
            print("\nReading genetic distances from provided matrix:", params.genmat)
            inmat = pd.read_csv(params.genmat, header=0, index_col=0, sep="\t")
            (gen, pop_gen) = ast.parse_input_genmat(params, inmat, point_coords, popmap)

        ast.report_genmats(params, gen, pop_gen, point_coords, pop_coords)

        # exit if this was the only step requested
        if params.run == "GENDIST":
            sys.exit(0)

    #IMPORTANT NOTE: This function will silently skip any loci for which calculations aren't possible (e.g., population having no data)
    if params.run == "RUNLOCI":
        genlist=list()
        popgenlist=list()
        for loc in ast.get_loc_data(seqs):
            try:
                (gen, pop_gen) = ast.get_gendist_mats(params, point_coords, popmap, loc)
                genlist.append(gen)
                popgenlist.append(pop_gen)
            except Exception:
                pass
        print("\nCalculating average of locus-wise genetic distance matrices...")
        if genlist[0] is not None:
            gen = np.mean(genlist, axis=0)
        else:
            gen = genlist[0]
        if popgenlist[0] is not None:
            popgen = np.mean(popgenlist, axis=0)
        else:
            popgen = popgenlist[0]

        ast.report_genmats_list(params, genlist, popgenlist, point_coords, pop_coords)

    #########################################################
    # Step 4: Constructing a minimal subgraph
    #########################################################

    # Extract subgraph if needed
    if params.run != "GENDIST":
        K = ast.parse_subgraph_from_points(params, point_coords, pop_coords, G)


    #########################################################
    # Step 5: Construct stream distance matrix
    # Optional: Mantel test
    #########################################################

    # calculate pairwise observed stream distances and indence matrix
    # calculate incidence matrix X, which takes the form:
    # nrows = rows in column vector form of D
    # ncols = number of collapsed branches in stream network K
    if params.run in ["STREAMDIST", "DISTANCES", "STREAMTREE", "IBD", "ALL", "RUNLOCI"]:
        if params.pop or params.geopop or params.clusterpop:
            points=pop_coords
            gen=pop_gen
        else:
            points=point_coords

        #calculate stream distances and incidence matrix
        (sdist, inc) = ast.get_stream_mats(points, K, params.length_col)
        print("\nStream distances:")
        print(sdist)
        sDF = pd.DataFrame(sdist, columns=list(points.keys()), index=list(points.keys()))
        sDF.to_csv((str(params.out) + ".streamDistMat.txt"), sep="\t", index=True)
        del sDF
        if params.run == "STREAMDIST":
            sys.exit(0)

        #if user requested testing isolation-by-distance:
        if params.run in ['ALL', 'IBD']:
            #HERE: Implement the IBD calculations and plots
            print("\nTesting for isolation-by-distance using Mantel test with",params.permutations,"permutations")
            ast.test_ibd(gen, sdist, params.out, params.permutations)

            #plot geographic x genetic DISTANCES
            ast.plot_gen_by_geo(gen, sdist, params.out)

            #if log request, re-do with logged geographic distances
            if params.and_log:
                print("\nTesting for isolation-by-distance with log geographic distances using Mantel test with",params.permutations,"permutations")
                out=params.out+".log"
                ast.test_ibd(gen, sdist, out, params.permutations, log=True)
                ast.plot_gen_by_geo(gen, sdist, out, log=True)

    #########################################################
    # Step 6: Least-squares fitting
    # Optional: Weighted least-squares
    # Optional: Fit for each locus
    #########################################################

    if params.run in ["STREAMTREE", "ALL", "RUNLOCI"]:
        if params.pop or params.geopop or params.clusterpop:
            gen=pop_gen
        print("\nIncidence matrix:")
        print(inc)
        ofh=params.out+".incidenceMatrix.txt"
        with np.printoptions(precision=0, suppress=True):
            np.savetxt(ofh, inc, delimiter="\t")
        print("Incidence matrix dimensions:")
        print(inc.shape)

        #fit least-squares branch lengths
        if params.run != "RUNLOCI":
            print()
            R = ast.fit_least_squares_distances(gen, inc.astype(int), params.iterative, params.out,params.weight)
            print("\nFitted least-squares distances:")
            print(R)
        else:
            print()
            if params.pop or params.geopop or params.clusterpop:
                genlist=popgenlist
            print("Fitting StreamTree distances on per-locus matrices...")
            Rlist = list()
            ast.block_print()
            for gen in genlist:
                r = ast.fit_least_squares_distances(gen, inc.astype(int), params.iterative, params.out,params.weight)
                Rlist.append(r)
            ast.enable_print()

            #aggregate locus distances (mean)
            print("\nCalculating average fitted distances across loci using:",params.loc_agg)
            averageR = np.array([agg.aggregate_dist(params.loc_agg, col) for col in zip(*Rlist)])
            print(averageR)
            #sys.exit()

            #get standard-deviation locus distances as well
            print("\nCalculating standard-deviation of fitted distances across loci using:",params.loc_agg)
            sdR = np.array([agg.aggregate_dist("SD", col) for col in zip(*Rlist)])
            print(sdR)


        #check observed versus fitted distances:
        if params.run == "RUNLOCI":
            R=averageR
        pred=ast.get_fitted_d(points, gen, inc, R)
        print("\nComparing observed versus predicted genetic distances:")
        print(pred)
        ast.output_fitted_d(pred, params.out)


        #########################################################
        # Step 7: Outputs
        #########################################################

        #Now, annotate originate geoDF with dissolved reach IDs
        #also, need to collect up the stream tree fitted D to each dissolved reach
        #finally, could add residuals of fitting D vs LENGTH_KM?
        #maybe include logDxlength, DxlogLength, logDxlogLength as well?

        #get list of all REACHIDs to extract from geoDF
        #edge_data = nx.get_edge_attributes(K,params.reachid_col)
        reach_to_edge = dict()
        i=0
        edges = list()
        #print("K:",len(K.edges())
        for e in K.edges():
            edge_data = K[e[0]][e[1]][params.reachid_col]
            #print(edge_data)
            for r in edge_data:
                reach_to_edge[r] = str(i)
            edges.append(i)
            i+=1
        #print("edges:",edges)
        #print(len(edges))
        del edge_data

        #save reach_to_edge table to file
        r2eDF = pd.DataFrame(list(reach_to_edge.items()), columns=[params.reachid_col,'EDGE_ID'])
        r2eDF.to_csv((str(params.out)+".reachToEdgeTable.txt"), sep="\t", index=False)

        #read in original shapefile as geoDF and subset it
        print("\nExtracting attributes from original dataframe...")
        geoDF = gpd.read_file(params.shapefile)
        mask = geoDF[params.reachid_col].isin(list(reach_to_edge.keys()))
        temp = geoDF.loc[mask]
        del mask
        del reach_to_edge

        #join EDGE_ID to geoDF
        geoDF = temp.merge(r2eDF, on=params.reachid_col, how='left')
        for col in geoDF.columns:
            if col.endswith('_x'):
                geoDF.drop(col, axis=1, inplace=True)
            elif col.endswith('_y'):
                geoDF.rename(columns={col: col.rstrip('_y')}, inplace=True)
        print(geoDF)
        del temp
        del r2eDF

        #plot by edge ID
        base = geoDF.plot(column="EDGE_ID", cmap = "prism")
        coords = clust.coords_to_dataframe(points)
        geo_coords = gpd.GeoDataFrame(coords, geometry=gpd.points_from_xy(coords.long, coords.lat))
        geo_coords.plot(ax=base, marker='o', color='black', markersize=10, zorder=10)

        plt.title("Stream network colored by EDGE_ID")
        plt.savefig((str(params.out)+".streamsByEdgeID.pdf"))

        #add in fitted distances & plot
        print("\nPlotting StreamTree fitted distances and writing new shapefile...")

        if params.run == "RUNLOCI":
            fittedD = pd.DataFrame({'EDGE_ID':list(edges), 'fittedD':R})
            sdD = pd.DataFrame({'EDGE_ID':list(edges), 'stdevD':sdR})
            i=1
            for locfit in Rlist:
                name="locD_" + str(i)
                fittedD[name] = locfit
                i+=1
        else:
            fittedD = pd.DataFrame({'EDGE_ID':list(edges), 'fittedD':R})
        geoDF['EDGE_ID'] = geoDF['EDGE_ID'].astype(int)
        #geoDF.drop(columns=["EDGE_ID_x", "EDGE_ID_y"], inplace=True)
        geoDF = geoDF.merge(fittedD, on='EDGE_ID', how="left")
        for col in geoDF.columns:
            if col.endswith('_x'):
                geoDF.drop(col, axis=1, inplace=True)
            elif col.endswith('_y'):
                geoDF.rename(columns={col: col.rstrip('_y')}, inplace=True)
        #print(geoDF)
        if params.run == "RUNLOCI":
            geoDF = geoDF.merge(sdD, on='EDGE_ID')
            geoDF.plot(column="stdevD", cmap = "RdYlGn_r", legend=True)
            plt.title("Stream network colored by standard deviation of StreamTree fitted distances")
            plt.savefig((str(params.out)+".streamsBystdevD.pdf"))
        geoDF.plot(column="fittedD", cmap = "RdYlGn_r", legend=True)
        plt.title("Stream network colored by StreamTree fitted distances")
        plt.savefig((str(params.out)+".streamsByFittedD.pdf"))

        #output a final annotated stream network layer
        geoDF.to_csv((str(params.out)+".streamTree.txt"), sep="\t", index=False)

        if geoDF.shape[1] > 2045:
            print("Too many columns to write shapefile (hard limit of 2046 columns). Only writing first 2046 columns to shapefile attribute table; full output can be left-joined from $out.streamtree.txt.")
            geoDF.iloc[:, : 2045].to_file((str(params.out)+".streamTree.shp"))
        else:
            geoDF.to_file((str(params.out)+".streamTree.shp"))

    refs = ref.fetch_references(params)
    print(refs)

    print("\nDone!\n")

#Call main function
if __name__ == '__main__':
    main()
