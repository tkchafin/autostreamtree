# autostreamtree

Implementation of Kalinowsky et al. 2008 'StreamTree' software for SNP datasets

### Table of Contents:
1. [Installation](#installation)
2. [Description](#ast)
    1. [Method](#ast_desc)
        1. [StreamTree Background](#ast_background)
    2. [Usage](#usage)
        1. [Options and Help Menu](#ast_help)
        2. [Input File](#ast_table)
        3. [Shapefile format](#ast_input)
        4. [Output files](#ast_output)
        5. [Genetic Distance Methods](#gen)
        6. [Defining Populations](#pops)
    3. [Example Workflows](#ast_workflow)
        1. [Large SNP datasets](#ast_example2)
        2. [Microhaplotypes](#ast_example3)
    4. [Runtimes and Benchmarking](#ast_benchmark)
    5. [References](#ast_refs)
3. [Scripts and Other Useful Features](#tools)
4. [Incorporating autoStreamTree into automated workflows](#snakemake)
5. [Contributing Guidelines](#contributing)


## 1. Installation <a name="installation"></a>

autoStreamTree is a Python package which relies upon a number of packages:

```
- pandas
- numpy 
- geopandas
- geopy
- pyogrio
- matplotlib
- seaborn
- shapely
- pyproj
- networkx
- scikit-learn
- momepy
- mantel
- sortedcontainers
```

The recommended method of installation is with conda or mamba. First, create and activate a new conda environment:
```
conda create -n streamtree_test python=3.10
conda activate streamtree_test
```

Then install the dependencies:
```
# most dependencies installable via conda/mamba
mamba install -c conda-forge -c bioconda pandas numpy geopandas geopy pyogrio matplotlib seaborn shapely pyproj networkx scikit-learn momepy pysam sortedcontainers

# then finally install mantel
pip install mantel
```

Note that if you are running Mac ARM, pysam will need to be installed via pip, so you can use:
```
# most dependencies installable via conda/mamba
mamba install -c conda-forge -c bioconda pandas numpy geopandas geopy pyogrio matplotlib seaborn shapely pyproj networkx scikit-learn momepy sortedcontainers

# then finally install mantel and pysam
pip install mantel pysam
```

autoStreamTree also uses code from 
http://jwcarr.github.io/MantelTest/
ResistanceGA

And methodology described in Kalinowski et al. 2008

### Dependencies <a name="conda"></a>



## autostreamtree <a name="ast"></a>

### Software Description <a name="ast_desc"></a>
autostreamtree is a Python software package providing various analyses aimed at analyzing patterns of genetic differentiation among aquatic stream-dwelling organisms. The intention is to take what was previously a tedious process involving multiple discrete steps and to integrate these all in one place.

Currently, it provides a companion library of functions for calculating various measures of genetic distances among individuals or populations, including model-corrected p-distances (e.g. Jukes-Cantor 1969, Kimura 2-parameter, Tamura-Nei 1993) as well as those based on allele frequencies (e.g. Theta-ST, linearized Fst, Jost's D -- full list of available distance models below). It also includes integrated functions for parsing an input vector shapefile of streams (see below 'Requirements for input shapefiles') for easy calculation of pairwise stream distances between sites, as well as the ordinary or weighted least-squares fitting of reach-wise genetic distances according to the "stream tree" model of Kalinowski et al. (2008). Various plotting functions are also provided for downstream analysis, including looking at patterns of isolation-by-distance. Outputs should also be directly importable into R, with additional outputs with annotated streamtree fitted distances provided for analysis in your GIS suite of choice.

If you use this package for analysis of fitted distances using the streamtree model, please cite the following:
* Kalinowski ST, MH Meeuwig, SR Narum, ML Taper (2008) Stream trees: a statistical method for mapping genetic differences between populations of freshwater organisms to the sections of streams that connect them. Canadian Journal of Fisheries and Aquatic Sciences (65:2752-2760).

#### StreamTree background <a name="ast_background"></a>

<fill in>

### Usage <a name="ast_usage"></a>

#### Options and Help Menu <a name="ast_help"></a>

To view all of the options for DistNet, call the program with the <-h> argument:
```
$ python3 autostreamtree.py -h

Exiting because help menu was called.

DistNet.py

Author: Tyler K Chafin, University of Arkansas
Contact: tkchafin@uark.edu
Description: Computes stream distances and genetic distances for georeferenced DNA sequences, performs tests for isolation-by-distance, and uses a least-squares method to fit distances to stream segments.

	Input file format:
		SampleName	Data	Lat	Long	seq1	[seq2]...[seqn]
		...
		...
	--NOTE: The "DATA" column can be anything- a population/ species identifier
	  (e.g. when used with --pop), or irrelevant data (e.g. GenBank accesson
	  number, if datafile produced by my autoFetcher script)

	Mandatory arguments:
		-s,--shp	: Path to shapefile containing cleaned, contiguous stream reaches
		-i,--input	: Input .tsv file containing sample coordinates and sequences

	General options:
		-o,--out	: Output prefix [default="out"]
		-n,--network	: Provide an already optimized network output from a previous run
			This will be the $out.network file written by autoStreamTree
		--overwrite	: Overwrite an input network (Only relevant with --network)
		-h,--help	: Displays help menu
		-r,--run	: Run which steps? Options: [all, gendist, ibd, streamdist, streamtree]
			ALL		: Run all steps
			GENDIST		: Only calculate genetic distance matrix
			STREAMDIST	: Only compute pairwise stream distances
			DISTANCES	: Only compute GENDIST + STREAMDIST
			IBD		: GENDIST + STREAMDIST + Mantel test
			STREAMTREE	: GENDIST + STREAMDIST + fit StreamTree model
			RUNLOCI	: Run STREAMTREE fitting on each locus
		-p,--pop		: Pool individuals based on column 2 of input file
			NOTE: The location will be taken as the centroid among individual samples
		-g,--geopop		: Pool individuals having identical coordinates
		-c,--clusterpop	: Use DBSCAN algorithm to automatically cluster populations
		--reachid_col	: Attribute name representing primary key in shapefile [default="REACH_ID"]
		--length_col	: Attribute name giving length in kilometers [default="LENGTH_KM"]

	Genetic distance options:
		-d,--dist	: Use which metric of distance? Options are:
			Substitution models (individual-based):
			  PDIST			: Uncorrected p-distances [# Differences / Length]
			  JC69 			: [default] Jukes-Cantor (1969) corrected p-distances
			  K2P			: Kimura 2-parameter distances
			  TN84			: Tajima and Nei's (1984) distance
			  TN93			: Tamura and Nei's (1993) distance
			Frequency models (when using --pop):
			  FST			: Weir and Cockerham's Fst formulation (=THETAst)
			  GST			: Hedrick's (2005) correction of Nei (1987) Gst [=G'st]
			  GSTPRIME		: Meirmans & Hedrick (2011) corrected G'st [=G''st]
			  LINFST		: [default] Rousset's (1997) Fst [=Fst/(1-Fst)]
			  JOST			: Jost's (2008) D
			  NEI72			: Nei's (1972) standard genetic distance
			  NEI83			: Nei and Chesser (1983) Da
			  EUCLID		: Euclidean distance
			  CHORD			: Cavalli-Sforza and Edwards (1967) chord distance
			  --NOTE: Individual-based metrics can also be computed for
		  	          populations. You can set how these are aggregated w/ --pop_agg
			  --NOTE: Multiple loci for PDIST, JC69, K2P, and EUCLID distances
		  	        will be reported using the method defined in --loc_agg
			  --NOTE: TN84 will use empirical base frequencies
		-G,--genmat	: Skip calculation and use the provided labeled .tsv matrix
		--coercemat	: [Boolean] Coerce negative values in input matrix to zero
		--locmatdir	: Directory of per-locus distance matrices
		--het		: [Boolean] Count partial differences [e.g. ind1=T, ind2=W]
		--snp		: [Boolean] Data represent concatenated SNPs
		--msat		: xxx[Boolean] Data represent msat alleles [not yet implemented]
		--global_het	: Estimate Ht using global frequencies (default is averaged over pops)

	DBSCAN options (only when --clusterpop):
		--min_samples	: Minimum samples per cluster [default=1]
		--epsilon		: Maximum distance (in km) within a cluster [default=20]

	Aggregation options:
		-P,--pop_agg	: Define aggregator function for certain genetic distances in pop samples
		-L,--loc_agg	: Define aggregator function for aggregating locus-wise distances
			All of these can take the following options:
			  ARITH		: [default] Use arithmetic mean
			  MEDIAN	: Use median distance
			  HARM		: Use harmonic mean
			  ADJHARM	: Adjusted harmonic mean (see docs)
			  GEOM		: Use geometric mean
			  MIN		: Use minimum distance
			  MAX		: Use maximum distance

	IBD options:
		--perm		: Number of permutations for mantel test [def=1000]
		--and_log	: Also perform IBD steps with log geographic distances

	StreamTree (see Kaliowski et al. 2008) options:
		--iterative	: Prevent negative distances using the iterative approach
		-w,--weight	: Desired weighting for least-squares fitting:
			Options:
			  FM67			: Fitch and Margoliash (1967) [w = 1/D^2]
			  BEYER74		: Beyer et al. (1974) weights [w = 1/D]
			  CSE67			: [default] Cavalli-Sforza and Edwards (1967) [w = 1]
```

#### Requirements for input shapefiles <a name="ast_input"></a>

I highly recommend using the existing global stream datasets provided by the [HydroLab group](https://wp.geog.mcgill.ca/hydrolab/) at McGill University, specifically the [HydroAtlas](https://www.hydrosheds.org/page/hydroatlas) or [free-flowing rivers dataset](https://wp.geog.mcgill.ca/hydrolab/free-flowing-rivers/) as these are already properly formatted for use, and the additional variables included will make downstream analysis very easy. Because of their size, I would recommend clipping them to the relevant scale first (e.g. the drainage encompassing all of your samples).

Note that a valid path is required between all sites in order to calculate pairwise stream distances. Thus, if you are analyzing data from multiple drainages which only share an oceanic connection, you will need to augment the shapefile. For example this could be accomplished by adding a vector representing the coastline to create an artificial connection among drainages.

If for some reason you cannot use the HydroRIVERS dataset, you will need to do some things first before loading your shapefile into autoStreamTree. First, you will need to include two variables in the attribute table of your shapefile: 1) REACH_ID (case sensitive) must provide a unique identifier to each stream reach; and 2) LENGTH_KM should give the length of each segment. Next, because sometime large stream layers will have small gaps in between streams, you will need to span any small gaps between streams which should be contiguous, and also dissolve any lines that overlap with one another so that any given section of river is represented by a single line. I will provide a tutorial for doing this in ArcMAP later, but for now there are some scripts in our complementary package that can help with these steps using the ArcPy API: https://github.com/stevemussmann/StreamTree_arcpy. Note that this package will also help you in running the original Stream Tree package on Windows, if you want to do so.

#### Input file format <a name="ast_table"></a>

Coordinates do not need to exactly match nodes in the input shapefile, as points will be 'snapped' to network nodes after parsing. autoStreamTree will output both a table($out.snapDistances.txt) and a histogram plot ($out.snapDistances.pdf) showing distances in kilometers that samples or populations had to be snapped:

![](https://raw.githubusercontent.com/tkchafin/autostreamtree/master/examples/plots/example.snapDistances.png)

#### Outputs <a name="ast_outputs"></a>

The first thing autoStreamTree will do upon reading your input shapefile is to calculate a minimally reduced sub-network which collapses the input river network into continuous reaches (="edges"), with nodes either representing sample localities or junctions. Because the full river network will likely contain many branches and contiguous reaches which do not contain samples, these are removed to speed up computation. The underlying metadata will be preserved, and the final output will consist of an annotated shapefile containing an EDGE_ID attribute which tells you how reaches were dissolved into contiguous edges in the graph, and a FittedD attribute giving the least-squares optimized distances.

The reduced sub-network will be plotted for you in a file called $OUT.subGraph.pdf:
![](https://raw.githubusercontent.com/tkchafin/autostreamtree/master/examples/plots/example.subGraph.png)

Here, the total cumulative stream length (in km) is plotted along edges (NOTE: Any natural curvature in the river is not preserved in this plot), with sample sites as blue dots and junctions as black dots. A geographically accurate representation, coloring individual streams to designate different dissolved edges, will be provided as $out.streamsByEdgeID.pdf:
![](https://raw.githubusercontent.com/tkchafin/autostreamtree/master/examples/plots/example.networkByEdgeID.png)

After fitting genetic distances, autoStreamTree will create several other outputs. First, a table called $out.reachToEdgeTable.txt will give a tab-delimited map of how REACH_ID attributes were dissolved into contiguous edges. Second, a tabular and graphical representation of how fitted pairwise distances compare to the raw calculates (or user-provided) pairwise distances: $out.obsVersusFittedD.txt and $out.obsVersusFittedD.pdf
![](https://raw.githubusercontent.com/tkchafin/autostreamtree/master/examples/plots/example.obsByFittedD.png)

Finally, the fitted distances per stream edge will be output both as an added column to the original shapefile attribute table ($out.streamTree.shp and $out.streamTree.txt), and also as a plot showing how distances compare across all streams:
![](https://raw.githubusercontent.com/tkchafin/autostreamtree/master/examples/plots/example.networkByStreamTree.png)

#### Genetic distance models <a name="gen"></a>

The currently recommended way to run autoStreamTree is to provide a labelled matrix of pairwise genetic distances, as there are many available packages for calculating these. This input matrix is provided using the --genmat argument, and should be tab-delimited with both column and row labels matching your population or individual identifiers.

Built-in genetic distance calculations are currently in beta, meaning they are provided conditionally and still require more extensive external testing. They appear to be by-and-large functioning as written, but I would like to have a much more comprehensive test of whether or not my implementation of each statistic produces identical (or very similar) results to those of other packages such as Arlequin or adegenet. If you use autoStreamTree for your research and have the ability to directly compare your genetic distance matrix with those of other packages, please do so and feel free to let me know how they compare. Please note that I also offer the ability to directly import a genetic distance matrix that has been externally calculated, skipping this step altogether.

Currently I provide options for individual-based distances (p-distance and various substitution-model corrected distances), and population-based distances which use allele frequency data (such as Fst):

```
	Genetic distance options:
		-d,--dist	: Use which metric of distance? Options are:
			Substitution models (individual-based):
			  PDIST			: Uncorrected p-distances [# Differences / Length]
			  JC69 			: [default] Jukes-Cantor (1969) corrected p-distances
			  K2P			: Kimura 2-parameter distances
			  TN84			: Tajima and Nei's (1984) distance
			  TN93			: Tamura and Nei's (1993) distance
			Frequency models (when using --pop):
			  FST			: Weir and Cockerham's Fst formulation (theta)
			  GST			: Hedrick's (2005) correction of Nei's (1987) Gst [=G'st]
			  GSTPRIME		: Meirmans & Hedrick (2011) corrected G'st [=G''st]
			  LINFST		: [default] Rousset's (1997) linearized Fst [=Fst/(1-Fst)]
			  JOST			: Jost's (2008) D
			  LINJOST		: 1/1-D, where D=Jost's (2008) D
			  NEI72			: Nei's (1972) standard genetic distance
			  NEI83			: Nei and Chesser (1983) Da
			  EUCLID		: Euclidean distance
			  CHORD			: Cavalli-Sforza and Edwards (1967) chord distance
```

Optionally, the user can also opt to aggregate individual-based distance measures (when using a priori population assignments or the --geopop option). This can be provided using the --pop_agg argument, with any of the following options available:

```
	Aggregation options:
		--pop_agg	: Define aggregator function for certain genetic distances w/ --pops:
			All of these can take the following options:
			  ARITH		: [default] Use arithmetic mean
			  MEDIAN	: Use median distance
			  HARM		: Use harmonic mean
			  ADJHARM	: Adjusted harmonic mean (see docs)
			  GEOM		: Use geometric mean
			  MIN		: Use minimum distance
			  MAX		: Use maximum distance
```

For datasets containing multiple non-concatenated loci, note that individual-based distances (e.g. PDIST or JC69) will also need to be aggregated among loci within each pairwise calculation. Any of the above options can again be used here, provided using the --loc_agg argument.

#### Defining populations <a name="pops"></a>
There are currently three ways in which you can define populations for population-wise analysis. The first (specified using --pop) assumes that the 2nd column in the input file contains population identifiers. These can take any form (e.g., integer or string). The second (--geopop) will group any samples into populations which "snap" to the same stream node (see below).

A third option (--clusterpop) will automatically cluster geographically similar individuals using the DBSCAN algorithm in scikit-learn, using great-circle geographic distances (i.e., this is not informed by stream distances calculated as a part of some workflows). Two relevant options are provided for manipulating the DBSCAN results:
```
DBSCAN options (only when --clusterpop):
	--min_samples	: Minimum samples per cluster [default=1]
	--epsilon	: Maximum distance (in km) within a cluster [default=20]
```

If using population labels, whether provided in the input file (--pop/--geopop) or calculating using DBSCAN (--clusterpop), autoStreamTree will output a plot showing cluster membership in 2-D space called $OUT.clusteredPoints.pdf:

![](https://raw.githubusercontent.com/tkchafin/autostreamtree/master/examples/plots/example.clusteredPoints.png)

In this example, DBSCAN was used (hence population IDs are formatted as "DB_"#). Population centroids, which are ultimately used to "snap" populations to the stream network are shown with an "x". Note that this means that the population will only be represented by a single point on the network!


### Example workflows <a name="ast_workflow"></a>

#### Fitting single-marker distances <a name="ast_example1"></a>

#### Working with large SNP datasets <a name="ast_example2"></a>

#### Microsatellites <a name="ast_example3"></a>

### Runtimes and benchmarking <a name="ast_benchmark"></a>

### References <a name="ast_refs"></a>
#### Citations for DistNet methods
Below is a full list of citations for the various methods used in DistNet. Apologies to anyone I missed - feel free to let me know if you notice any discrepancies.
* Beyer WM, Stein M, Smith T, Ulam S. 1974. A molecular sequence metric and evolutionary trees. Mathematical Biosciences. 19: 9-25.
* Cavalli-Sforza LL, Edwards AWF. 1967. Phylogenetic analysis: model and estimation procedures. American Journal of Human Genetics. 19: 233-257.
* Ester M, Kriegel HP, Sander J, Xu X. 1996. A density-based algorithm for discovering  clusters in large spatial databases with noise. IN: Simoudis E, Han J, Fayyad UM. (eds.). Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231.
* Felsenstein J. 2004. Inferring Phylogenies: Chapter 11. Sunderland: Sinauer.
* Fitch WM, Margloiash E. 1967. Construction of phylogenetic trees. Science. 155: 279-84.
* Hagberg A, Swart P, S Chult D. 2008. Exploring network structure, dynamics, and function using NetworkX. Los Alamos National Lab.(LANL), Los Alamos, NM
* Hedrick PW. 2005. A standardized genetic differentiation measure. Evolution. 59: 1633–1638
* Jordahl K. 2014. GeoPandas: Python tools for geographic data. URL: https://github.com/geopandas/geopandas.
* Jost L. 2008. Gst and its relatives do not measure differentiation. Molecular Ecology. 17: 4015-4026.
* Jukes TH, Cantor CR. 1969. Evolution of protein molecules. New York: Academic Press.
* Kalinowski ST, MH Meeuwig, SR Narum, ML Taper (2008) Stream trees: a statistical method for mapping genetic differences between populations of freshwater organisms to the sections of streams that connect them. Canadian Journal of Fisheries and Aquatic Sciences (65:2752-2760)
* Kimura M. 1980. A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. Journal of Molecular Evolution. 16(2): 111-120.
* Mantel N. 1967. The detection of disease clustering and a generalized regression approach. Cancer Research 27(2): 209-220.
* Meirmans PG, Hedrick PW. 2011. Assessing population structure: Fst and related measures. Molecular Ecology Resources. 11: 5-18.
* Nei M. 1972. Genetic distance between populations. American Naturalist. 106: 283-292.
* Nei M. 1987. Molecular Evolutionary Genetics. Columbia University Press, New York
* Nei M, Chesser RK. 1983. Estimation of fixation indices and gene diversities. Annals of Human Genetics 47(3): 253-259.
* Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O, Blondel M, Prettenhofer P, Weiss R, Dubourg V, Vanderplas J. 2011. Scikit-learn: Machine learning in Python. The Journal of machine Learning research. 1(12):2825-30
* Rossmann LA. DFLOW User's Manual. U.S. Environmental Protection Agency.[For description of zero-adjusted harmonic mean]
* Rousset F. 1997. Genetic differentiation and estimation of gene flow from F-statistics under isolation by distance. Genetics. 145: 1219-28.
* Tajima F, Nei M. 1984. Estimation of evolutionary distance between nucleotide sequences. Molecular Biology and Evolution 1:269-285
* Tamura K, Nei M. 1993. Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. Molecular Biology and Evolution. 10(3):512-526.
* Weir BS, Cockerham CC. 1984. Estimating F-statistics for the analysis of population structure. Evolution. 38: 1358-1370.

#### Other reading
Here are some recommended readings and resources:
* Comte L, Olden JD. 2018. Fish dispersal in flowing waters: A synthesis of movement- and genetic-based studies. Fish and Fisheries. 19(6): 1063-1077.
* Comte L, Olden JD. 2018. Evidence for dispersal syndromes in freshwater fishes. Proceedings Royal Society: B. 285(1871):  
* Grill, G., Lehner, B., Thieme, M. et al. 2019. Mapping the world’s free-flowing rivers. Nature. 569:215–221.
* Linke S, Lehner B, Ouellet Dallaire C. et al. 2019. Global hydro-environmental sub-basin and river reach characteristics at high spatial resolution. Sci Data 6, 283
* Meffe GK, Vrijenhoek RC. 1988. Conservation genetics in the management of desert fishes. Conservation Biology. 2(2):157-69.
* Meirmans PG. 2012. The trouble with isolation by distance. Molecular Ecology 21(12): 2839-46.
* Sere M, Thevenon S, Belem AMG, De Meeus T. 2017. Comparison of different genetic distances to test isolation by distance between populations. 2017. 119(2):55-63.
* Thomaz AT, Christie MR, Knowles LL. 2016. The architecture of river networks can drive the evolutionary dynamics of aquatic populations. Evolution. 70(3): 731-739.
* Tonkin JD, Altermatt F, Finn DS, Heino J, Olden JD, Pauls SU, Lytle DA. 2017. The role of dispersal in river network metacommunities: Patterns, processes, and pathways. Freshwater Biology. 61(1): 141-163.
* Wright S. 1965. Isolation by distance. Genetics. 28: 114-138.

## Scripts and Tools

#### Re-plotting StreamTree outputs

The default $out.streamdByFittedD plot may not be exactly what you wanted. To prevent cluttering the help menu of the main program too much, we've provided a separate script for loading up autoStreamTree outputs to re-make the plot, which has some added options for customization: scripts/plotStreamTree.py

'''
$ python3 ./scripts/plotStreamTree -h
plotStreamTree.py

Author: Tyler K Chafin, University of Arkansas
Contact: tkchafin@uark.edu
Description: Script for re-plotting StreamTree results after running DistNet

		-p,--prefix	: Prefix for autoStreamTree output
		-m,--min	: Minimum genetic distance
		-M,--max	: Maximum genetic distance
		-c,--cmap	: Colormap (any valit matplotlib cmap value)
			see: https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html
		-o,--out	: Output prefix (if not overwriting original)
		-h,--help	: Displays this help menu)
'''

For example, to re-plot values with a distance ceiling of 0.2 and a viridis color scale:
'''
python3 ./scripts/plotStreamTree.py -p out2 -m 0.0 -M 0.2 -c "viridis"
'''

#### Clustering populations using any distance matrix

For example, if you wanted to cluster individuals using their stream distances, I've provided a script called clusterPopsDB.py which will use a DBSCAN clustering algorithm to output a tab-delimited population map given any arbitrary distance matrix:
'''
$ python3 ./scripts/clusterPopsDB.py -h
'''
