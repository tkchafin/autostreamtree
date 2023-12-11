# autostreamtree

`autoStreamTree` is a Python package meant for analysing genomic variant (e.g., SNPs, microhaplotypes) differentiation in geospatial networks such as riverscapes. 

An input geospatial network (as shapefile, geodatabase, or .gpkg) is first used to compute a minimally representative graph structure for all least-cost paths between your sampling locations. Pairwise genomic distances (e.g., Fst) are then 'fitted' to this graph by calculating the contribution of each graph segment, using the least-squares algorithm from Kalinowsky et al. (2008).

`autoStreamTree` automates what was previously a manual process, while also integrating new functionality such as per-locus or per-site fitted distances, which can be used for downstream outlier analysis.

### Table of Contents:
1. [Installation](#installation)
2. [Usage](#ast)
    1. [Interface](#ast_cli)
    2. [Inputs](#ast_inputs)
    3. [Genetic distances](#gendist)
    4. [Outputs](#ast_outputs)
3. [Example](#tutorial)
4. [Downstream analysis](#downstream)
    1. [Utility scripts](#utils)
    2. [Incorporating into workflows](#workflow)
5. [References](#ast_refs)
6. [Contributing Guidelines](#contrib)


## 1. Installation <a name="installation"></a>

### Dependencies

`autoStreamTree` is a Python package which relies upon a number of dependencies:

```
- pandas
- numpy 
- geopandas
- pyogrio
- matplotlib
- seaborn
- pyproj
- networkx
- scikit-learn
- momepy
- mantel
- pysam
- sortedcontainers
```

`autoStreamTree` also uses methodology described in [Kalinowski et al. 2008](https://www.montana.edu/kalinowski/documents/2008_Stream_Trees%20_CJFAS.pdf), if you use this software you should cite: 
```
Kalinowski ST, MH Meeuwig, SR Narum, ML Taper (2008) Stream trees: a statistical method for mapping genetic differences between populations of freshwater organisms to the sections of streams that connect them. Canadian Journal of Fisheries and Aquatic Sciences (65:2752-2760)
```

### Installation with GitHub

Simple installation of the development version may be done simply as:
```
git clone https://github.com/tkchafin/autostreamtree.git

cd autostreamtree

pip install . 
```

The command-line interface is then accesible as `autostreamtree`.

### Installation with conda/mamba

The recommended method of installation is with conda or mamba. First, create \
and activate a new conda environment:
```
conda create -n streamtree python=3.10
conda activate streamtree
```

Then install the package:
```
mamba install -c ecoevoinfo -c conda-forge -c bioconda autostreamtree
```

If you require the development branch, you can then install from GitHub as \
well:

```
git clone https://github.com/tkchafin/autostreamtree.git
cd autostreamtree 
pip install .
```

You should now be ready to go!

#### Note for Apple silicon (ARM) Macs

Currently there are some dependencies which can't be easily installed on Macs 
running the Apple ARM CPU architectures. If you are using an Apple Silicon (M1/M2) Mac there are some extra steps to run x86_64 (Intel-based) applications, which you can do via Rosetta 2. Here are the steps to set up and launch a terminal in Rosetta mode:

1. **Install Rosetta**:
   - Open the Terminal app (found in `Applications` > `Utilities`).
   - Enter the following command and press Return:

     ```
     /usr/sbin/softwareupdate --install-rosetta --agree-to-license
     ```

   - If Rosetta is already installed, this command will have no effect. If it's not installed, you'll be prompted to agree to the license terms, and the installation will proceed.

2. **Find Terminal in Finder**:
   - Open Finder.
   - Navigate to the `Applications` folder, then `Utilities`.
   - Find `Terminal`.

3. **Duplicate Terminal**:
   - Right-click (or Ctrl-click) on `Terminal`.
   - Select `Duplicate` to create a copy of the Terminal application. This ensures you have a separate instance for Rosetta without altering your default Terminal.

4. **Enable Rosetta for the Duplicated Terminal**:
   - Right-click (or Ctrl-click) on the duplicated Terminal.
   - Select `Get Info`.
   - In the `Info` window, check the box `Open using Rosetta`.

5. **Launch the Rosetta Terminal**:
   - Open the duplicated Terminal application. This instance of Terminal will run under Rosetta, allowing you to run x86_64 architecture applications.


Once this is done, you can create a conda environment:

```
conda create -n streamtree
conda activate streamtree 
conda config --env --set subdir osx-64
conda install python=3.10 mamba
```

Then install `autostreamtree`:

```
mamba install -c ecoevoinfo -c conda-forge -c bioconda autostreamtree
```

## 2. Usage <a name="ast"></a>

### Command-line interface <a name="ast_cli"></a>

After installation, you can view the up-to-date help menu by calling the script with `autostreamtree -h`:

```

autostreamtree

Author: Tyler K Chafin, Biomathematics and Statistics Scotland
Description: Methods for analysing genetic distances in networks.

Mandatory arguments:
    -s, --shp       : Path to shapefile containing cleaned, contiguous stream reaches
                       (can also support geodatabase or GPKG files)
    -i, --input     : Input .tsv file containing sample coordinates
    -v, --vcf       : Input VCF file containing genotypes

General options:
    -o, --out       : Output prefix [default="out"]
    -O, --gdf_out   : Output driver for annotated geodataframe (options "SHP", "GPKG", "GDB")
    -C, --concat    : Concatenate all SNPs ("all"), by locus ("loc"), or not at all ("none")
    -n, --network   : Provide an already optimized network output from a previous run
    --overwrite     : Overwrite an input network (Only relevant with --network)
    -h, --help      : Displays help menu
    -r, --run       : Run which steps? Options:
        ALL         : Run all steps
        GENDIST     : Only calculate genetic distance matrix
        STREAMDIST  : Only compute pairwise stream distances
        DISTANCES   : Only compute GENDIST + STREAMDIST
        IBD         : GENDIST + STREAMDIST + Mantel test
        STREAMTREE  : GENDIST + STREAMDIST + fit StreamTree model
        RUNLOCI     : Run STREAMTREE fitting on each locus
    -p, --pop       : Pool individuals based on an input population map tsv file
        NOTE: The location will be taken as the centroid among individual samples
    -g, --geopop    : Pool individuals having identical coordinates
    -c, --clusterpop: Use DBSCAN algorithm to automatically cluster populations
    --reachid_col   : Attribute name representing primary key in shapefile [default="HYRIV_ID"]
    --length_col    : Attribute name giving length in kilometers [default="LENGTH_KM"]
    --seed          : Seed for RNG

Genetic distance options:
    -d, --dist      : Use which metric of distance? Options:
        Individual-based:
          PDIST     : Uncorrected p-distances [# Differences / Length]
        Frequency models (when using --pop):
          FST       : Weir and Cockerham's Fst formulation (=THETAst)
          LINFST    : [default] Rousset's (1997) Fst [=Fst/(1-Fst)]
          JOST      : Jost's (2008) D
          CHORD     : Cavalli-Sforza and Edwards (1967) chord distance
          --NOTE: Individual-based metrics can also be computed for
                    populations. You can set how these are aggregated w/ --pop_agg
          --NOTE: Multiple loci for PDIST
                  will be reported using the method defined in --loc_agg
    -G, --genmat    : Skip calculation and use the provided labeled .tsv matrix
    --coercemat     : [Boolean] Coerce negative values in input matrix to zero
    --het           : [Boolean] Count partial differences [e.g. ind1=T, ind2=W]
    --global_het    : Estimate Ht using global frequencies (default is averaged over pops)

DBSCAN options (only when --clusterpop):
    --min_samples   : Minimum samples per cluster [default=1]
    --epsilon       : Maximum distance (in km) within a cluster [default=20]

Aggregation options:
    -P, --pop_agg   : Define aggregator function for certain genetic distances in pop samples
    -L, --loc_agg   : Define aggregator function for aggregating locus-wise distances
        All of these can take the following options:
          ARITH     : [default] Use arithmetic mean
          MEDIAN    : Use median distance
          HARM      : Use harmonic mean
          ADJHARM   : Adjusted harmonic mean (see docs)
          GEOM      : Use geometric mean
          MIN       : Use minimum distance
          MAX       : Use maximum distance

IBD options:
    --perm          : Number of permutations for mantel test [def=1000]
    --and_log       : Also perform IBD steps with log geographic distances

StreamTree options (see Kalinowski et al. 2008) :
    -w, --weight    : Desired weighting for least-squares fitting:
        Options:
          FM67      : Fitch and Margoliash (1967) [w = 1/D^2]
          BEYER74   : Beyer et al. (1974) weights [w = 1/D]
          CSE67     : [default] Cavalli-Sforza & Edwards (1967) [w = 1]

```

The main workflow control argument is `-r/--run`, which tells which steps of the pipeline you want to run. The options for this are:
- -r GENDIST: Only calculates a pairwise genetic distance matric for input samples 
- -r STREAMDIST: Only calculates a pairwise matrix of hydrologic (stream) distances
- -r DISTANCES: Computes both of the above 
- -r IBD: Computes both distance matrices and performes a simple Mantel test for isolation-by-distance
- -r STREAMTREE: Computes distance matrices and fits the StreamTree model 
- -r ALL: All of the above 
- -r RUNLOCI: All of the above, but fitting StreamTree models to each individual locus 

Some of these are included only for convenience -- for example, there are much more comprehensive packages out there for computing genetic distance matrices, and I don't recommend the Mantel test as the sole means of testing for isolation-by-distance -- and in most cases I think `-r ALL` or `-r STREAMTREE` will be what you want. 

Many of the runtime options relate to specific components of the required inputs, which I will discuss below.

### Input files <a name="ast_inputs"></a>

Depending on the workflow you are running, some of these files may not be required (for example, computing genetic distances only will not require an input geodatabase). 

#### Genotype data 

The intended input format for the genotype data is the widely used VCF format. Some examples are provided in the example `data/` directory, which were produced by the program ipyrad:
```
##fileformat=VCFv4.0
##fileDate=2021/03/21
##source=ipyrad_v.0.9.68
##reference=pseudo-reference (most common base at site)
##phasing=unphased
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=CATG,Number=1,Type=String,Description="Base Counts (CATG)">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  63burk03        63cburk01       >
RAD_0   28      loc0_pos27      G       A       13      PASS    .       GT:DP:CATG      0/0:81:0,0,0,81 >
RAD_0   76      loc0_pos75      C       T       13      PASS    .       GT:DP:CATG      0/0:81:81,0,0,0 >
...
...
```

Parsing is faster with compressed and indexed files (and will be more space efficient), so I recommend keeping these compressed with bgzip and indexed with tabix:
```
# compress
bgzip test_sub100.vcf 

# index
tabix test_sub100.vcf.gz
```

This file should be provided using `-v/--vcf`.
#### Individual metadata

Individual metadata will be needed in the form of coordinates, which can be passed as a simple tab-delimited text file with header (see example in `data/test.coords`):
```
sample  lat     long
63burk03        26.925414       90.396139
63cburk01       26.92383        90.39815
63cburk02       26.92383        90.39815
63cdikc01       27.2676 90.04778
63cdikc02       27.2676 90.04778
63cdikc03       27.2676 90.04778
63cdolk01       26.88353        90.34272
63cdolk03       26.88353        90.34272
63cdort03       26.86776        89.37437
...
...
```

For population-level analysis, there are currently three ways in which you can define populations for population-wise analysis. One option (`-g/--geopop`) will group any samples into populations which "snap" to the same stream node (see below).

A third option (`-c,--clusterpop`) will automatically cluster geographically similar individuals using the DBSCAN algorithm in scikit-learn, using great-circle geographic distances (i.e., this is not informed by stream distances calculated as a part of some workflows). Two relevant options are provided for manipulating the DBSCAN results:
```
	--min_samples	: Minimum samples per cluster [default=1]
	--epsilon	    : Maximum distance (in km) within a cluster [default=20]
```

or using an input file providing population assignments for each individual, again as a simple tab-delimited text file (see example `data/test.popmap`), passed via `-p/--pops`:
```
63burk03        burk
63cburk01       burk
63cburk02       burk
63dakp01        dakp
63dakp02        dakp
63dakp03        dakp
63dakp04        dakp
63dakp05        dakp
...
...
```
#### Geodatabase 
The input stream network can be provided as a shapefile, geodatabase, or GPKG file, all passed uing the `-s/--shp` option. There are a number of requirements for this file in order for the result to create a valid network. I highly recommend using the existing global stream datasets provided by the [HydroLab group](https://wp.geog.mcgill.ca/hydrolab/) at McGill University, specifically the [HydroAtlas](https://www.hydrosheds.org/page/hydroatlas) or [free-flowing rivers dataset](https://wp.geog.mcgill.ca/hydrolab/free-flowing-rivers/) as these are already properly formatted for use, and the additional variables included will make downstream analysis very easy. 

If for some reason you cannot use the HydroRIVERS dataset, you will need to do some things first before loading your shapefile into autoStreamTree. First, you will need to include two variables in the attribute table of your shapefile: 1) `REACH_ID` (case sensitive) must provide a unique identifier to each stream reach; and 2) `LENGTH_KM` should give the length of each segment. Next, because sometime large stream layers will have small gaps in between streams, you will need to span any small gaps between streams which should be contiguous, and also dissolve any lines that overlap with one another so that any given section of river is represented by a single line. I will provide a tutorial for doing this in `ArcMAP` later, but for now there are some scripts in our complementary package that can help with these steps using the `ArcPy` API: https://github.com/stevemussmann/StreamTree_arcpy. Note that this package will also help you in running the original Stream Tree package on Windows, if you want to do so.

Note that a valid path is required between all sites in order to calculate pairwise stream distances. Thus, if you are analyzing data from multiple drainages which only share an oceanic connection, you will need to augment the shapefile. For example this could be accomplished by adding a vector representing the coastline to create an artificial connection among drainages.

### Genetic distances <a name="gendist"></a>

There are a number of existing (and more comprehensive) packages out there for computing pairwise genetic distance matrices. If you already have this (either at individual or population level), you can provide a labelled matrix as a tab-delimited text file with `--genmat`. An example of how this should be formatted may be found in `data/test.popGenDistMat.txt`.

For convenience, a number of options are built-in, which can be selected using the `-d,--dist` argument:
```
    -d, --dist      : Use which metric of distance? Options:
        Individual-based:
          PDIST     : Uncorrected p-distances [# Differences / Length]
        Frequency models (when using --pop):
          FST       : Weir and Cockerham's Fst formulation (=THETAst)
          LINFST    : [default] Rousset's (1997) Fst [=Fst/(1-Fst)]
          JOST      : Jost's (2008) D
          CHORD     : Cavalli-Sforza and Edwards (1967) chord distance
          --NOTE: Individual-based metrics can also be computed for
                    populations. You can set how these are aggregated w/ --pop_agg
          --NOTE: Multiple loci for PDIST, and JC69 distances
                  will be reported using the method defined in --loc_agg
```
For most use cases, I would suggest `-d LINFST` which will compute a "linearised" version of Weir and Cockerham's Theta-ST. 

Optionally, the user can also opt to aggregate individual-based distance measures, either those provided (p-distances) or from an input matrix that is only at the individual level. This can be provided using the `--pop_agg` argument, with any of the following options available:

```
Aggregation options:
    -P, --pop_agg   : Define aggregator function for certain genetic distances in pop samples
    -L, --loc_agg   : Define aggregator function for aggregating locus-wise distances
        All of these can take the following options:
          ARITH     : [default] Use arithmetic mean
          MEDIAN    : Use median distance
          HARM      : Use harmonic mean
          ADJHARM   : Adjusted harmonic mean (see docs)
          GEOM      : Use geometric mean
          MIN       : Use minimum distance
          MAX       : Use maximum distance
```

Another useful feature is the ability to concatenate the input variant data -- this can be done either globally (for example if you want to only compute a global p-distance) with `-C all`, or if you have phased data used to form "pseudo" microhaplotypes using `-C loc`, which will group variants using the CHROM field in the input VCF. Note that if using a program such as ipyrad (e.g., for RADseq data), the entries in this field will represent 'independent' loci -- if you have data aligned to a genome and want to define loci in some way you will need to insert this information into the CHROM field as a pre-processing step. 

### Outputs <a name="ast_outputs"></a>

If running the full workflow, the first thing `autoStreamTree` will do upon reading your input geodatabase is to calculate a minimally reduced sub-network which collapses the input river network into continuous reaches (="edges"), with nodes either representing sample localities or junctions. Because the full river network will likely contain many branches and contiguous reaches which do not contain samples, these are removed to speed up computation. The underlying metadata will be preserved, and the final output will consist of an annotated shapefile containing an `EDGE_ID` attribute which tells you how reaches were dissolved into contiguous edges in the graph, and a FittedD attribute giving the least-squares optimized distances.

The reduced sub-network will be plotted for you in a file called *out*.subGraph.pdf.

Here, the total cumulative stream length (in km) is plotted along edges (NOTE: Any natural curvature in the river is not preserved in this plot), with sample sites as blue dots and junctions as black dots. A geographically accurate representation, coloring individual streams to designate different dissolved edges, will be provided as *out*.streamsByEdgeID.pdf.

After fitting genetic distances, `autoStreamTree` will create several other outputs. First, a table called *out*.reachToEdgeTable.txt will give a tab-delimited map of how `REACH_ID` attributes were dissolved into contiguous edges. Second, a tabular and graphical representation of how fitted pairwise distances compare to the raw calculates (or user-provided) pairwise distances: *out*.obsVersusFittedD.txt and *out*.obsVersusFittedD.pdf

Finally, the fitted distances per stream edge will be output both as an added column to the original shapefile attribute table (*out*.streamTree.shp and *out*.streamTree.txt), and also as a plot showing how distances compare across all streams.

## 3. Example <a name="ast_example"></a>

If you have just installed `autoStreamTree` and are running it for the first time, I recommend you first run the example analysis using the provided files in the `autostreamtree/data` directory. These include all necessary inputs, and a geodatabase. To run the full workflow on the example data, simply use:

```
# Change directories to the autoStreamTree repository
cd autostreamtree 

# run the analysis
autostreamtree -s autostreamtree/data/test.shp -i autostreamtree/data/test.coords -v autostreamtree/data/test_sub100.vcf.gz -p autostreamtree/data/test.popmap -r ALL --reachid_col "HYRIV_ID" --length_col "LENGTH_KM" -o test
```

This will produce a number of output text files and plots using the prefix provided with `-o`.


## 4. Downstream analysis <a name="downstream"></a>

`autoStreamTree` is designed to aid in the downstream analysis of genetic differentiation, adaptation, and ecological aspects of aquatic species, although it is versatile enough to be applicable to any dataset structured as a network. 

A key use case involves the examination of locus-wise genetic distances across network segments. This is efficiently done by importing the out.streamTree.tsv file, which is an output of autoStreamTree using the `-r LOC` option, into R or a similar analytical platform. Moreover, if autoStreamTree is run with an annotated geodatabase, such as HydroRIVERS, the resultant fitted distances are enriched with additional attributes. These attributes are particularly useful in conducting in-depth analyses to understand how environmental covariates influence genetic patterns, for instance, through methods like redundancy analysis. This multifaceted approach makes autoStreamTree a robust tool for exploring and interpreting complex genetic and ecological data in aquatic species.

### Utility scripts <a name="utils"></a>

`autoStreamTree` is currently distributed with two utility scripts: `networkDimensions.py` and `streeToDendrogram.py`. We expect this list to grow, and if you need help writing any specific utilities using `autoStreamTree` outputs, please don't hesitate to get in touch! 

As with the `autoStreamTree` CLI, when installed via `pip` or `conda/mamba`, these can be accessed via front-ends `networkDimensions` and `streeToDendrogram`, both having help menus displayed by calling with `-h`

#### networkDimensions

`networkDimensions` is a simple script to capture the spatial dimensions of your minimized network, including the total length, and the area of a the minimally bounding convex hull. 

The options are: 
```
(autostreamtree) tyler@Tylers-MacBook-Pro-2 autostreamtree % networkDimensions -h
usage: networkDimensions [-h] network_file

Calculate total lineage length and bounding box area of a network.

positional arguments:
  network_file  Path to the pickled Networkx graph file.

options:
  -h, --help    show this help message and exit
```

Simply provide the output `.network` file created by `autoStreamTree`.

#### streeToDendrogram

`streeToDendrogram` makes it easier to visualise the fitted distances as a dendrogram (rather than a colored spatial network). The options are:

```
(autostreamtree) tyler@Tylers-MacBook-Pro-2 autostreamtree % streeToDendrogram -h
usage: streeToDendrogram [-h] --shp SHP --points POINTS [--pop POP] [--edge_id EDGE_ID]
                         [--dist DIST] [--out OUT]

Plot a stream tree with fittedD attribute as dendrogram

options:
  -h, --help         show this help message and exit
  --shp SHP          Network as shapefile
  --points POINTS    Sample coordinates
  --pop POP          Population map file (optional)
  --edge_id EDGE_ID  Edge ID attribute.
  --dist DIST        Dist attribute.
  --out OUT          Output prefix
```

### Incorporating into workflows <a name="workflow"></a>
All utilities from `autoStreamTree` are intended to be easy to run within an automated pipeline. An example of building these into a [`SnakeMake`](https://snakemake.github.io) workflow can be found in the [Open Science Framework](https://osf.io/4uxfj/) repository (doi: 10.17605/OSF.IO/4UXFJ).

For example, in the `run_spd.smk` file, you will find an example of fitting distances for regional subsets of samples for an empirical ddRADseq dataset, and then performing downstream environmental resistance modelling on each using [`ResistNet`](https://github.com/tkchafin/resistnet)

We hope to provide some more `SnakeMake` example pipelines soon!
## 5. References <a name="ast_refs"></a>
### Citations for autoStreamTree methods
Below is a full list of citations for the various methods used. Apologies to anyone I missed - please let me know if you notice any discrepancies.

* Beyer WM, Stein M, Smith T, Ulam S. 1974. A molecular sequence metric and evolutionary trees. Mathematical Biosciences. 19: 9-25.
* Cavalli-Sforza LL, Edwards AWF. 1967. Phylogenetic analysis: model and estimation procedures. American Journal of Human Genetics. 19: 233-257.
* Ester M, Kriegel HP, Sander J, Xu X. 1996. A density-based algorithm for discovering  clusters in large spatial databases with noise. IN: Simoudis E, Han J, Fayyad UM. (eds.). Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231.
* Felsenstein J. 2004. Inferring Phylogenies: Chapter 11. Sunderland: Sinauer.
* Fitch WM, Margloiash E. 1967. Construction of phylogenetic trees. Science. 155: 279-84.
* Hagberg A, Swart P, S Chult D. 2008. Exploring network structure, dynamics, and function using NetworkX. Los Alamos National Lab.(LANL), Los Alamos, NM
* Hedrick PW. 2005. A standardized genetic differentiation measure. Evolution. 59: 1633–1638
* Jordahl K. 2014. GeoPandas: Python tools for geographic data. URL: https://github.com/geopandas/geopandas.
* Jost L. 2008. Gst and its relatives do not measure differentiation. Molecular Ecology. 17: 4015-4026.
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

### Other reading
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

## 5. Contributing <a name="contrib"></a>

If you are interested in contributing to autoStreamTree, first of all thank you! Any form of contribution is always welcome, be that bug reports, code review, bug fixes, or even adding new features! 

If you encounter a problem in using or installing autoStreamTree, please post to the [GitHub Issues](https://github.com/tkchafin/autostreamtree/issues) page. If you do so, please be sure to share the full steps which led to your problem, and the full output (including any error messages). In some instances it may be necessary to see your input files. 

If you would like to contribute any changes to the code, just follow these steps!
1. Fork the repository 
2. Make changes
3. Submit a pull request

That's it! If you have an idea for a feature that would be helpful for your research, but aren't sure how to implement it, this can be logged using the GitHub Issues page linked above and I'll do my best to help! 
