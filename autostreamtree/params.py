import sys
import getopt


# Object to parse command-line arguments
class parseArgs():
    def __init__(self):
        # Define options
        try:
            options, _ = getopt.getopt(
                sys.argv[1:], 'hs:i:r:p:d:a:lw:o:gP:L:Scn:G:v:C:O:',
                ["shp=", "help", "input=", "run=", "pop=", "popmap=", "dist=",
                 "agg_method=", "het", "genmat=", "snp", "snps", "msat",
                 "msats", "log", "and_log", "iterative", "weight=", "out=",
                 "method=", "plots", "plot", "perm=", "phased", "median",
                 "diploid", "geopop", "geopops", "global_het", "haploid",
                 "loc_agg=", "pop_agg=", "sdist_agg=", "clusterpop",
                 "epsilon=", "min_samples=", "sclusterpop", "network=",
                 "overwrite", "reachid_col=", "length_col=", "coercemat",
                 "locmatdir=", "vcf=", "concat=", "edge_list=", "gdf_out="]
            )
        except getopt.GetoptError as err:
            print(err)
            self.display_help(
                "\nExiting because getopt returned non-zero exit status.")
        # Default values for params
        # Input params
        self.shapefile = None
        self.locmatdir = None
        self.geodb = None
        self.vcf = None
        self.concat = "none"
        self.run = "ALL"
        self.network = None
        self.pop = False
        self.overwrite = False
        self.geopop = False
        self.clusterpop = False
        self.sclusterpop = False
        self.output_driver = "GPKG"
        self.dist = "LINFST"
        self.het = False
        self.genmat = None
        self.log = False
        self.and_log = False
        self.iterative = True
        self.weight = "CSE67"
        self.permutations = 1000
        self.coercemat = False
        self.method = "PEARSON"
        self.edge_list = None
        self.plots = False
        self.out = "out"
        self.median = False
        self.ploidy = 2
        self.global_het = False
        self.loc_agg = "ARITH"
        self.pop_agg = "ARITH"
        self.sdist_agg = "ARITH"

        self.reachid_col = "HYRIV_ID"
        self.length_col = "LENGTH_KM"

        # dbscan Options
        self.min_samples = 1
        self.epsilon = 20

        # First pass to see if help menu was called
        for o, a in options:
            if o in ("-h", "-help", "--help"):
                self.display_help("Exiting because help menu was called.")

        # Second pass to set all args.
        for opt, arg in options:
            arg = arg.strip()
            opt = opt.replace("-", "")

            if opt in ('s', 'shp'):
                self.shapefile = arg
            elif opt in ('h', 'help'):
                pass
            elif opt in ('i', 'coords'):
                self.geodb = arg
            elif opt in ("C", "concat"):
                arg_lower = str(arg).lower()
                if arg_lower in ("none", "all", "loc"):
                    self.concat = arg_lower
                else:
                    self.display_help(f"Invalid option {arg_lower} \
                                      for option <--concat>")
            elif opt in ('v', 'vcf'):
                self.vcf = arg
            elif opt in ('r', 'run'):
                self.run = arg.upper()
                valid_runs = ["ALL", "GENDIST", "IBD", "STREAMDIST",
                              "STREAMTREE", "DISTANCES", "RUNLOCI"]
                if self.run not in valid_runs:
                    self.display_help(f"Invalid option {self.run} \
                                      for option <-r/--run>")
            elif opt in ('p', 'pop', "popmap"):
                self.pop = arg
            elif opt in ("g", "geopop", "geopops"):
                self.geopop = True
            elif opt in ('d', 'dist'):
                self.dist = arg.upper()
                valid_dists = ["PDIST", "FST",
                               "LINFST", "JOST", "NEI83",
                               "CHORD", "GST", "HARMD"]
                if self.dist not in valid_dists:
                    self.display_help(f"Invalid option {self.dist} for \
                                      option <-d/--dist>")
            elif opt == "het":
                self.het = True
            elif opt == "locmatdir":
                self.locmatdir = arg
            elif opt == "reachid_col":
                self.reachid_col = arg
            elif opt == "length_col":
                self.length_col = arg
            elif opt in ("genmat", "G"):
                self.genmat = arg
            elif opt in ('l', 'log'):
                self.log = True
            elif opt == "and_log":
                self.and_log = True
            elif opt == "edge_list":
                self.edge_list = arg
            elif opt in ("clusterpop", "c"):
                self.clusterpop = True
            elif opt == "epsilon":
                self.epsilon = float(arg)
            elif opt == "min_samples":
                self.min_samples = int(arg)
            elif opt in ('w', 'weight'):
                self.weight = arg.upper()
                weight_mappings = {
                    "FM": "FM67", "1/D": "FM67",
                    "BEYER": "BEYER74", "1/D^2": "BEYER74",
                    "1": "CSE67", "CSE": "CSE67"
                }
                self.weight = weight_mappings.get(self.weight, self.weight)
                valid_weights = ["CSE67", "FM67", "BEYER74"]
                if self.weight not in valid_weights:
                    self.display_help(f"Invalid option {self.weight} \
                                      for option <-w/--weight>")
            elif opt in ('o', 'out'):
                self.out = arg
            elif opt == "stream_fit":
                self.stream_fit = True
            elif opt == "overwrite":
                self.overwrite = True
            elif opt == "perm":
                self.permutations = int(arg)
            elif opt == "method":
                print("Sorry: Option --method is not yet implemented.")
                sys.exit(0)
            elif opt in ("plot", "plots"):
                self.plots = True
            elif opt == "phased":
                self.phased = True
                print("WARNING: Option <--snp> not yet implemented")
            elif opt == "median":
                self.median = True
            elif opt == "diploid":
                self.ploidy = 2
            elif opt == "haploid":
                self.ploidy = 1
            elif opt in ("O", "gdf_out"):
                arg_upper = str(arg).upper()
                if arg_upper not in ["GPKG", "SHP", "GDB"]:
                    self.display_help(f"Invalid option {arg_upper} for \
                                      option <--gdf_out>")
                self.output_driver = arg_upper
            elif opt == "global_het":
                self.global_het = True
            elif opt == "coercemat":
                self.coercemat = True
            elif opt in ("network", "n"):
                self.network = arg
            elif opt in ("pop_agg", "P"):
                self.pop_agg = arg.upper()
                valid_agg_methods = ["HARM", "ADJHARM", "ARITH",
                                     "GEOM", "MEDIAN", "MAX", "MIN"]
                if self.pop_agg not in valid_agg_methods:
                    self.display_help(f"Invalid option {self.pop_agg} for \
                                      option <--pop_agg>")
            elif opt in ("loc_agg", "L"):
                self.loc_agg = arg.upper()
                valid_agg_methods = ["HARM", "ADJHARM", "ARITH", "GEOM",
                                     "MEDIAN", "MAX", "MIN"]
                if self.loc_agg not in valid_agg_methods:
                    self.display_help(f"Invalid option {self.loc_agg} for \
                                      option <--loc_agg>")
            else:
                assert False, f"Unhandled option {opt}"

        if not self.geodb:
            self.display_help("No input provided <-i/--coords>")
        if not self.shapefile and self.run != "GENDIST":
            self.display_help("No shapefile provided <-s/--shp>")
        if self.ploidy > 2 or self.ploidy < 1:
            self.display_help(
                f"Ploidy of {self.ploidy} not currently \
                    allowable. Please choose 1 (haploid) or 2 (diploid)")

        # sanity checks
        if self.dist not in ["PDIST"]:
            if not self.pop and not self.geopop:
                self.display_help(
                    f"ERROR: Distance metric {self.dist} not possible \
                        without --pop or --geopop data."
                )

    def display_help(self, message=None):
        if message is not None:
            print("\n", message)

        print("\nautostreamtree\n")
        print("Author: Tyler K Chafin, Biomathematics and Statistics Scotland")
        print("Description: Methods for analysing genetic distances in \
networks.")
        print(
            "\nMandatory arguments:\n"
            "    -s, --shp       : Path to shapefile containing cleaned, \
contiguous stream reaches\n"
            "                       (can also support geodatabase or GPKG \
files)\n"
            "    -i, --input     : Input .tsv file containing sample \
coordinates\n"
            "    -v, --vcf       : Input VCF file containing genotypes\n\n"

            "General options:\n"
            "    -o, --out       : Output prefix [default=\"out\"]\n"
            "    -O, --gdf_out   : Output driver for annotated geodataframe \
(options \"SHP\", \"GPKG\", \"GDB\")\n"
            "    -C, --concat    : Concatenate all SNPs (\"all\"), by locus \
(\"loc\"), or not at all (\"none\")\n"
            "    -n, --network   : Provide an already optimized network \
output from a previous run\n"
            "    --overwrite     : Overwrite an input network (Only relevant \
with --network)\n"
            "    -h, --help      : Displays help menu\n"
            "    -r, --run       : Run which steps? Options:\n"
            "        ALL         : Run all steps\n"
            "        GENDIST     : Only calculate genetic distance matrix\n"
            "        STREAMDIST  : Only compute pairwise stream distances\n"
            "        DISTANCES   : Only compute GENDIST + STREAMDIST\n"
            "        IBD         : GENDIST + STREAMDIST + Mantel test\n"
            "        STREAMTREE  : GENDIST + STREAMDIST + fit StreamTree \
model\n"
            "        RUNLOCI     : Run STREAMTREE fitting on each locus\n"
            "    -p, --pop       : Pool individuals based on an input \
population map tsv file\n"
            "        NOTE: The location will be taken as the centroid among \
individual samples\n"
            "    -g, --geopop    : Pool individuals having identical \
coordinates\n"
            "    -c, --clusterpop: Use DBSCAN algorithm to automatically \
cluster populations\n"
            "    --reachid_col   : Attribute name representing primary key in \
shapefile [default=\"HYRIV_ID\"]\n"
            "    --length_col    : Attribute name giving length in kilometers \
[default=\"LENGTH_KM\"]\n\n"

            "Genetic distance options:\n"
            "    -d, --dist      : Use which metric of distance? Options:\n"
            "        Individual-based:\n"
            "          PDIST     : Uncorrected p-distances [# Differences / \
Length]\n"
            "        Frequency models (when using --pop):\n"
            "          FST       : Weir and Cockerham's Fst formulation \
(=THETAst)\n"
            "          LINFST    : [default] Rousset's (1997) Fst \
[=Fst/(1-Fst)]\n"
            "          JOST      : Global estimate Jost's (2008) D\n"
            "          HARMD     : Harmonic mean of Jost's D per-locus\n"
            "          NEI83     : Nei and Chesser (1983) Da\n"
            "          CHORD     : Cavalli-Sforza and Edwards (1967) chord \
distance\n"
            "          --NOTE: Individual-based metrics can also be computed \
for\n"
            "                    populations. You can set how these are \
aggregated w/ --pop_agg\n"
            "          --NOTE: Multiple loci for PDIST\n"
            "                  will be reported using the method defined in \
--loc_agg\n"
            "    -G, --genmat    : Skip calculation and use the provided \
labeled .tsv matrix\n"
            "    --coercemat     : [Boolean] Coerce negative values in input \
matrix to zero\n"
            "    --het           : [Boolean] Count partial differences [e.g. \
ind1=T, ind2=W]\n"
            "    --global_het    : Estimate Ht using global frequencies \
(default is averaged over pops)\n\n"

            "DBSCAN options (only when --clusterpop):\n"
            "    --min_samples   : Minimum samples per cluster [default=1]\n"
            "    --epsilon       : Maximum distance (in km) within a cluster \
[default=20]\n\n"

            "Aggregation options:\n"
            "    -P, --pop_agg   : Define aggregator function for certain \
genetic distances in pop samples\n"
            "    -L, --loc_agg   : Define aggregator function for aggregating \
locus-wise distances\n"
            "        All of these can take the following options:\n"
            "          ARITH     : [default] Use arithmetic mean\n"
            "          MEDIAN    : Use median distance\n"
            "          HARM      : Use harmonic mean\n"
            "          ADJHARM   : Adjusted harmonic mean (see docs)\n"
            "          GEOM      : Use geometric mean\n"
            "          MIN       : Use minimum distance\n"
            "          MAX       : Use maximum distance\n\n"

            "IBD options:\n"
            "    --perm          : Number of permutations for mantel test \
[def=1000]\n"
            "    --and_log       : Also perform IBD steps with log geographic \
distances\n\n"

            "StreamTree options (see Kalinowski et al. 2008) :\n"
            "    -w, --weight    : Desired weighting for least-squares \
fitting:\n"
            "        Options:\n"
            "          FM67      : Fitch and Margoliash (1967) [w = 1/D^2]\n"
            "          BEYER74   : Beyer et al. (1974) weights [w = 1/D]\n"
            "          CSE67     : [default] Cavalli-Sforza & Edwards (1967) \
[w = 1]\n\n"
        )
        print()
        sys.exit()
