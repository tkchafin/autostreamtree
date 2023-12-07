import sys
from unittest.mock import patch
from autostreamtree.params import parseArgs


def test_parse_args_basic():
    test_args = ["autostreamtree", "-o", "my_prefix", "-v",
                 "vcf_file", "-s", "shapefile.shp", "-i", "coords_file",
                 "-p", "popmap_file"]
    sys.argv = test_args
    args = parseArgs()
    assert args.out == "my_prefix"
    assert args.geodb == "coords_file"
    assert args.vcf == "vcf_file"
    assert args.shapefile == "shapefile.shp"
    assert args.pop == "popmap_file"


def test_no_shapefile():
    test_args = ["autostreamtree", "-o", "my_prefix", "-v",
                 "vcf_file", "-i", "coords_file", "-r", "ALL",
                 "-p", "popmap_file"]
    sys.argv = test_args
    with patch('sys.exit') as mock_exit, patch('builtins.print') as mock_print:
        _ = parseArgs()
        mock_exit.assert_called()
        mock_print.assert_called()


def test_no_shapefile_ok():
    test_args = ["autostreamtree", "-o", "my_prefix", "-v",
                 "vcf_file", "-i", "coords_file", "-r", "GENDIST",
                 "-p", "popmap_file"]
    sys.argv = test_args
    args = parseArgs()
    assert args.shapefile is None
    assert args.run == "GENDIST"


def test_invalid_option():
    test_args = ["autostreamtree", "-invalid_option="]
    sys.argv = test_args

    with patch('builtins.print') as mock_print, patch('sys.exit') as mock_exit:
        _ = parseArgs()
        mock_exit.assert_called()
        mock_print.assert_called()


def test_invalid_dist():
    test_args = ["autostreamtree", "-o", "my_prefix", "-v",
                 "vcf_file", "-i", "coords_file", "-r", "GENDIST",
                 "-d", "FST"]
    sys.argv = test_args
    with patch('sys.exit') as mock_exit, patch('builtins.print') as mock_print:
        _ = parseArgs()
        mock_exit.assert_called()
        mock_print.assert_called()


def test_display_help_called():
    test_args = ["autostreamtree", "-h"]
    sys.argv = test_args
    with patch('sys.exit') as mock_exit, patch('builtins.print') as mock_print:
        _ = parseArgs()
        mock_exit.assert_called()
        mock_print.assert_called()


# default values
def test_parse_args_defaults():
    test_args = ["autostreamtree", "-v", "vcf_file", "-s", "shapefile.shp",
                 "-i", "coords_file", "-p", "popmap_file"]
    sys.argv = test_args
    args = parseArgs()

    # Assertions to check if defaults are set correctly
    assert args.out == "out"
    assert args.concat == "none"
    assert args.run == "ALL"
    assert args.network is None
    assert args.overwrite is False
    assert args.geopop is False
    assert args.clusterpop is False
    assert args.sclusterpop is False
    assert args.output_driver == "GPKG"
    assert args.dist == "LINFST"
    assert args.het is False
    assert args.genmat is None
    assert args.log is False
    assert args.and_log is False
    assert args.iterative is True
    assert args.weight == "CSE67"
    assert args.permutations == 1000
    assert args.coercemat is False
    assert args.method == "PEARSON"
    assert args.edge_list is None
    assert args.plots is False
    assert args.median is False
    assert args.ploidy == 2
    assert args.global_het is False
    assert args.loc_agg == "ARITH"
    assert args.pop_agg == "ARITH"
    assert args.sdist_agg == "ARITH"


# non-default values
def test_parse_args_nondefaults_v1():
    test_args = [
        "autostreamtree", "-o", "custom_prefix", "-v", "vcf_file", "-s",
        "shapefile.shp", "-i", "coords_file", "-p", "popmap_file",
        "--epsilon", "30", "-C", "LOC", "-n", "network_file", "--overwrite",
        "-r", "IBD", "--geopop", "--reachid_col", "R", "--length_col", "L",
        "--genmat", "genmat_file", "-d", "JOST", "--coercemat", "--het",
        "--min_samples", "2", "-L", "MEDIAN", "-P", "HARM", "--perm", "100",
        "--and_log", "-w", "FM67"
    ]
    sys.argv = test_args
    args = parseArgs()

    # Assertions to check if non-defaults are set correctly
    assert args.out == "custom_prefix"
    assert args.epsilon == 30
    assert args.concat == "loc"
    assert args.network == "network_file"
    assert args.overwrite is True
    assert args.run == "IBD"
    assert args.geopop is True
    assert args.reachid_col == "R"
    assert args.length_col == "L"
    assert args.genmat == "genmat_file"
    assert args.dist == "JOST"
    assert args.coercemat is True
    assert args.het is True
    assert args.min_samples == 2
    assert args.loc_agg == "MEDIAN"
    assert args.pop_agg == "HARM"
    assert args.permutations == 100
    assert args.and_log is True
    assert args.weight == "FM67"


# non-default values 2
def test_parse_args_nondefaults_v2():
    test_args = [
        "autostreamtree", "-o", "custom_prefix", "-v", "vcf_file", "-s",
        "shapefile.shp", "-i", "coords_file", "-p", "popmap_file",
        "--epsilon", "30", "-n", "network_file", "--overwrite",
        "-r", "RUNLOCI", "--geopop", "--reachid_col", "R", "--length_col", "L",
        "--genmat", "genmat_file", "-d", "PDIST", "--global_het",
        "--min_samples", "2", "-L", "MEDIAN", "-P", "ADJHARM",
        "--and_log", "-w", "CSE"
    ]
    sys.argv = test_args
    args = parseArgs()

    # Assertions to check if non-defaults are set correctly
    assert args.out == "custom_prefix"
    assert args.epsilon == 30
    assert args.concat == "none"
    assert args.network == "network_file"
    assert args.overwrite is True
    assert args.run == "RUNLOCI"
    assert args.geopop is True
    assert args.reachid_col == "R"
    assert args.length_col == "L"
    assert args.genmat == "genmat_file"
    assert args.dist == "PDIST"
    assert args.coercemat is False
    assert args.global_het is True
    assert args.min_samples == 2
    assert args.loc_agg == "MEDIAN"
    assert args.pop_agg == "ADJHARM"
    assert args.permutations == 1000
    assert args.and_log is True
    assert args.weight == "CSE67"
