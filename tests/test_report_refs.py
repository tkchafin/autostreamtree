import sys
from autostreamtree.params import parseArgs
from autostreamtree.report_refs import fetch_references


def test_fetch_references_default():
    test_args = ["autostreamtree", "-v", "vcf_file", "-s", "shapefile.shp",
                 "-i", "coords_file", "-p", "popmap_file"]
    sys.argv = test_args
    args = parseArgs()

    refs = fetch_references(args)
    assert refs


def test_fetch_references_nondefault():
    test_args = ["autostreamtree", "-v", "vcf_file", "-s", "shapefile.shp",
                 "-i", "coords_file", "-p", "popmap_file", "-w", "FM67",
                 "-d", "JC69", "--clusterpop", "--loc_agg", "ADJHARM"]
    sys.argv = test_args
    args = parseArgs()

    refs = fetch_references(args)
    assert refs


def test_fetch_references_nondefault_2():
    test_args = ["autostreamtree", "-v", "vcf_file", "-s", "shapefile.shp",
                 "-i", "coords_file", "-p", "popmap_file", "-w", "BEYER74",
                 "-d", "JOST", "-r", "IBD"]
    sys.argv = test_args
    args = parseArgs()

    refs = fetch_references(args)
    assert refs


def test_fetch_references_nondefault_3():
    test_args = ["autostreamtree", "-v", "vcf_file", "-s", "shapefile.shp",
                 "-i", "coords_file", "-p", "popmap_file", "-w", "CSE67",
                 "-d", "CHORD", "-r", "STREAMTREE"]
    sys.argv = test_args
    args = parseArgs()

    refs = fetch_references(args)
    assert refs


def test_fetch_references_nondefault_4():
    test_args = ["autostreamtree", "-v", "vcf_file", "-s", "shapefile.shp",
                 "-i", "coords_file", "-p", "popmap_file", "-w", "CSE67",
                 "-d", "NEI83", "--pop_agg", "ADJHARM", "-r", "ALL", 
                 "--loc_agg", "ARITH", "-g", "genmat_file"]
    sys.argv = test_args
    args = parseArgs()

    refs = fetch_references(args)
    assert refs