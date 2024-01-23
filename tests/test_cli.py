import os
import sys
import pytest
import tempfile
import autostreamtree
import autostreamtree.cli as cli


pytestmark = pytest.mark.filterwarnings(
    "ignore:Geometry is in a geographic CRS. Results from 'length' are likely incorrect. Use 'GeoSeries.to_crs()' to re-project geometries to a projected CRS before this operation.:UserWarning",
    "ignore:Creating a 256th field, but some DBF readers might only support 255 fields:RuntimeWarning"
)

@pytest.fixture
def output_directory_path():
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir


@pytest.fixture
def shapefile_path():
    base_path = os.path.dirname(autostreamtree.__file__)
    file_path = os.path.join(base_path, 'data', 'test.shp')
    return file_path


@pytest.fixture
def coords_file_path():
    base_path = os.path.dirname(autostreamtree.__file__)
    file_path = os.path.join(base_path, 'data', 'test.coords')
    return file_path


@pytest.fixture
def vcf_file_path():
    base_path = os.path.dirname(autostreamtree.__file__)
    file_path = os.path.join(base_path, 'data', 'test_sub100.vcf.gz')
    return file_path


@pytest.fixture
def popmap_file_path():
    base_path = os.path.dirname(autostreamtree.__file__)
    file_path = os.path.join(base_path, 'data', 'test.popmap')
    return file_path


def test_autostreamtree_outputs(shapefile_path, coords_file_path,
                                vcf_file_path, popmap_file_path):
    expected_files = [
        "clusterCentroids.txt", "clusteredPoints.pdf", "full.network",
        "incidenceMatrix.txt", "isolationByDistance.pdf",
        "isolationByDistance.txt", "leastSquaresConstrained.txt", "network",
        "obsVersusFittedD.pdf", "obsVersusFittedD.txt", "pointCoords.txt",
        "popGenDistMat.txt", "reachToEdgeTable.txt", "snapDistances.pdf",
        "snapDistances.txt", "streamDistMat.txt", "streamTree.gpkg",
        "streamTree.txt", "streamsByEdgeID.pdf", "streamsByFittedD.pdf",
        "subGraph.pdf"
    ]

    with tempfile.TemporaryDirectory() as temp_dir:
        output_prefix = os.path.join(temp_dir, "test")
        args = [
            "autostreamtree", "-s", shapefile_path, "-i", coords_file_path,
            "-v", vcf_file_path, "-p", popmap_file_path, "-r", "ALL",
            "--reachid_col", "HYRIV_ID", "--length_col", "LENGTH_KM",
            "-o", output_prefix, "--seed", "1234", "--gdf_out", "GPKG"
        ]

        # Temporarily replace sys.argv
        original_argv = sys.argv
        sys.argv = args

        # Call the main function directly
        cli.main()

        # Restore original sys.argv
        sys.argv = original_argv

        # Check for the existence of each expected output file
        for file_suffix in expected_files:
            assert os.path.exists(f"{output_prefix}.{file_suffix}")


def test_autostreamtree_output_shp(shapefile_path, coords_file_path,
                                vcf_file_path, popmap_file_path):
    expected_files = [
        "clusterCentroids.txt", "clusteredPoints.pdf", "full.network",
        "incidenceMatrix.txt", "isolationByDistance.pdf",
        "isolationByDistance.txt", "leastSquaresConstrained.txt", "network",
        "obsVersusFittedD.pdf", "obsVersusFittedD.txt", "pointCoords.txt",
        "popGenDistMat.txt", "reachToEdgeTable.txt", "snapDistances.pdf",
        "snapDistances.txt", "streamDistMat.txt", "streamTree.shp",
        "streamTree.txt", "streamsByEdgeID.pdf", "streamsByFittedD.pdf",
        "subGraph.pdf"
    ]

    with tempfile.TemporaryDirectory() as temp_dir:
        output_prefix = os.path.join(temp_dir, "test")
        args = [
            "autostreamtree", "-s", shapefile_path, "-i", coords_file_path,
            "-v", vcf_file_path, "-p", popmap_file_path, "-r", "ALL",
            "--reachid_col", "HYRIV_ID", "--length_col", "LENGTH_KM",
            "-o", output_prefix, "--seed", "1234", "--gdf_out", "SHP"
        ]

        # Temporarily replace sys.argv
        original_argv = sys.argv
        sys.argv = args

        # Call the main function directly
        cli.main()

        # Restore original sys.argv
        sys.argv = original_argv

        # Check for the existence of each expected output file
        for file_suffix in expected_files:
            assert os.path.exists(f"{output_prefix}.{file_suffix}")


def test_autostreamtree_outputs_gendist(shapefile_path, coords_file_path,
                                        vcf_file_path, popmap_file_path):
    expected_files = [
        "popGenDistMat.txt"
    ]

    with tempfile.TemporaryDirectory() as temp_dir:
        output_prefix = os.path.join(temp_dir, "test")
        args = [
            "autostreamtree", "-s", shapefile_path, "-i", coords_file_path,
            "-v", vcf_file_path, "-p", popmap_file_path, "-r", "ALL",
            "--reachid_col", "HYRIV_ID", "--length_col", "LENGTH_KM",
            "-o", output_prefix, "--seed", "1234", "-r", "GENDIST"
        ]

        # Temporarily replace sys.argv
        original_argv = sys.argv
        sys.argv = args

        # Call the main function directly and catch SystemExit
        with pytest.raises(SystemExit) as exit_info:
            cli.main()
        assert exit_info.value.code == 0, "SystemExit with non-zero exit code"

        # Restore original sys.argv
        sys.argv = original_argv

        # Check for the existence of each expected output file
        for file_suffix in expected_files:
            assert os.path.exists(f"{output_prefix}.{file_suffix}")


def test_autostreamtree_outputs_sdist(shapefile_path, coords_file_path,
                                        vcf_file_path, popmap_file_path):
    expected_files = [
        "full.network", "network", "snapDistances.pdf",
        "snapDistances.txt", "streamDistMat.txt", "subGraph.pdf"
    ]

    with tempfile.TemporaryDirectory() as temp_dir:
        output_prefix = os.path.join(temp_dir, "test")
        args = [
            "autostreamtree", "-s", shapefile_path, "-i", coords_file_path,
            "-v", vcf_file_path, "-p", popmap_file_path, "-r", "ALL",
            "--reachid_col", "HYRIV_ID", "--length_col", "LENGTH_KM",
            "-o", output_prefix, "--seed", "1234", "-r", "STREAMDIST"
        ]

        # Temporarily replace sys.argv
        original_argv = sys.argv
        sys.argv = args

        # Call the main function directly and catch SystemExit
        with pytest.raises(SystemExit) as exit_info:
            cli.main()
        assert exit_info.value.code == 0, "SystemExit with non-zero exit code"

        # Restore original sys.argv
        sys.argv = original_argv

        # Check for the existence of each expected output file
        for file_suffix in expected_files:
            assert os.path.exists(f"{output_prefix}.{file_suffix}")