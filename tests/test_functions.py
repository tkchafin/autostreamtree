import os
import pytest
from io import StringIO
import sys
import tempfile
import numpy as np
import pandas as pd
import pickle
import networkx as nx
from scipy import stats
from sortedcontainers import SortedDict
from unittest.mock import patch, mock_open, MagicMock
import geopandas as gpd
import autostreamtree
import autostreamtree.functions as funcs

pytestmark = pytest.mark.filterwarnings("ignore:Geometry is in a geographic CRS. Results from 'length' are likely incorrect. Use 'GeoSeries.to_crs()' to re-project geometries to a projected CRS before this operation.:UserWarning")


@pytest.fixture
def sample_distances():
    return np.array([[0, 1, 1],
                     [1, 0, 2],
                     [3, 1, 1]])


@pytest.fixture
def sample_incidence_matrix():
    return np.array([[0, 0, 1],
                     [0, 1, 0],
                     [1, 0, 0]])


@pytest.fixture
def empty_graph():
    return nx.Graph()


@pytest.fixture
def sample_graph():
    G = nx.Graph()
    G.add_edge(1, 2, weight=10, reachid_col='a')
    G.add_edge(2, 3, weight=20, reachid_col='b')
    G.add_edge(3, 4, weight=30, reachid_col='c')
    return G


@pytest.fixture
def sample_gen_matrix():
    return np.array([[0, 2, 3], [2, 0, 1], [3, 1, 0]])


@pytest.fixture
def sample_pop_gen_matrix():
    return np.array([[0, 1], [1, 0]])

@pytest.fixture
def sample_geo_matrix():
    return np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])


@pytest.fixture
def mock_network_path():
    G = nx.Graph()
    G.add_edge('A', 'B')
    G.add_edge('B', 'C')
    with tempfile.NamedTemporaryFile(mode='wb', suffix='.pickle',
                                     delete=False) as temp_file:
        pickle.dump(G, temp_file)
        return temp_file.name


@pytest.fixture
def point_coords():
    return {'A': (0, 0), 'B': (1, 1), 'C': (2, 2)}


@pytest.fixture
def pop_coords():
    return {'Pop1': (0, 0), 'Pop2': (1, 1)}


@pytest.fixture
def shapefile_path():
    base_path = os.path.dirname(autostreamtree.__file__)
    file_path = os.path.join(base_path, 'data', 'test.shp')
    return file_path


@pytest.fixture
def vcf_file_path():
    with tempfile.NamedTemporaryFile(
            mode='w', suffix='.vcf', delete=False) as temp_file:
        temp_file.write("""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tInd1\tInd2\tInd3\tInd4
1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t0/1\t1/1
1\t150\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\t0/1\t0/1
2\t200\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0\t1/1\t1/1
""")
        return temp_file.name


@pytest.fixture
def mock_points_df():
    data = {'sample': ['Ind1', 'Ind2', 'Ind3', 'Ind4'],
            'lat': [26.925414, 26.925414, 27.21041, 27.21041],
            'long': [90.396139, 90.396139, 91.20795, 91.20795]}
    return pd.DataFrame(data)


@pytest.fixture
def mock_network_graph():
    G = nx.Graph()
    G.add_edge((26.92, 90.39), (27.21, 91.20))
    return G


@pytest.fixture
def mock_popmap():
    return SortedDict({
        'Pop1': ['Ind1', 'Ind2'],
        'Pop2': ['Ind3', 'Ind4']
    })


@pytest.fixture
def mock_coords():
    return SortedDict({
        'Ind1': (90.396139, 26.925414),
        'Ind2': (90.396139, 26.925414),
        'Ind3': (91.20795, 27.21041),
        'Ind4': (91.20795, 27.21041)
    })


@pytest.fixture
def mock_params():
    class MockParams:
        def __init__(self):
            self.coercemat = False
            self.pop = False
            self.geopop = False
            self.clusterpop = False
            self.pop_agg = "ARITH"
            self.loc_agg = "ARITH"
            self.ploidy = 2
            self.global_het = False
            self.run = "ALL"
            self.out = os.path.join(tempfile.gettempdir(), "test_output")

    return MockParams()


def test_read_vcf(vcf_file_path):
    result = funcs.read_vcf(vcf_file_path, concat="none")
    expected = {
        'Ind1': ['A/A', 'A/T', 'A/A'],
        'Ind2': ['A/T', 'A/T', 'A/A'],
        'Ind3': ['A/T', 'A/T', 'T/T'],
        'Ind4': ['T/T', 'A/T', 'T/T']
    }
    assert result == expected


@pytest.fixture
def mock_popmap_file(mock_popmap):
    # Create a temporary popmap file and write mock_popmap data to it
    with tempfile.NamedTemporaryFile(mode='w', delete=False,
                                     suffix='.txt') as temp_file:
        for pop, samples in mock_popmap.items():
            for sample in samples:
                temp_file.write(f"{sample}\t{pop}\n")
        return temp_file.name


@pytest.fixture
def mock_indmat():
    data = np.array([[0, 1, 2, 3], [1, 0, 4, 5], [2, 4, 0, 6], [3, 5, 6, 0]])
    samples = ['Ind1', 'Ind2', 'Ind3', 'Ind4']
    return pd.DataFrame(data, index=samples, columns=samples)


@pytest.fixture
def mock_gdf():
    return gpd.GeoDataFrame({'geometry': []})


@pytest.fixture
def output_prefix():
    return "test_output"


def mock_to_csv(*args, **kwargs):
    print("Mocking to_csv call")


def test_read_network(mock_network_path, shapefile_path):
    # Test reading from a saved network file
    G_pickle = funcs.read_network(network=mock_network_path, shapefile=None)
    assert isinstance(G_pickle, nx.Graph)
    assert len(G_pickle.nodes()) > 0  # Ensure it has nodes
    assert len(G_pickle.edges()) > 0  # Ensure it has edges

    # Test building a network from a shapefile
    G_shapefile = funcs.read_network(network=None, shapefile=shapefile_path)
    assert isinstance(G_shapefile, nx.Graph)
    assert len(G_shapefile.nodes()) > 0  # Ensure it has nodes
    assert len(G_shapefile.edges()) > 0  # Ensure it has edges


def test_get_loc_data(vcf_file_path):
    seqs = funcs.read_vcf(vcf_file_path, concat="none")
    loc_data_generator = funcs.get_loc_data(seqs)
    loc_data_list = list(loc_data_generator)

    # Expected output for each locus
    expected_loci_data = [
        {'Ind1': ['A/A'], 'Ind2': ['A/T'], 'Ind3': ['A/T'], 'Ind4': ['T/T']},
        {'Ind1': ['A/T'], 'Ind2': ['A/T'], 'Ind3': ['A/T'], 'Ind4': ['A/T']},
        {'Ind1': ['A/A'], 'Ind2': ['A/A'], 'Ind3': ['T/T'], 'Ind4': ['T/T']}
    ]

    # Check if the generated data matches the expected output for each locus
    assert len(loc_data_list) == len(expected_loci_data)
    for generated, expected in zip(loc_data_list, expected_loci_data):
        assert generated == expected


def test_block_and_enable_print():
    original_stdout = sys.stdout
    captured_output = StringIO()
    sys.stdout = captured_output

    # block_print and try to print something
    funcs.block_print()
    print("This should not be printed")
    assert captured_output.getvalue() == ""
    sys.stdout = original_stdout

    # enable print and test again
    funcs.enable_print()
    captured_output = StringIO()
    sys.stdout = captured_output

    print("This should be printed")
    assert captured_output.getvalue() == "This should be printed\n"
    sys.stdout = original_stdout


def test_process_samples(mock_params, mock_points_df, mock_network_graph,
                         mock_popmap_file):
    mock_params.pop = mock_popmap_file

    # Call the function with the mocked data
    point_coords, pop_coords, popmap = funcs.process_samples(
        mock_params, mock_points_df, mock_network_graph)

    # Assertions to verify the outputs are as expected
    # Verify point_coords
    assert len(point_coords) == len(mock_points_df)
    for name, coords in point_coords.items():
        assert name in mock_points_df['sample'].values
        assert isinstance(coords, tuple)

    # Verify pop_coords and popmap if applicable
    if mock_params.pop:
        assert len(pop_coords) == len(popmap)
        for pop, members in popmap.items():
            assert pop in pop_coords
            for member in members:
                assert member in mock_points_df['sample'].values


def test_process_samples_v2(mock_params, mock_points_df, mock_network_graph):
    mock_params.pop = False
    mock_params.geopop = True

    # Call the function with the mocked data
    point_coords, pop_coords, popmap = funcs.process_samples(
        mock_params, mock_points_df, mock_network_graph)

    # Assertions to verify the outputs are as expected
    # Verify point_coords
    assert len(point_coords) == len(mock_points_df)
    for name, coords in point_coords.items():
        assert name in mock_points_df['sample'].values
        assert isinstance(coords, tuple)

    # Verify pop_coords and popmap if applicable
    if mock_params.pop:
        assert len(pop_coords) == len(popmap)
        for pop, members in popmap.items():
            assert pop in pop_coords
            for member in members:
                assert member in mock_points_df['sample'].values


def test_process_samples_clusterpop(mock_params, mock_points_df,
                                    mock_network_graph):
    mock_params.pop = False
    mock_params.clusterpop = True
    mock_params.min_samples = 1
    mock_params.epsilon = 10

    # Call the function with the mocked data
    point_coords, pop_coords, popmap = funcs.process_samples(
        mock_params, mock_points_df, mock_network_graph)

    # Assertions to verify the outputs are as expected
    # Clusterpop should result in some clustering
    assert len(pop_coords) <= len(point_coords)
    assert len(popmap) <= len(point_coords)
    for name, coords in point_coords.items():
        assert name in mock_points_df['sample'].values
        assert isinstance(coords, tuple)

    # Verify that each cluster has at least min_samples
    for _, members in popmap.items():
        assert len(members) >= mock_params.min_samples


def test_get_lower_tri():
    mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    expected = np.array([4, 7, 8])
    assert np.array_equal(funcs.get_lower_tri(mat), expected)


def test_great_circle():
    # Define two known points (Latitude, Longitude)
    coord1 = (48.8566, 2.3522)  # Paris
    coord2 = (51.5074, -0.1278)  # London

    # Known distance between Paris and London in kilometers
    known_distance = 344
    calculated_distance = funcs.great_circle(coord1[1], coord1[0],
                                             coord2[1], coord2[0])
    assert abs(calculated_distance - known_distance) < 1  # 1 km tolerance


def test_nCr():
    assert funcs.nCr(5, 2) == 10
    assert funcs.nCr(10, 5) == 252
    assert funcs.nCr(0, 0) == 1


def test_get_point_table(mock_coords):
    result_df = funcs.get_point_table(mock_coords)

    assert isinstance(result_df, pd.DataFrame)
    assert list(result_df.columns) == ['sample', 'lat', 'long']
    assert len(result_df) == len(mock_coords)

    for index, row in result_df.iterrows():
        sample = row['sample']
        long, lat = mock_coords[sample]
        assert row['lat'] == lat
        assert row['long'] == long


def test_r2():
    # Create mock data
    x = np.array([1, 2, 3, 4, 5])
    y = np.array([2, 4, 6, 8, 10])

    expected_r2 = 1 ** 2
    result = funcs.r2(x, y)
    assert result == expected_r2

    x2 = np.array([1, 2, 3, 4, 5])
    y2 = np.array([2, 3, 2, 5, 7])
    expected_r2_less_correlated = stats.pearsonr(x2, y2)[0] ** 2
    result2 = funcs.r2(x2, y2)

    assert result2 == pytest.approx(expected_r2_less_correlated)


@pytest.mark.parametrize("output_driver, expected_extension, expected_driver",
[
    ("SHP", ".shp", "ESRI SHAPEFILE"),
    ("GPKG", ".gpkg", "GPKG"),
    ("GDB", ".gdb", "GDB")
])
def test_write_geodataframe_extension(mock_gdf, output_prefix, output_driver,
                                      expected_extension, expected_driver):
    with patch("geopandas.GeoDataFrame.to_file") as mock_to_file:
        funcs.write_geodataframe(mock_gdf, output_prefix, output_driver)
        output_path = f"{output_prefix}{expected_extension}"
        mock_to_file.assert_called_with(output_path,
                                        driver=expected_driver.upper())


def test_write_geodataframe_gdb_directory_creation(mock_gdf, output_prefix):
    output_driver = "GDB"
    output_path = f"{output_prefix}.gdb"

    with patch("os.path.exists", return_value=False) as mock_exists:
        with patch("os.makedirs") as mock_makedirs:
            with patch("geopandas.GeoDataFrame.to_file") as mock_to_file:
                funcs.write_geodataframe(mock_gdf, output_prefix,
                                         output_driver)
                mock_exists.assert_called_once_with(output_path)
                mock_makedirs.assert_called_once_with(output_path)
                mock_to_file.assert_called_with(output_path,
                                                driver=output_driver)


def test_write_geodataframe_file_writing(mock_gdf, output_prefix):
    output_driver = "GPKG"
    output_path = f"{output_prefix}.gpkg"

    with patch("geopandas.GeoDataFrame.to_file") as mock_to_file:
        funcs.write_geodataframe(mock_gdf, output_prefix, output_driver)
        mock_to_file.assert_called_with(output_path,
                                        driver=output_driver)


def test_find_pair():
    lst = [1, 2, 3, 4, 5]
    assert funcs.find_pair(lst, 1, 2) is True
    assert funcs.find_pair(lst, 3, 4) is True
    assert funcs.find_pair(lst, 4, 3) is True
    assert funcs.find_pair(lst, 1, 3) is False
    assert funcs.find_pair(lst, 6, 7) is False


@pytest.mark.parametrize("d, weight, expected_diagonal", [
    (np.array([1, 2, 3]), "CSE67", np.array([1, 1, 1])),
    (np.array([1, 2, 0]), "BEYER74", np.array([1, 0.5, 0])),
    (np.array([1, 2, 0]), "FM67", np.array([1, 0.25, 0])),
    (np.array([1, 2, 3]), "UNKNOWN", np.array([1, 1, 1]))
])
def test_generate_weights_matrix(d, weight, expected_diagonal):
    W = funcs.generate_weights_matrix(d, weight)
    np.testing.assert_array_equal(np.diag(W), expected_diagonal)


# Tests for vectorize_mat
@pytest.mark.parametrize("mat, expected_length", [
    (np.array([[0, 1], [1, 0]]), 1),
    (np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]]), 3),
])
def test_vectorize_mat(mat, expected_length):
    vec = funcs.vectorize_mat(mat)
    assert len(vec) == expected_length


@pytest.mark.parametrize("iterative, expected_output", [
    (False, np.array([2., 1., 1.])),
    (True, np.array([2., 1., 1.]))
])
def test_fit_least_squares_distances(sample_distances,
                                     sample_incidence_matrix,
                                     output_prefix, iterative,
                                     expected_output):
    with patch("pandas.DataFrame.to_csv", new=mock_to_csv):
        result = funcs.fit_least_squares_distances(sample_distances,
                                                   sample_incidence_matrix,
                                                   iterative, output_prefix)
        np.testing.assert_array_almost_equal(result, expected_output)


def test_test_ibd(sample_gen_matrix, sample_geo_matrix):
    with patch('mantel.test') as mock_mantel, patch('autostreamtree.functions.pd.DataFrame.to_csv') as mock_to_csv:
        mock_mantel.return_value = MagicMock(r=0.5, p=0.01, z=2.5)
        try:
            funcs.test_ibd(sample_gen_matrix, sample_geo_matrix,
                           'output', 999, log=False)
            funcs.test_ibd(sample_gen_matrix, sample_geo_matrix,
                           'output', 999, log=True)
            no_errors = True
        except Exception:
            no_errors = False

        assert no_errors, "test_ibd function raised an error"


def test_output_fitted_d():
    pred = pd.DataFrame({'observed_D': [1, 2, 3], 'predicted_D': [1.1, 2.1, 3.1]})
    with patch('autostreamtree.functions.pd.DataFrame.to_csv') as mock_to_csv, patch('autostreamtree.functions.plt.savefig') as mock_savefig:
        funcs.output_fitted_d(pred, 'output')


def test_prune_graph(sample_graph):
    pruned = funcs.prune_graph(sample_graph, ['a', 'c'], 'reachid_col')
    assert len(pruned.edges) == 2
    assert (1, 2) in pruned.edges
    assert (3, 4) in pruned.edges
    assert len(list(nx.isolates(pruned))) == 0
    with pytest.raises(ValueError):
        _ = funcs.prune_graph(sample_graph, ['d'], 'reachid_col')
    assert len(sample_graph.edges) == 3


def test_report_genmats(sample_gen_matrix, sample_pop_gen_matrix,
                        point_coords, pop_coords, mock_params):
    with patch('builtins.print') as mock_print, patch('autostreamtree.functions.pd.DataFrame.to_csv') as mock_to_csv:
        funcs.report_genmats(mock_params, sample_gen_matrix, sample_pop_gen_matrix, point_coords, pop_coords)

        assert mock_print.call_count == 4 
        assert mock_to_csv.call_count == 2 
    with patch('builtins.print') as mock_print, patch('autostreamtree.functions.pd.DataFrame.to_csv') as mock_to_csv:
        funcs.report_genmats(mock_params, None, None, point_coords, pop_coords)
        assert mock_print.call_count == 0
        assert mock_to_csv.call_count == 0


def test_report_genmats_list(sample_gen_matrix, sample_pop_gen_matrix,
                             point_coords, pop_coords, mock_params):
    genlist = [sample_gen_matrix, None, sample_gen_matrix]
    popgenlist = [sample_pop_gen_matrix, sample_pop_gen_matrix]

    with patch('os.makedirs') as mock_makedirs, patch('autostreamtree.functions.pd.DataFrame.to_csv') as mock_to_csv:
        funcs.report_genmats_list(mock_params, genlist, popgenlist, point_coords, pop_coords)
        mock_makedirs.assert_called_once_with(mock_params.out + "_locmats", exist_ok=True)
        assert mock_to_csv.call_count == 4


def test_read_popmap():
    mock_file_content = "ind1 pop1\nind2 pop2\n\nind3 pop3"
    with patch(
            "builtins.open",
            mock_open(read_data=mock_file_content)) as mock_file:
        result = funcs.read_popmap("fake_popmap.txt")
        mock_file.assert_called_once_with("fake_popmap.txt", "r")
        expected_result = {"ind1": "pop1", "ind2": "pop2", "ind3": "pop3"}
        assert result == expected_result


def test_replace_zeroes():
    data_with_zeroes = np.array([0, 1, 2, 0, 3])
    expected_result = np.array([1, 1, 2, 1, 3])
    assert np.array_equal(funcs.replace_zeroes(data_with_zeroes),
                          expected_result)

    data_without_zeroes = np.array([1, 2, 3])
    assert np.array_equal(funcs.replace_zeroes(data_without_zeroes),
                          data_without_zeroes)

    data_all_zeroes = np.array([0, 0, 0])
    with pytest.raises(ValueError):
        _ = funcs.replace_zeroes(data_all_zeroes)
