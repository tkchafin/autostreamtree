import pytest
from unittest.mock import mock_open, patch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import autostreamtree.cluster_pops as cluster_pops
from sortedcontainers import SortedDict


# Test for basic functionality
def test_dbscan_cluster_basic():
    coords = {'A': (33.4484, -112.0740), 'B': (37.7749, -122.4194),
              'C': (40.7128, -74.0060)}
    epsilon = 300
    min_samples = 2
    expected_result = SortedDict({'DB_-1': ['A', 'B', 'C']})
    assert cluster_pops.dbscan_cluster(coords, epsilon, min_samples) == \
        expected_result


# Test for invalid coordinates
@pytest.mark.parametrize("coords, epsilon, min_samples", [
    ({"A": (33.4484, None), "B": (37.7749, -122.4194)}, 300, 2),
    ({"A": "invalid", "B": (37.7749, -122.4194)}, 300, 2),
    # Add more cases as needed
])
def test_dbscan_cluster_invalid_coords(coords, epsilon, min_samples):
    with pytest.raises(ValueError):
        cluster_pops.dbscan_cluster(coords, epsilon, min_samples)


# Test for output file writing
@patch("builtins.open", new_callable=mock_open)
def test_dbscan_cluster_output_file(mock_file):
    coords = {'A': (33.4484, -112.0740), 'B': (37.7749, -122.4194)}
    epsilon = 300
    min_samples = 2
    output_file = "test_output.txt"
    cluster_pops.dbscan_cluster(coords, epsilon, min_samples, output_file)
    mock_file.assert_called_with(output_file, 'w')


def test_coords_to_dataframe():
    coords = {'A': (33.4484, -112.0740), 'B': (37.7749, -122.4194)}
    expected_df = pd.DataFrame(
        {'long': [33.4484, 37.7749], 'lat': [-112.0740, -122.4194]},
        index=['A', 'B']
    )
    pd.testing.assert_frame_equal(
        cluster_pops.coords_to_dataframe(coords), expected_df)


def test_coords_to_dataframe_empty():
    coords = {}
    # Create an empty DataFrame with an object index
    expected_df = pd.DataFrame(
        {'long': [], 'lat': []}, index=pd.Index([], dtype='object'))
    expected_df = expected_df.astype({'long': 'float64', 'lat': 'float64'})
    pd.testing.assert_frame_equal(
        cluster_pops.coords_to_dataframe(coords), expected_df)


def test_coords_to_matrix():
    coords = {'A': (33.4484, -112.0740), 'B': (37.7749, -122.4194)}
    expected_matrix = np.array([[33.4484, -112.0740], [37.7749, -122.4194]])
    assert np.array_equal(
        cluster_pops.coords_to_matrix(coords), expected_matrix)


def test_coords_to_matrix_empty():
    coords = {}
    expected_matrix = np.array([], dtype=float).reshape(0, 2)
    assert cluster_pops.coords_to_matrix(coords).shape == expected_matrix.shape


def test_get_cluster_centroid_basic():
    coords = SortedDict({'A': (33.4484, -112.0740), 'B': (37.7749, -122.4194)})
    popmap = SortedDict({'cluster1': ['A'], 'cluster2': ['B']})
    expected_centroids = SortedDict(
        {'cluster1': (33.4484, -112.0740), 'cluster2': (37.7749, -122.4194)})
    assert cluster_pops.get_cluster_centroid(coords, popmap) == \
        expected_centroids


def test_get_cluster_centroid_no_coords():
    coords = SortedDict({'A': (33.4484, -112.0740)})
    popmap = SortedDict({'cluster1': ['B']})  # No coordinate for 'B'
    with pytest.raises(ValueError):
        cluster_pops.get_cluster_centroid(coords, popmap)


@patch("builtins.open", new_callable=mock_open)
def test_get_cluster_centroid_output_file(mock_file):
    coords = SortedDict({'A': (33.4484, -112.0740)})
    popmap = SortedDict({'cluster1': ['A']})
    output_file = "test_output"
    cluster_pops.get_cluster_centroid(coords, popmap, output_file)
    mock_file.assert_called_with(output_file + ".clusterCentroids.txt", "w")


@pytest.mark.parametrize("coords, subset, expected, raises_exception, \
                         exception_type", [
    # Valid subsets
    (SortedDict({'A': (33.4484, -112.0740),
                 'B': (37.7749, -122.4194), 'C': (40.7128, -74.0060)}),
     ['A', 'C'],
     np.array([[33.4484, -112.0740], [40.7128, -74.0060]]),
     False, None),
    # Empty subset
    (SortedDict({'A': (33.4484, -112.0740), 'B': (37.7749, -122.4194)}),
     [],
     np.array([], dtype=float).reshape(0, 2),
     False, None),
    # Nonexistent subset elements
    (SortedDict({'A': (33.4484, -112.0740), 'B': (37.7749, -122.4194)}),
     ['X', 'Y'],
     None,
     True, ValueError),
    # Invalid dictionary input
    ("NotADict",
     ['A', 'C'],
     None,
     True, TypeError),
    # Invalid subset input
    (SortedDict({'A': (33.4484, -112.0740)}),
     "NotAList",
     None,
     True, TypeError),
    # Invalid coordinate format
    (SortedDict({'A': (33.4484,), 'B': (37.7749, -122.4194)}),
     ['A', 'B'],
     None,
     True, ValueError),
])
def test_get_pop_coords_matrix(coords, subset, expected, raises_exception,
                               exception_type):
    if raises_exception:
        with pytest.raises(exception_type):
            cluster_pops.get_pop_coords_matrix(coords, subset)
    else:
        result = cluster_pops.get_pop_coords_matrix(coords, subset)
        if expected is not None:
            assert np.array_equal(result, expected)
        else:
            # Check for empty array of shape (0, 2) when expected is None
            assert result.shape == (0, 2)

@pytest.mark.parametrize("point_coords, popmap, centroids", [
    # Test with valid inputs and no centroids
    (SortedDict({'A': (33.4484, -112.0740), 'B': (37.7749, -122.4194)}),
     SortedDict({'pop1': ['A'], 'pop2': ['B']}),
     None),

    # Test with valid inputs and centroids
    (SortedDict({'A': (33.4484, -112.0740), 'B': (37.7749, -122.4194)}),
     SortedDict({'pop1': ['A'], 'pop2': ['B']}),
     SortedDict({'pop1': (33.4484, -112.0740), 'pop2': (37.7749, -122.4194)})),
])
@patch("matplotlib.pyplot.savefig")
def test_plot_clustered_points(mock_savefig, point_coords, popmap, centroids):
    out = "test_output"
    cluster_pops.plot_clustered_points(point_coords, popmap, out, centroids)
    mock_savefig.assert_called_with(out + ".clusteredPoints.pdf")
    plt.clf()


@pytest.mark.parametrize("dat", [
    # Test with a list of distances
    [0.5, 1.2, 0.8, 1.5],

    # Test with an empty list
    [],
])
@patch("matplotlib.pyplot.savefig")
def test_plot_histogram(mock_savefig, dat):
    out = "test_output"
    cluster_pops.plot_histogram(dat, out)
    mock_savefig.assert_called_with(f"{out}.snapDistances.pdf")
    plt.clf()


@pytest.mark.parametrize("popmap, expected", [
    # Typical case
    (SortedDict({'pop1': ['A', 'B'], 'pop2': ['C']}),
     {'A': 'pop1', 'B': 'pop1', 'C': 'pop2'}),

    # Empty input
    (SortedDict({}),
     {}),

    # Case with individual in multiple populations (if applicable)
    (SortedDict({'pop1': ['A', 'B'], 'pop2': ['B', 'C']}),
     {'A': 'pop1', 'B': 'pop2', 'C': 'pop2'}),
])
def test_flatten_popmap(popmap, expected):
    assert cluster_pops.flatten_popmap(popmap) == expected
