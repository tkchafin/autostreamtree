import numpy as np
import pytest
from sortedcontainers import SortedDict
from unittest.mock import patch
import autostreamtree.genetic_distances as gendist


dist_options = ["JC69", "PDIST", "JOST", "HARMD", "FST", "LINFST", "CHORD"]


def mock_agg_dist(*args, **kwargs):
    # Mock behavior of agg.aggregate_dist
    return 0.5  # Return a constant value or based on input for testing


@pytest.fixture
def mock_popmap():
    # Example mock population map
    return SortedDict(
        {'Pop1': ['Ind1', 'Ind2'], 'Pop2': ['Ind3', 'Ind4']}
    )


@pytest.fixture
def mock_dat():
    # Example mock data
    return ['Ind1', 'Ind2', 'Ind3', 'Ind4']


@pytest.fixture
def mock_points():
    return SortedDict(
        {'Ind1': (0, 0), 'Ind2': (1, 0), 'Ind3': (0, 1), 'Ind4': (1, 1)}
    )


@pytest.fixture
def mock_snps():
    return {'Ind1': ['A/A', 'A/T', 'A/A'], 'Ind2': ['A/T', 'A/T', 'A/A'],
            'Ind3': ['A/T', 'A/T', 'T/T'], 'Ind4': ['T/T', 'A/T', 'T/T']}


@pytest.fixture
def mock_seqs():
    return {
        'Ind1': ['AAAAAAA', 'AAAAAAT', 'AAAAAAA'],
        'Ind2': ['AAAAAAT', 'AAAAAAT', 'AAAAAAA'],
        'Ind3': ['AAAAAAT', 'AAAAAAT', 'TTTTTTT'],
        'Ind4': ['TTTTTTT', 'AAAAAAT', 'TTTTTTT']
    }


@pytest.fixture
def mock_seqs_het():
    return {
        'Ind1': ['AAAAAAA/AAAAAAA', 'AAAAAAA/AAAAAAA'],
        'Ind2': ['AAAAAAA/AAAAAAT', 'AAAAAAA/AAAAAAT'],
        'Ind3': ['AAAAAAA/TTTTTTT', 'AAAAAAA/TTTTTTT'],
        'Ind4': ['TTTTTTT/TTTTTTT', 'TTTTTTT/TTTTTTT']
    }


@pytest.fixture
def mock_indmat():
    # Example mock individual distance matrix
    return np.array([[0, 1, 2, 3], [1, 0, 4, 5], [2, 4, 0, 6], [3, 5, 6, 0]])


@patch('autostreamtree.genetic_distances.agg.aggregate_dist',
       side_effect=mock_agg_dist)
def test_get_pop_genmat_unsupported_option(mock_agg, mock_popmap, mock_dat,
                                           mock_snps, mock_indmat):
    with pytest.raises(ValueError):
        gendist.get_pop_genmat("UNSUPPORTED_DIST", mock_indmat, mock_popmap,
                               mock_dat, mock_snps)


@pytest.mark.parametrize("dist", dist_options)
@patch('autostreamtree.genetic_distances.agg.aggregate_dist',
       side_effect=mock_agg_dist)
def test_get_pop_genmat_valid_options(mock_agg, mock_popmap, mock_dat,
                                      mock_snps, mock_indmat, dist):
    result = gendist.get_pop_genmat(dist, mock_indmat, mock_popmap, mock_dat,
                                    mock_snps)

    assert isinstance(result, np.ndarray)
    assert result.shape == (2, 2)
    assert np.all(np.diag(result) == 0)


def test_fst_calculation(mock_snps, mock_popmap):
    fst_result = gendist.get_pop_genmat("FST", mock_indmat, mock_popmap,
                                        mock_dat, mock_snps)
    expected_fst_values = np.array([[0.0, 0.5], [0.5, 0.0]])
    np.testing.assert_allclose(fst_result, expected_fst_values, atol=0.05)


def test_linfst_calculation(mock_snps, mock_popmap):
    fst_result = gendist.get_pop_genmat("LINFST", mock_indmat, mock_popmap,
                                        mock_dat, mock_snps)
    expected_fst_values = np.array([[0.0, 1.125], [1.125, 0.0]])
    np.testing.assert_allclose(fst_result, expected_fst_values, atol=0.05)


def test_chord_calculation(mock_snps, mock_popmap):
    fst_result = gendist.get_pop_genmat("CHORD", mock_indmat, mock_popmap,
                                        mock_dat, mock_snps)
    expected_fst_values = np.array([[0.0, np.sqrt(2.5)], [np.sqrt(2.5), 0.0]])
    np.testing.assert_allclose(fst_result, expected_fst_values, atol=0.01)


def test_jost_calculation_v1(mock_snps, mock_popmap):
    fst_result = gendist.get_pop_genmat(dist="JOST", indmat=mock_indmat,
                                        popmap=mock_popmap, dat=mock_dat,
                                        seqs=mock_snps, loc_agg="ARITH",
                                        global_het=True)
    expected_fst_values = np.array([[0.0, 0.588], [0.588, 0.0]])
    np.testing.assert_allclose(fst_result, expected_fst_values, atol=0.01)


def test_jost_calculation_v2(mock_snps, mock_popmap):
    fst_result = gendist.get_pop_genmat(dist="JOST", indmat=mock_indmat,
                                        popmap=mock_popmap, dat=mock_dat,
                                        seqs=mock_snps, loc_agg="HARM",
                                        global_het=True)
    expected_fst_values = np.array([[0.0, 0.25], [0.25, 0.0]])
    np.testing.assert_allclose(fst_result, expected_fst_values, atol=0.01)


def test_jost_calculation_v3(mock_snps, mock_popmap):
    fst_result = gendist.get_pop_genmat(dist="JOST", indmat=mock_indmat,
                                        popmap=mock_popmap, dat=mock_dat,
                                        seqs=mock_snps, loc_agg="ADJHARM",
                                        global_het=True)
    expected_fst_values = np.array([[0.0, 0.88], [0.88, 0.0]])
    np.testing.assert_allclose(fst_result, expected_fst_values, atol=0.01)


def test_jost_calculation_v4(mock_snps, mock_popmap):
    fst_result = gendist.get_pop_genmat(dist="HARMD", indmat=mock_indmat,
                                        popmap=mock_popmap, dat=mock_dat,
                                        seqs=mock_snps, loc_agg="ARITH",
                                        global_het=True)
    expected_fst_values = np.array([[0.0, 0.326], [0.326, 0.0]])
    np.testing.assert_allclose(fst_result, expected_fst_values, atol=0.01)


def test_jost_calculation_v5(mock_snps, mock_popmap):
    fst_result = gendist.get_pop_genmat(dist="HARMD", indmat=mock_indmat,
                                        popmap=mock_popmap, dat=mock_dat,
                                        seqs=mock_snps, loc_agg="HARM",
                                        global_het=True)
    expected_fst_values = np.array([[0.0, 0.326], [0.326, 0.0]])
    np.testing.assert_allclose(fst_result, expected_fst_values, atol=0.01)


def test_get_genmat_pdist(mock_seqs, mock_points):
    expected_matrix = np.array([
        [0.0, 1/21, (1/7 + 1)/3, (1 + 1)/3],
        [1/21, 0.0, 1/3, (6/7 + 1)/3],
        [(1/7 + 1)/3, 1/3, 0.0, (6/7)/3],
        [(1 + 1)/3, (6/7 + 1)/3, (6/7)/3, 0.0]
    ])
    result_matrix = gendist.get_genmat("PDIST", mock_points, mock_seqs,
                                       ploidy=2, het=False, loc_agg="ARITH")
    np.testing.assert_allclose(result_matrix, expected_matrix, atol=1e-6)


def test_get_genmat_pdist_het(mock_seqs_het, mock_points):
    expected_matrix = np.array([
        [0.0, 0.071429, 0.5, 1.0],
        [0.071429, 0.0, 0.428571, 0.928571],
        [0.5, 0.428571, 0.0, 0.5],
        [1.0, 0.928571, 0.5, 0.0]
    ])
    result_matrix = gendist.get_genmat("PDIST", mock_points, mock_seqs_het,
                                       ploidy=2, het=True, loc_agg="ARITH")
    np.testing.assert_allclose(result_matrix, expected_matrix, atol=1e-3)


def test_get_genmat_jc69_het(mock_seqs_het, mock_points):
    expected_matrix = np.array([
        [0.0, 0.075, 0.8239, 1.0],
        [0.071429, 0.0, 0.428571, 0.928571],
        [0.5, 0.428571, 0.0, 0.5],
        [1.0, 0.928571, 0.5, 0.0]
    ])
    result_matrix = gendist.get_genmat("JC69", mock_points, mock_seqs_het,
                                       ploidy=2, het=True, loc_agg="ARITH")
    np.testing.assert_allclose(result_matrix, expected_matrix, atol=1e-3)