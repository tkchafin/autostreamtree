import numpy as np
import pytest
from sortedcontainers import SortedDict
from unittest.mock import patch
import autostreamtree.genetic_distances as gendist


dist_options = ["PDIST", "JOST", "HARMD", "FST", "LINFST", "CHORD"]


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
    # Per locus (Jost 2008, n = 2 subpopulations): 0.4, 0 and 1; their
    # harmonic mean (without the 0) is 4/7
    expected_fst_values = np.array([[0.0, 4 / 7], [4 / 7, 0.0]])
    np.testing.assert_allclose(fst_result, expected_fst_values, atol=0.01)


def test_jost_calculation_v5(mock_snps, mock_popmap):
    fst_result = gendist.get_pop_genmat(dist="HARMD", indmat=mock_indmat,
                                        popmap=mock_popmap, dat=mock_dat,
                                        seqs=mock_snps, loc_agg="HARM",
                                        global_het=True)
    # Per locus (Jost 2008, n = 2 subpopulations): 0.4, 0 and 1; their
    # harmonic mean (without the 0) is 4/7
    expected_fst_values = np.array([[0.0, 4 / 7], [4 / 7, 0.0]])
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


def test_pop_pdist(mock_indmat, mock_popmap, mock_seqs, mock_dat):
    # will fil in expected later 
    expected_matrix = np.array([
        [0.0, 3.5],
        [3.5, 0.0]
    ])
    result_matrix = gendist.get_pop_genmat(dist="PDIST", indmat=mock_indmat,
                                        popmap=mock_popmap, dat=mock_dat,
                                        seqs=mock_seqs, pop_agg="ARITH")
    np.testing.assert_allclose(result_matrix, expected_matrix, atol=1e-3)


def test_jost_pairwise_with_more_populations(mock_snps):
    # Pairwise D uses Jost's n = 2 (the two populations compared), whatever
    # the number of populations: two populations fixed for different
    # alleles have D = 1, and identical ones D = 0
    snps = {'Ind1': ['A/A'], 'Ind2': ['A/A'], 'Ind3': ['T/T'],
            'Ind4': ['T/T'], 'Ind5': ['A/A'], 'Ind6': ['A/A']}
    popmap = SortedDict({'Pop1': ['Ind1', 'Ind2'], 'Pop2': ['Ind3', 'Ind4'],
                         'Pop3': ['Ind5', 'Ind6']})
    d = gendist.get_pop_genmat(dist="JOST", indmat=None, popmap=popmap,
                               dat=None, seqs=snps, loc_agg="ARITH",
                               global_het=True)
    np.testing.assert_allclose(d, [[0, 1, 0], [1, 0, 1], [0, 1, 0]],
                               atol=1e-12)


def test_weir_cockerham_unequal_sample_sizes():
    # Weir & Cockerham's (1984) theta for populations of unequal sizes (2
    # and 6 individuals) equals scikit-allel's weir_cockerham_fst (1.3.13)
    # on the same genotypes: 0.4924657534246576 (with the old n_c, which
    # halved the sample size correction, it was underestimated)
    gts = {
        "a1": ["A/A", "A/T", "T/T", "A/A", "A/T"],
        "a2": ["A/T", "A/A", "T/T", "A/A", "T/T"],
        "b1": ["T/T", "A/T", "A/A", "A/T", "A/A"],
        "b2": ["T/T", "T/T", "A/T", "A/A", "A/A"],
        "b3": ["A/T", "T/T", "A/A", "T/T", "A/T"],
        "b4": ["T/T", "A/T", "A/A", "A/T", "A/A"],
        "b5": ["T/T", "T/T", "A/T", "A/A", "A/A"],
        "b6": ["A/T", "T/T", "A/A", "T/T", "A/T"],
    }
    popmap = SortedDict({"A": ["a1", "a2"],
                         "B": ["b1", "b2", "b3", "b4", "b5", "b6"]})
    fst = gendist.get_pop_genmat("FST", None, popmap, None, gts)
    assert fst[0, 1] == pytest.approx(0.4924657534246576, abs=1e-12)
    lin = gendist.get_pop_genmat("LINFST", None, popmap, None, gts)
    assert lin[0, 1] == pytest.approx(
        0.4924657534246576 / (1 - 0.4924657534246576), abs=1e-12)


def test_pop_genmat_warns_only_for_negative_distances(capsys):
    # Negative FST is set to 0 with a warning; no warning otherwise
    popmap = SortedDict({"a": ["a1", "a2"], "b": ["b1", "b2"],
                         "c": ["c1", "c2"]})
    same = {"a1": ["A/A"], "a2": ["A/T"], "b1": ["A/A"], "b2": ["A/T"],
            "c1": ["T/T"], "c2": ["T/T"]}
    mat = gendist.get_pop_genmat("FST", None, popmap, None, same)
    assert mat[0, 1] == 0.0 and mat[0, 2] > 0
    assert "Coercing" in capsys.readouterr().out
    apart = {"a1": ["A/A"], "a2": ["A/A"], "b1": ["T/T"], "b2": ["T/T"],
             "c1": ["A/T"], "c2": ["A/T"]}
    mat = gendist.get_pop_genmat("FST", None, popmap, None, apart)
    assert (mat[np.triu_indices(3, 1)] > 0).all()
    assert "Coercing" not in capsys.readouterr().out
