import numpy as np
import scipy.stats

import autostreamtree.aggregators as agg


def test_aggregate_dist():
    # Test data
    data = [1, 2, 3, 4, 5]

    # Expected results
    expected_harm = scipy.stats.hmean(data)
    expected_arith = np.mean(data)
    expected_geom = scipy.stats.mstats.gmean(data)
    expected_median = np.median(data)
    expected_max = np.max(data)
    expected_min = np.min(data)
    expected_adjharm = agg.adjusted_harmonic_mean(data)
    expected_sd = np.std(data)
    expected_var = np.var(data)

    # Assertions
    assert agg.aggregate_dist("HARM", data) == expected_harm
    assert agg.aggregate_dist("ARITH", data) == expected_arith
    assert agg.aggregate_dist("GEOM", data) == expected_geom
    assert agg.aggregate_dist("MEDIAN", data) == expected_median
    assert agg.aggregate_dist("MAX", data) == expected_max
    assert agg.aggregate_dist("MIN", data) == expected_min
    assert agg.aggregate_dist("ADJHARM", data) == expected_adjharm
    assert agg.aggregate_dist("SD", data) == expected_sd
    assert agg.aggregate_dist("VAR", data) == expected_var

    # Test with zero and negative values for harmonic mean
    data_with_zero = [0, 1, 2, 3, 4, 5]
    assert agg.aggregate_dist("HARM", data_with_zero) == expected_harm
