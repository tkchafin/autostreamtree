import os
import sys
import numpy as np
import scipy


def aggregate_dist(method, stuff):
    """
    Aggregates the given array according to the specified method.

    Args:
        method (str): The aggregation method to use (e.g. "HARM", "ARITH", "GEOM", etc.).
        stuff (array-like): The input array to be aggregated.

    Returns:
        float: The aggregated result.
    """
    if method == "HARM":
        try:
            stuff = [x for x in stuff if x > 0]
            return scipy.stats.hmean(stuff) if len(stuff) > 0 else 0.0
        except ValueError as e:
            print(e)
            print("ERROR (DivideByZero): Harmonic mean cannot be calculated using a zero distance. Try recomputing using the \"ADJHARM\" option.")
            print("")
            sys.exit(1)
    elif method == "ARITH":
        return(np.mean(stuff))
    elif method == "GEOM":
        return(scipy.stats.mstats.gmean(stuff))
    elif method == "MEDIAN":
        return(np.median(stuff))
    elif method == "MAX":
        return(np.max(stuff))
    elif method == "MIN":
        return(np.min(stuff))
    elif method == "ADJHARM":
        return(adjusted_harmonic_mean(stuff))
    elif method == "SD":
        return(np.std(stuff))
    elif method == "VAR":
        return(np.var(stuff))

def adjusted_harmonic_mean(stuff):
    """
    Computes an adjusted harmonic mean that is corrected for non-positive values.

    Args:
        stuff (array-like): The input array.

    Returns:
        float: The adjusted harmonic mean.
    """
    s=np.array(stuff)
    vals = s[s>0.0]
    bads = s[s<=0.0]
    mu = (1.0 / (np.sum([1.0/x for x in vals]) / (len(vals)-len(bads)))) * ((len(vals)-len(bads))/len(vals))
    return(mu)
