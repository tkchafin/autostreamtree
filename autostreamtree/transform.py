import sys
import itertools
import numpy as np
import pandas as pd

from typing import List, Tuple, Dict, Any, Union, Optional

"""
Transformations from Peterman et al. (2018) ResistanceGA package
Full Citation:
Peterman, W. E. 2018. ResistanceGA: An R package for the optimization 
of resistance surfaces using genetic algorithms. Methods in Ecology 
and Evolution doi:10.1111/2041-210X.12984
GitHub:
https://github.com/wpeterman/ResistanceGA/tree/3be50a6fa5515fe23953fca39135c5b340415607
"""

def rescale_cols(df: pd.DataFrame, m: float, M: float) -> pd.DataFrame:
    """
    Rescale all columns in a Pandas DataFrame to a new range between m and M.

    Args:
        df (pd.DataFrame): The DataFrame to rescale.
        m (float): The minimum value of the new range.
        M (float): The maximum value of the new range.

    Returns:
        rescaled_df (pd.DataFrame): A new DataFrame with all columns rescaled to the new range.
    """
    df -= df.min()
    df /= df.max()
    rescaled_df = (df * (M - m)) + m
    return rescaled_df

def ricker(dat: pd.Series, shape: float, ceiling: float) -> pd.Series:
    """
    Compute Ricker function values for a given range of data points.

    Args:
        dat (pd.Series): A pandas Series object containing the data points to compute Ricker values for.
        shape (float): The "shape" parameter of the Ricker function.
        ceiling (float): The "ceiling" parameter of the Ricker function.

    Returns:
        r (pd.Series): A pandas Series object containing the Ricker values for the input data points.
    """
    r = ceiling * dat * np.exp(-1 * dat / shape) + 1
    return r

def inv_ricker(dat: pd.Series, shape: float, ceiling: float) -> pd.Series:
    """
    Compute the inverse of the Ricker function for a given range of data points.

    Args:
        dat (pd.Series): A pandas Series object containing the data points to compute inverse Ricker values for.
        shape (float): The "shape" parameter of the Ricker function.
        ceiling (float): The "ceiling" parameter of the Ricker function.

    Returns:
        inv_r (pd.Series): A pandas Series object containing the inverse Ricker values for the input data points.
    """
    inv_r = (-1 * ceiling) * dat * np.exp(-1 * dat / shape) - 1
    return inv_r

def rev_inv_ricker(dat: pd.Series, shape: float, ceiling: float) -> pd.Series:
    """
    Compute the reverse inverse of the Ricker function for a given range of data points.

    Args:
        dat (pd.Series): A pandas Series object containing the data points to compute reverse inverse Ricker values for.
        shape (float): The "shape" parameter of the Ricker function.
        ceiling (float): The "ceiling" parameter of the Ricker function.

    Returns:
        rescaled_inv_r (pd.Series): A pandas Series object containing the reverse inverse Ricker values for the input data points, rescaled to a new range between the minimum and maximum values of the original values.
    """
    inv_r = inv_ricker(dat, shape, ceiling)
    rescaled_inv_r = rescale_cols((-1 * inv_r), min(inv_r), max(inv_r))
    return rescaled_inv_r

def inv_ricker(dat: pd.Series, shape: float, ceiling: float) -> pd.Series:
    """
    Compute the reverse of the Ricker function for a given range of data points.

    Args:
        dat (pd.Series): A pandas Series object containing the data points to compute reverse Ricker values for.
        shape (float): The "shape" parameter of the Ricker function.
        ceiling (float): The "ceiling" parameter of the Ricker function.

    Returns:
        rescaled_r (pd.Series): A pandas Series object containing the reverse Ricker values for the input data points, rescaled to a new range between the minimum and maximum values of the original values.
    """
    d = rescale_cols((-1 * dat), min(dat), max(dat))
    rescaled_r = ricker(d, shape, ceiling)
    return rescaled_r

def monomolecular(dat: pd.Series, shape: float, ceiling: float) -> pd.Series:
    """
    Compute monomolecular function values for a given range of data points.

    Args:
        dat (pd.Series): A pandas Series object containing the data points to compute monomolecular values for.
        shape (float): The "shape" parameter of the monomolecular function.
        ceiling (float): The "ceiling" parameter of the monomolecular function.

    Returns:
        m (pd.Series): A pandas Series object containing the monomolecular values for the input data points.
    """
    m = ceiling * (1 - np.exp(-1 * dat / shape)) + 1
    return m
    
def inv_monomolecular(dat: pd.Series, shape: float, ceiling: float) -> pd.Series:
    """
    Compute the inverse of the monomolecular function for a given range of data points.

    Args:
        dat (pd.Series): A pandas Series object containing the data points to compute inverse monomolecular values for.
        shape (float): The "shape" parameter of the monomolecular function.
        ceiling (float): The "ceiling" parameter of the monomolecular function.

    Returns:
        inv_m (pd.Series): A pandas Series object containing the inverse monomolecular values for the input data points.
    """
    d = ceiling * np.exp(-1 * dat / shape)
    inv_m = (d - min(d)) + 1
    return inv_m

def rev_inv_monomolecular(dat: pd.Series, shape: float, ceiling: float) -> pd.Series:
    """
    Compute the reverse of the inverse monomolecular function for a given range of data points.

    Args:
        dat (pd.Series): A pandas Series object containing the data points to compute reverse inverse monomolecular values for.
        shape (float): The "shape" parameter of the inverse monomolecular function.
        ceiling (float): The "ceiling" parameter of the inverse monomolecular function.

    Returns:
        rescaled_inv_m (pd.Series): A pandas Series object containing the reverse inverse monomolecular values for the input data points, rescaled to a new range between the minimum and maximum values of the original values.
    """
    d = rescale_cols((-1 * dat), min(dat), max(dat))
    rescaled_inv_m = inv_monomolecular(d, shape, ceiling)
    return rescaled_inv_m

def rev_monomolecular(dat: pd.Series, shape: float, ceiling: float) -> pd.Series:
    """
    Compute the reverse of the monomolecular function for a given range of data points.

    Args:
        dat (pd.Series): A pandas Series object containing the data points to compute reverse monomolecular values for.
        shape (float): The "shape" parameter of the monomolecular function.
        ceiling (float): The "ceiling" parameter of the monomolecular function.

    Returns:
        rescaled_m (pd.Series): A pandas Series object containing the reverse monomolecular values for the input data points, rescaled to a new range between the minimum and maximum values of the original values.
    """
    d = rescale_cols((-1 * dat), min(dat), max(dat))
    rescaled_m = monomolecular(d, shape, ceiling)
    return rescaled_m

    