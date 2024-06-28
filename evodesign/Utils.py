from typing import List
import pandas as pd
import numpy as np
import numpy.typing as npt





def dominates(a: npt.NDArray[np.float64], 
              b: npt.NDArray[np.float64]
              ) -> bool:
    """
    Returns `True` if vector `a` dominates vector `b`, or `False` 
    otherwise. 
    
    We say vector `a` _dominates_ vector `b` when, for all
    components _x_ in `a` and all components _y_ in `b`, it is the case that
    _x_ >= _y_, and there exists at least one component _z_ in `a` such that 
    _z_ > _y_ for all components _y_.

    Note that if `dominates(a, b)` returns `False`, that does not imply that 
    `dominates(b, a)` will return `True`.
    
    Parameters
    ----------
    a : numpy.typing.NDArray[numpy.float64]
        A vector of real numbers to be compared.
    b : numpy.typing.NDArray[numpy.float64]
        The other vector of real numbers to be compared.

    Returns
    -------
    bool
        Returns `True` if vector `a` dominates over vector `b`, or `False` 
        otherwise.
    """
    return np.all(a >= b) and np.any(a > b)



def is_sorted_before(a: pd.Series,
                     b: pd.Series,
                     columns: List[str],
                     ascending: List[bool]
                     ) -> bool:
    """
    Returns `True` if series `a` must be sorted before series `b`, or `False` 
    otherwise, according to the values of the columns specified by `columns`
    in the specified by `ascending`.

    It is _assumed_ that both `a` and `b` contain the same columns.

    Parameters
    ----------
    a, b : pd.Series
        The series to be compared to determine their relative order.
    columns : List[str]
        The columns which values will be considered to determine the relative
        order of `a` and `b`.
    ascending : List[bool]
        If `True`, it indicates that the values in the corresponding column will 
        be sorted in _ascending_ order. And if `False`, indicates that the values 
        will be sorted in _descending_ order. It is _assumed_ that both `columns`
        and `ascending` have the same length.

    Returns
    -------
    bool
        Returns `True` if series `a` must be sorted before series `b`, or 
        `False` otherwise.
    """
    for k, asc in zip(columns, ascending):
        if a[k] == b[k]: continue
        return a[k] < b[k] if asc else a[k] > b[k]
    return True



def merge(a: pd.DataFrame,
          b: pd.DataFrame,
          columns: List[str],
          ascending: List[bool],
          size: int = 0
          ) -> pd.DataFrame:
    """
    Merges the sorted data frames `a` and `b` into a new sorted data frame.

    It is _assumed_ that both `a` and `b` are sorted and that they both contain
    the same columns. 

    Parameters
    ----------
    a, b : pd.DataFrame
        The data frames to be merged.
    columns : List[str]
        The columns which values will be considered to determine the relative
        order of `a` and `b`.
    ascending : List[bool]
        If `True`, it indicates that the values in the corresponding column will 
        be sorted in _ascending_ order. And if `False`, indicates that the values 
        will be sorted in _descending_ order. It is _assumed_ that both `columns`
        and `ascending` have the same length.
    size : int, optional
        If a value greater than zero is provided, then only the first `size` rows
        in both `a` and `b` will be merged. If the value zero is provided, then
        all the rows in both `a` and `b` will be merged. Default value is zero.

    Returns
    -------
    pd.DataFrame
        A new sorted data frame that contains the same rows as `a` and `b`.
    """
    if not size: size = len(a) + len(b)
    merged = pd.DataFrame(columns=a.columns)
    i, j, n, m = 0, 0, len(a), len(b)
    while i < n and j < m and len(merged) < size:
        row_a, row_b = a.iloc[i], b.iloc[j]
        if is_sorted_before(row_a, row_b, columns, ascending):
            merged = pd.concat([ merged, a.iloc[i].to_frame().T ],
                               ignore_index=True)
            i += 1
        else:
            merged = pd.concat([ merged, b.iloc[j].to_frame().T ],
                               ignore_index=True)
            j += 1
    unmerged = size - len(merged)
    if i < n and unmerged > 0:
        merged = pd.concat([ merged, a.iloc[i:i + unmerged] ], 
                          axis=0, 
                          ignore_index=True)
    if j < m and unmerged > 0:
        merged = pd.concat([ merged, b.iloc[j:j + unmerged] ], 
                          axis=0, 
                          ignore_index=True)
    return merged
