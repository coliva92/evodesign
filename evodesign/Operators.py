from typing import List
import pandas as pd
import numpy as np





def dominates(a: np.array, 
              b: np.array
              ) -> bool:
  return np.all(a >= b) and np.any(a > b)



def comes_first(a: pd.Series,
                b: pd.Series,
                sortColumns: List[str],
                sortAcending: List[bool]
                ) -> bool:
  for col, ascending in zip(sortColumns, sortAcending):
    if a[col] == b[col]:
      continue
    return a[col] < b[col] if ascending else a[col] > b[col]
  return True



def merge(a: pd.DataFrame,
          b: pd.DataFrame,
          sortColumns: List[str],
          sortAcending: List[bool],
          size: int = 0
          ) -> pd.DataFrame:
  if not size:
    size = len(a) + len(b)
  
  # we assume that both a and b are sorted and have the same columns
  merged = pd.DataFrame(columns=a.columns)
  i, j, n, m = 0, 0, len(a), len(b)
  while i < n and j < m and len(merged) < size:
    row_a, row_b = a.iloc[i], b.iloc[j]
    if comes_first(row_a, row_b, sortColumns, sortAcending):
      merged = pd.concat([ merged, a.iloc[i].to_frame().T ],
                         ignore_index=True)
      i += 1
    else:
      merged = pd.concat([ merged, b.iloc[j].to_frame().T ],
                         ignore_index=True)
      j += 1
  missing = size - len(merged)
  if i < n and missing > 0:
    merged = pd.concat([ merged, a.iloc[i:i + missing] ], 
                       axis=0, 
                       ignore_index=True)
  if j < m and missing > 0:
    merged = pd.concat([ merged, b.iloc[j:j + missing] ], 
                       axis=0, 
                       ignore_index=True)
  return merged
