import numpy as np





def dominates(a: np.array, 
              b: np.array
              ) -> bool:
  return np.all(a >= b) and np.any(a > b)
