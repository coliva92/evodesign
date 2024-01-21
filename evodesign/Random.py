import numpy as np
from numpy.random import Generator





class Random:

  _rng = None



  @classmethod
  def get_rng(cls) -> Generator:
    if not cls._rng:
      cls._rng = np.random.default_rng()
    return cls._rng
