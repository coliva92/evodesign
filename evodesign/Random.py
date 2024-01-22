import numpy as np
from numpy.random import Generator
from typing import Tuple





class Random:

  _rng = None



  @classmethod
  def get_rng(cls) -> Generator:
    if not cls._rng:
      cls._rng = np.random.default_rng()
    return cls._rng
  


  @classmethod
  def coin_toss(cls, weights: Tuple[float] = (0.5, 0.5)) -> bool:
    rng = cls.get_rng()
    return rng.choice([ True, False ], weights)
