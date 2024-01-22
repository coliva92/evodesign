import numpy as np
from typing import Tuple





class Random:

  _rng = None



  @classmethod
  def get_rng(cls) -> np.random.Generator:
    """
    Returns
    -------
    np.random.Generator
        The singleton RNG instance.
    """
    if not cls._rng:
      cls._rng = np.random.default_rng()
    return cls._rng
  


  @classmethod
  def coin_toss(cls, bias: float = 0.5) -> bool:
    """
    Randomly produce a boolean value.

    Parameters
    ----------
    bias : float, optional
        The probability for producing `True`; the default is 0.5. 
        The probability for producing `False` will be computed as `1.0 - bias`.

    Returns
    -------
    bool
        The produced value.
    """
    rng = cls.get_rng()
    return rng.choice([ True, False ], [ bias, 1. - bias ])
