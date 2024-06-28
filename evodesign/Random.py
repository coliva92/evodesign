import numpy as np



_rng = None





def generator() -> np.random.Generator:
    """
    Returns a singleton RNG instance.

    Returns
    -------
    numpy.random.Generator
        The singleton RNG instance.
    """
    global _rng
    if not _rng: _rng = np.random.default_rng()
    return _rng
  


def coin_toss(bias: float = 0.5) -> bool:
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
        A randomly produced value.
    """
    rng = generator()
    return rng.choice([ True, False ], p=[ bias, 1. - bias ])
