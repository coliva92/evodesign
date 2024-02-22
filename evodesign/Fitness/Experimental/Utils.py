from scipy.stats import norm

_MEAN = 1.3248119
_STANDARD_DEV = 0.10498072





def normalize_rmsd(rmsd: float) -> float:
  return 1.0 / (1.0 + rmsd)



def cyclization_probability(cyclization: float) -> float:
  z_score = (cyclization - _MEAN) / _STANDARD_DEV
  return norm.cdf(z_score)
