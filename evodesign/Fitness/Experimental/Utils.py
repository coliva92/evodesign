from scipy.stats import norm

_CYC_MEAN = 1.3248119
_CYC_STD_DEV = 0.10498072





def normalize_rmsd(rmsd: float) -> float:
  return 1.0 / (1.0 + rmsd)



def cyclization_probability(cyclization: float) -> float:
  z_score = (cyclization - _CYC_MEAN) / _CYC_STD_DEV
  return norm.cdf(z_score)



def normalize_cyclization(cyclization: float) -> float:
  return 1.0 / (1.0 + abs(_CYC_MEAN - cyclization))



def normalize_plddt(plddt: float) -> float:
  return plddt / 100 if plddt > 1.0 else plddt
