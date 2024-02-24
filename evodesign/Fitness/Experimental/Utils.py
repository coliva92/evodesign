_CYC_MEAN = 1.3248119
_CYC_STD_DEV = 0.10498072





def normalize(x: float) -> float:
  return 1.0 / (1.0 + x)



def cyclization_z_score(cyclization: float) -> float:
  return (cyclization - _CYC_MEAN) / _CYC_STD_DEV



def normalize_plddt(plddt: float) -> float:
  return plddt / 100 if plddt > 1.0 else plddt
