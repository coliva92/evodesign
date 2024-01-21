from .Random import Random





AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'



def random_sequence(length: int) -> str:
  rng = Random.get_rng()
  indices = rng.choice(len(AMINO_ACIDS), size=length)
  return ''.join([ AMINO_ACIDS[i] for i in indices ])



def switch_residue(residue: str) -> str:
  rng = Random.get_rng()
  while True:
    new_residue = AMINO_ACIDS[rng.choice(len(AMINO_ACIDS))]
    if new_residue != residue:
      return new_residue
