import random





AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'



def random_sequence(length: int) -> str:
  return ''.join(random.choices(AMINO_ACIDS, k=length))



def switch_residue(residue: str) -> str:
  while True:
    new_residue = random.choice(AMINO_ACIDS)
    if new_residue != residue:
      return new_residue
