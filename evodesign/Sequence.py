from numpy.random import Generator





AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'



def random_sequence(rng: Generator, 
                    length: int
                    ) -> str:
  indices = rng.choice(len(AMINO_ACIDS), size=length)
  return ''.join([ AMINO_ACIDS[i] for i in indices ])



def switch_residue(rng: Generator,
                   residue: str
                   ) -> str:
  while True:
    new_residue = AMINO_ACIDS[rng.choice(len(AMINO_ACIDS))]
    if new_residue != residue:
      return new_residue
