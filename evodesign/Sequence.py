from .Random import Random





AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'



def create_random(length: int) -> str:
  """
  Generates randomly a sequence of amino acids.

  Parameters
  ----------
  length : int
      The number of residues in the produced sequence.

  Returns
  -------
  str
      The generated sequence. Each residue is represented as a
      single letter corresponding to one of the 20 essential amino acids.
  """
  rng = Random.get_rng()
  indices = rng.choice(len(AMINO_ACIDS), size=length)
  return ''.join([ AMINO_ACIDS[i] for i in indices ])



def switch_residue(residue: str) -> str:
  """
  Randomly produce a single amino acid letter that's different from the one
  specified by `residue`.

  Parameters
  ----------
  residue : str
      A single letter representing one of the 20 essential amino acids.

  Returns
  -------
  str
      A different letter, which also represents one of the 20 essential amino 
      acids.
  """
  rng = Random.get_rng()
  while True:
    new_residue = AMINO_ACIDS[rng.choice(len(AMINO_ACIDS))]
    if new_residue != residue:
      return new_residue
