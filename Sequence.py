import random





AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'



def random_sequence(length: int) -> str:
  return ''.join(random.choices(AMINO_ACIDS, k=length))



def switch_residue(old_letter: str) -> str:
  new_letter = random.choice(AMINO_ACIDS)
  while new_letter == old_letter:
    new_letter = random.choice(AMINO_ACIDS)
  return new_letter
