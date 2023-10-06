import random





AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
AMINO_ACIDS_SET = { 
  'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 
  'S', 'T', 'V', 'W', 'Y' 
}



def random_sequence(length: int) -> str:
  return ''.join(random.choices(AMINO_ACIDS, k=length))



def switch_residue(old_letter: str) -> str:
  new_letter = random.choice(AMINO_ACIDS)
  while new_letter == old_letter:
    new_letter = random.choice(AMINO_ACIDS)
  return new_letter
