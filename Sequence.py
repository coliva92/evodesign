import random





AMINOACIDS = 'ACDEFGHIKLMNPQRSTVWY'



def random_sequence(length: int) -> str:
  return ''.join(random.choices(AMINOACIDS, k=length))



def switch_residue(old_letter: str) -> str:
  new_letter = random.choice(AMINOACIDS)
  while new_letter == old_letter:
    new_letter = random.choice(AMINOACIDS)
  return new_letter
