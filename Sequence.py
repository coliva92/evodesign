import random





AMINOACIDS = 'ACDEFGHIKLMNPQRSTVWY'



def create_random_sequence(length: int) -> str:
  return ''.join(random.choices(AMINOACIDS, k=length))
