import random





"""
Las letras válidas para representar los veinte aminoácidos esenciales.
"""
AMINOACIDS = 'ACDEFGHIKLMNPQRSTVWY'



def create_random_sequence(length: int) -> str:
  """
  Genera de manera aleatoria una cadena de aminoácidos de longitud especificada 
  por `length`.
  """
  return ''.join(random.choices(AMINOACIDS, k=length))
