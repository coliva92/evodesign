"""Colección de funciones auxiliares para trabajar con secuencias de 
aminoácidos representadas como cadenas de texto donde cada letra corresponde a 
una cadena lateral particular.
"""
import random





"""Las letras válidas para representar los veinte aminoácidos esenciales.
"""
AMINOACIDS = 'ACDEFGHIKLMNPQRSTVWY'



"""La posición ordinal de cada letra que representa un aminoácido, según su 
orden alfabético."""
AMINOACID_ORDINALS = {
  'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8,
  'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 
  'V': 17, 'W': 18, 'Y': 19
}



TEST_SEQUENCE = 'GREETINGSHEREAREELEVENSIDECHAINNAMES' + \
  'ALANINECYSTEINEASPARTATEPHENYLALANINEGLYCINEHISTIDINELYSINEASPARAGINE' + \
  'ARGININESERINEVALINE'



def create_random_sequence(length: int) -> str:
  """Retorna una una cadena de aminoácidos aleatoria de longitud especificada 
  por `length`.
  """
  return ''.join(random.choices(AMINOACIDS, k=length))
