"""Colección de funciones auxiliares para trabajar con secuencias de 
aminoácidos representadas como cadenas de texto donde cada letra corresponde a 
una cadena lateral particular.
"""
import random





"""Las letras válidas para representar los veinte aminoácidos esenciales.
"""
AMINOACIDS = 'ACDEFGHIKLMNPQRSTVWY'



_TEST_SEQUENCE = 'GREETINGSHEREAREELEVENSIDECHAINNAMES' + \
  'ALANINECYSTEINEASPARTATEPHENYLALANINEGLYCINEHISTIDINELYSINEASPARAGINE' + \
  'ARGININESERINEVALINE'



def create_random_sequence(length: int) -> str:
  """Retorna una una cadena de aminoácidos aleatoria de longitud especificada 
  por `length`.
  """
  return ''.join(random.choices(AMINOACIDS, k=length))
