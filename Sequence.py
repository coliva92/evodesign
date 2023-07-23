"""Colección de funciones auxiliares para trabajar con secuencias de 
aminoácidos representadas como cadenas de texto donde cada letra corresponde a 
una cadena lateral particular.
"""
import random





"""Las letras válidas para representar los veinte aminoácidos esenciales.
"""
AMINOACIDS = 'ACDEFGHIKLMNPQRSTVWY'



AMINOACIDS_BY_CATEGORY = {
  'negative': 'DE',
  'positive': 'HKR',
  'unchargedPolar': 'NQSTY',
  'nonpolar': 'ACFGILMPVW'
}



AMINOACID_CATEGORIES = {
  'A': 'nonpolar',
  'C': 'nonpolar',
  'D': 'negative',
  'E': 'negative',
  'F': 'nonpolar',
  'G': 'nonpolar',
  'H': 'positive',
  'I': 'nonpolar',
  'K': 'positive',
  'L': 'nonpolar',
  'M': 'nonpolar',
  'N': 'unchargedPolar',
  'P': 'nonpolar',
  'Q': 'unchargedPolar',
  'R': 'positive',
  'S': 'unchargedPolar',
  'T': 'unchargedPolar',
  'V': 'nonpolar',
  'W': 'nonpolar',
  'Y': 'unchargedPolar'
}



TEST_SEQUENCE = 'GREETINGSHEREAREELEVENSIDECHAINNAMES' + \
  'ALANINECYSTEINEASPARTATEPHENYLALANINEGLYCINEHISTIDINELYSINEASPARAGINE' + \
  'ARGININESERINEVALINE'



def create_random_sequence(length: int) -> str:
  """Retorna una una cadena de aminoácidos aleatoria de longitud especificada 
  por `length`.
  """
  return ''.join(random.choices(AMINOACIDS, k=length))
