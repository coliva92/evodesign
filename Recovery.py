"""Colección de funciones auxiliares que se encargan de almacenar y restaurar 
los datos respaldados del algoritmo evolutivo y de las poblaciones procesadas 
durante la ejecución de ese mismo algoritmo.
"""
import sys
import json
from .Individual import Individual
from typing import List, Tuple
from .Algorithm import Algorithm
from .Prediction import *
from .Fitness import *
from .GA import *
from .GA.Selection import *
from .GA.Recombination import *
from .GA.Mutation import *
from .GA.Replacement import *





def load_dict_from_json(filename: str) -> dict:
  """Carga el diccionario contenido en el archivo JSON especificado por 
  `filename`.
  """
  with open(filename, 'rt', encoding='utf-8') as the_file:
    data = json.load(the_file)
  return data



def create_algorithm_from_settings(settings: dict) -> Algorithm:
  """Usa los datos contenidos en el diccionario especificado por `memento` para 
  construir una instancia del algoritmo evolutivo indicado en ese mismo 
  diccionario.
  """
  algorithm_class = getattr(sys.modules[__name__], settings['algorithmType'])
  i = settings['algorithmType'].find('_')
  algorithm_name = settings['algorithmType'][:i]
  algorithm_params, operation_params = {}, {}
  other_algorithm_params = {}
  for key, value in settings['algorithmParams'].items():
    if key.endswith('Params'):
      name = key.replace('Params', '')
      operation_params[name] = value
      continue
    if type(value) == str and value.startswith(algorithm_name):
      algorithm_params[key] = getattr(sys.modules[__name__], value)
      continue
    if type(value) == str and value.startswith('Predictor'):
      other_algorithm_params[key] = getattr(sys.modules[__name__], value)()
      continue
    if type(value) == str and value.startswith('Fitness'):
      other_algorithm_params[key] = getattr(sys.modules[__name__], value)()
      continue
    other_algorithm_params[key] = value
  for key, value in algorithm_params.items():
    algorithm_params[key] = value(**operation_params[key])
  for key, value in other_algorithm_params.items():
    algorithm_params[key] = value
  algorithm = algorithm_class(**algorithm_params)
  if '__savedPopulations' in settings:
    algorithm.append_and_cache_memento(settings['__savedPopulations'])
  else:
    algorithm.append_and_cache_memento()
  return algorithm



def load_latest_population_from_settings(settings: dict
                                         ) -> Tuple[int, List[Individual]]:
  """Carga la colección de individuos del archivo JSON más reciente descrito en 
  el diccionario de datos de respaldo especificado por `memento`.
  """
  if '__savedPopulations' not in settings or \
      len(settings['__savedPopulations']) == 0:
    return 0, []
  filename = settings['__savedPopulations'][-1]
  with open(filename, 'rt', encoding='utf-8') as json_file:
    individuals_data = json.load(json_file)
  iterationId = len(settings['__savedPopulations']) - 1
  population = [ Individual(**params) for params in individuals_data ]
  return iterationId, population
