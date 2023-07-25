"""Colección de funciones auxiliares que se encargan de almacenar y restaurar 
los datos respaldados del algoritmo evolutivo y de las poblaciones procesadas 
durante la ejecución de ese mismo algoritmo.
"""
import sys
import json
from .Population import Individual
from typing import List, Tuple
from .Algorithm import Algorithm
from .Prediction import *
from .Fitness import *
from .GA import *
from .GA.Selection import *
from .GA.Recombination import *
from .GA.Mutation import *
from .GA.Replacement import *





def load_memento_from_json_file(filename: str) -> dict:
  """Carga el diccionario contenido en el archivo JSON especificado por 
  `filename`.
  """
  with open(filename, 'rt', encoding='utf-8') as the_file:
    memento = json.load(the_file)
  return memento



def load_algorithm_from_memento(memento: dict) -> Algorithm:
  """Usa los datos contenidos en el diccionario especificado por `memento` para 
  construir una instancia del algoritmo evolutivo indicado en ese mismo 
  diccionario.
  """
  algorithm_class = getattr(sys.modules[__name__], memento['algorithmType'])
  algorithm_name = memento['algorithmType'][:memento['algorithmType'].find('_')]
  algorithm_params, operation_params = {}, {}
  other_algorithm_params = {}
  for key, value in memento['algorithmParams'].items():
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
  if '__savedPopulations' in memento:
    algorithm.append_and_cache_memento(memento['__savedPopulations'])
  else:
    algorithm.append_and_cache_memento()
  return algorithm



def load_population_from_memento(memento: dict) -> Tuple[int, List[Individual]]:
  """Carga la colección de individuos del archivo JSON más reciente descrito en 
  el diccionario de datos de respaldo especificado por `memento`.
  """
  if '__savedPopulations' not in memento or \
      len(memento['__savedPopulations']) == 0:
    return 0, []
  with open(memento['__savedPopulations'][-1], 'rt', encoding='utf-8') as the_file:
    population_memento = json.load(the_file)
  iterationId = len(memento['__savedPopulations']) - 1
  population = [ Individual(**params) for params in population_memento ]
  return iterationId, population
