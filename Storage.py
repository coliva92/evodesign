import sys
import json
from .Population import Individual
from typing import List, Tuple
from .Algorithm import Algorithm
from .Prediction import Testing as Predictor_Testing
from .Prediction import ESMFold as Predictor_ESMFold
from .Fitness import Testing as Fitness_Testing
from .Fitness import NegativeRmsd as Fitness_NegativeRmsd
from .Fitness import NegativeGdt as Fitness_NegativeGdt
from .GA import GeneticAlgorithm as GA
from .GA.Selection import Overselection as GA_Selection_Overselection
from .GA.Selection import Tournament as GA_Selection_Tournament
from .GA.Selection import UniformRandom as GA_Selection_UniformRandom
from .GA.Selection import BestAgainstPopulation as GA_Selection_BestAgainstPopulation
from .GA.Recombination import CenterPointCrossover as GA_Recombination_CenterPointCrossover
from .GA.Recombination import SinglePointCrossover as GA_Recombination_SinglePointCrossover
from .GA.Recombination import TwoPointCrossover as GA_Recombination_TwoPointCrossover
from .GA.Recombination import UniformCrossover as GA_Recombination_UniformCrossover
from .GA.Mutation import SingleSwitch as GA_Mutation_SingleSwitch
from .GA.Mutation import MultipleSwitches as GA_Mutation_MultipleSwitches
from .GA.Replacement import WorseOut as GA_Replacement_WorseOut
from .GA.Replacement import TotalReplacement as GA_Replacement_TotalReplacement





def load_memento_from_json_file(filename: str) -> dict:
  with open(filename, 'rt', encoding='utf-8') as file:
    memento = json.load(file)
  return memento



def load_algorithm_from_memento(memento: dict) -> Algorithm:
  algorithm_class = getattr(sys.modules[__name__], memento['algorithmType'])
  algorithm_params, operation_params = {}, {}
  other_algorithm_params = {}
  for key, value in memento['algorithmParams'].items():
    if key.endswith('Params'):
      name = key.replace('Params', '')
      operation_params[name] = value
      continue
    if type(value) == str and value.startswith(memento['algorithmType']):
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
    algorithm.append_memento_and_cache(memento['__savedPopulations'])
  else:
    algorithm.append_memento_and_cache()
  return algorithm



def load_population_from_memento(memento: dict) -> Tuple[int, List[Individual]]:
  if '__savedPopulations' not in memento or \
      len(memento['__savedPopulations']) == 0:
    return 0, []
  with open(memento['__savedPopulations'][-1], 'rt', encoding='utf-8') as file:
    population_memento = json.load(file)
  iterationId = len(memento['__savedPopulations']) - 1
  population = [ Individual(**params) for params in population_memento ]
  return iterationId, population
