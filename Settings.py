import sys
import json
from .Algorithm import Algorithm
from .Prediction import *
from .Fitness import *
from .GA import *
from .GA.Selection import *
from .GA.Recombination import *
from .GA.Mutation import *
from .GA.Replacement import *





def load_algorithm_from_settings(filename: str) -> Algorithm:  
  with open(filename, 'rt', encoding='utf-8') as json_file:
    settings = json.load(json_file)
  i = settings['algorithmType'].find('_')
  algorithm_name = settings['algorithmType'][:i]

  def is_class_key(key: str) -> bool:
    if key.endswith('Params') or type(settings[key]) != str:
      return False
    return settings[key].startswith('Predictor') or \
           settings[key].startswith('Fitness') or \
           settings[key].startswith(algorithm_name)
  
  def is_non_class_key(key: str) -> bool:
    if key.endswith('Params'):
      return False
    if type(settings[key]) != str:
      return True
    return not settings[key].startswith('Predictor') and \
           not settings[key].startswith('Fitness') and \
           not settings[key].startswith(algorithm_name)

  non_class_params = {
    key: settings[key] for key in filter(is_non_class_key, 
                                         settings['algorithmParams'])
  }
  classes = {
    key: getattr(sys.modules[__name__], settings[key]) \
      for key in filter(is_class_key, settings['algorithmParams'])
  }
  class_params = {
    key: the_class(**settings[f'{key}Params']) \
         if f'{key}Params' in settings else the_class() \
         for key, the_class in classes.items()
  }

  if '__savedPopulations' in settings:
    non_class_params['populationFilenames'] = settings['__savedPopulations']
  
  algorithm_class = getattr(sys.modules[__name__], settings['algorithmType'])
  return algorithm_class(**{ **non_class_params, **class_params })
