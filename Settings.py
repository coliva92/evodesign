import sys
import json
from .Algorithm import Algorithm
from .Prediction import *
from .Fitness import *
from .GA import *
from .GA.Selection import *
from .GA.Recombination import *
from .GA.Mutation import *
from .GA.ChildSelection import *
from .GA.Replacement import *





def _set_auto_params_GA(params: dict, 
                        nonClassParams: dict, 
                        classes: dict) -> dict:
  population_size = nonClassParams['populationSize']
  numOffspringsPerPair = 1 \
                         if isinstance(classes['recombination'], list) \
                         else 2
  if 'childrenSelection' not in params:
    params['childrenSelection'] = 'GA_ChildrenSelection_SingleBest'
    classes['childrenSelection'] = GA_ChildrenSelection_SingleBest
  params['childrenSelectionParams'] = {}
  params['selectionParams']['selectionSize'] = \
    params['childrenSelectionParams']['maxInputSize'] = \
    numOffspringsPerPair * population_size
  return params



_AUTO_PARAMS_SETTER = {
  'GA_Simple': _set_auto_params_GA
}



def load_algorithm_from_settings(filename: str, 
                                 additionalIterations: int = 0 
                                 ) -> Algorithm:  
  with open(filename, 'rt', encoding='utf-8') as json_file:
    settings = json.load(json_file)
  i = settings['algorithmType'].find('_')
  algorithm_name = settings['algorithmType'][:i]
  params = settings['algorithmParams']

  def is_class_key(key: str) -> bool:
    if key.endswith('Params') or type(params[key]) != str:
      return False
    return params[key].startswith('Predictor') or \
           params[key].startswith('Fitness') or \
           params[key].startswith(algorithm_name)
  
  def is_non_class_key(key: str) -> bool:
    if key.endswith('Params'):
      return False
    if type(params[key]) != str:
      return True
    return not params[key].startswith('Predictor') and \
           not params[key].startswith('Fitness') and \
           not params[key].startswith(algorithm_name)

  non_class_params = {
    key: params[key] for key in filter(is_non_class_key, params)
  }
  classes = {
    key: getattr(sys.modules[__name__], params[key])
    for key in filter(is_class_key, params)
  }
  set_auto_params = _AUTO_PARAMS_SETTER[settings['algorithmType']]
  params = set_auto_params(params, non_class_params, classes)
  class_params = {
    key: the_class(**params[f'{key}Params'])
    if f'{key}Params' in params else the_class()
    for key, the_class in classes.items()
  }
  
  algorithm_params = { **non_class_params, **class_params }
  algorithm_params['numIterations'] += additionalIterations
  algorithm_class = getattr(sys.modules[__name__], settings['algorithmType'])
  return algorithm_class(**algorithm_params)
