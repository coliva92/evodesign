import sys
from .Algorithms.GA2 import GA2 as Algorithms_GA2
from .Algorithms.NSGA2 import NSGA2 as Algorithms_NSGA2
from .Algorithms.GA1 import GA1 as Algorithms_GA1
from .Algorithms.GA3 import GA3 as Algorithms_GA3
from .Metrics.Plddt import Plddt as Metrics_Plddt
from .Metrics.Rmsd import Rmsd as Metrics_Rmsd
from .Metrics.Gdt import Gdt as Metrics_Gdt
from .Metrics.ContactMapRms import ContactMapRms as Metrics_ContactMapRms
from .Metrics.Cyclization import Cyclization as Metrics_Cyclization
from .Metrics.Normalization.Reciprocal import Reciprocal \
  as Metrics_Normalization_Reciprocal
from .Metrics.Normalization.Cyclization import Cyclization \
  as Metrics_Normalization_Cyclization
from .Fitness.WeightedMean import WeightedMean as Fitness_WeightedMean
from .Fitness.LinearCombination import LinearCombination \
  as Fitness_LinearCombination
from .Fitness.Experimental.Rastrigin import Rastrigin as Fitness_Rastrigin
from .GA.Mutation.RandomResetting import RandomResetting \
  as GA_Mutation_RandomResetting
from .GA.Mutation.Switch import Switch as GA_Mutation_Switch
from .GA.Recombination.MiddlePointCrossover import MiddlePointCrossover \
  as GA_Recombination_MiddlePointCrossover
from .GA.Recombination.SinglePointCrossover import SinglePointCrossover \
  as GA_Recombination_SinglePointCrossover
from .GA.Recombination.TwoPointsCrossover import TwoPointsCrossover \
  as GA_Recombination_TwoPointsCrossover
from .GA.Recombination.UniformCrossover import UniformCrossover \
  as GA_Recombination_UniformCrossover
from .GA.Selection.Overselection import Overselection \
  as GA_Selection_Overselection
from .GA.Selection.Tournament import Tournament as GA_Selection_Tournament
from .GA.Selection.Uniform import Uniform as GA_Selection_Uniform
from .Prediction.AlphaFold import AlphaFold as Prediction_AlphaFold
from .Prediction.ESMFold import ESMFold as Prediction_ESMFold_Local
from .Prediction.ESMFoldColab import ESMFoldColab as Prediction_ESMFold_Colab
from .Prediction.ESMFoldRemoteApi import ESMFoldRemoteApi \
  as Prediction_ESMFold_RemoteApi
from .Prediction.Null import Null as Prediction_None
import copy





class Settings:

  @classmethod
  def parse(cls, settings: dict):
    class_name = list(settings.keys())[0]
    params = copy.deepcopy(settings[class_name])
    actual_class = getattr(sys.modules[__name__], class_name.replace('.', '_'))
    for key, item in params.items():
      if type(item) == dict:
        params[key] = cls.parse(item)
      if type(item) == list and type(item[0]) == dict:
        params[key] = [ cls.parse(s) for s in item ] \
                      if type(item[0]) == dict \
                      else item
    return actual_class(**params)
