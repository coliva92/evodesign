import sys
from .Algorithms.GA2 import GA2 as Algorithms_GA2
from .Algorithms.NSGA2 import NSGA2 as Algorithms_NSGA2
from .Algorithms.GA1 import GA1 as Algorithms_GA1
from .Algorithms.GA3 import GA3 as Algorithms_GA3
from .Fitness.Rmsd import Rmsd as Fitness_Rmsd
from .Fitness.Gdt import Gdt as Fitness_Gdt
from .Fitness.Plddt import Plddt as Fitness_Plddt
from .Fitness.Rastrigin import Rastrigin as Fitness_Rastrigin
from .Fitness.Experimental.RmsdCyclization import RmsdCyclization \
  as Fitness_Experimental_RmsdCyclization
from .Fitness.Experimental.RmsdPlddt import RmsdPlddt \
  as Fitness_Experimental_RmsdPlddt
from .Fitness.Experimental.RmsdCyclizationPlddt import RmsdCyclizationPlddt \
  as Fitness_Experimental_RmsdCyclizationPlddt
from .Fitness.Experimental.GdtCyclization import GdtCyclization \
  as Fitness_Experimental_GdtCyclization
from .Fitness.Experimental.GdtPlddt import GdtPlddt \
  as Fitness_Experimental_GdtPlddt
from .Fitness.Experimental.GdtCyclizationPlddt import GdtCyclizationPlddt \
  as Fitness_Experimental_GdtCyclizationPlddt
from .Fitness.Experimental.SideChainPacking import SideChainPacking \
  as Fitness_Experimental_SideChainPacking
from .Fitness.Experimental.CyclizationPlddt import CyclizationPlddt \
  as Fitness_Experimental_CyclizationPlddt
from .Fitness.Experimental.GdtContactMapCyclizationPlddt import GdtContactMapCyclizationPlddt \
  as Fitness_Experimental_GdtContactMapCyclizationPlddt
from .Fitness.Experimental.DescGdtContactMapCyclizationPlddt import DescGdtContactMapCyclizationPlddt \
  as Fitness_Experimental_DescriptorsGdtContactMapCyclizationPlddt
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
