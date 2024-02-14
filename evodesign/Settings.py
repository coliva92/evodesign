import sys
from .Algorithms.PGPD import PGPD as Algorithms_PGPD
from .Algorithms.NSGA2 import NSGA2 as Algorithms_NSGA2
from .Fitness.Rmsd import Rmsd as Fitness_Rmsd
from .Fitness.Gdt import Gdt as Fitness_Gdt
from .Fitness.Plddt import Plddt as Fitness_Plddt
from .Fitness.Rastrigin import Rastrigin as Fitness_Rastrigin
from .Fitness.Experimental.Cyclization import Cyclization \
  as Fitness_Experimental_Cyclization
from .Fitness.Experimental.SideChainPacking import SideChainPacking \
  as Fitness_Experimental_SideChainPacking
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
      if type(item) == list:
        params[key] = [ cls.parse(s) for s in item ]
    return actual_class(**params)
