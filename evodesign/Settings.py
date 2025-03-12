import sys
import copy
import json
from .Algorithms.GA2 import GA2 as Algorithms_GA2
from .Algorithms.NSGA2 import NSGA2 as Algorithms_NSGA2
from .Algorithms.GASteadyState import GASteadyState as Algorithms_GASteadyState
from .Algorithms.GASimpleElitism import GASimpleElitism as Algorithms_GASimpleElitism
from .Algorithms.GASimple import GASimple as Algorithms_GASimple
from .Algorithms.RandomSearch import RandomSearch as Algorithms_RandomSearch
from .Metrics.RMSD import RMSD as Metrics_RMSD
from .Metrics.GDT import GDT as Metrics_GDT
from .Metrics.DistanceMap import DistanceMap as Metrics_DistanceMap
from .Metrics.Cyclization import Cyclization as Metrics_Cyclization
from .Metrics.TMScore import TMScore as Metrics_TMScore
from .Metrics.RosettaEnergyFunction import (
    RosettaEnergyFunction as Metrics_PyRosettaRef2015,
)
from .Metrics.iLearnDescriptors import iLearnDescriptors as Metrics_iLearnDescriptors
from .Metrics.ESM2Descriptors import ESM2Descriptors as Metrics_ESM2Descriptors
from .Metrics.ESM2Descriptors import ESM2Descriptors as Metrics_ESM2Descriptors
from .Metrics.ESM2DescriptorsRemoteAPI import (
    ESM2DescriptorsRemoteAPI as Metrics_ESM2DescriptorsRemoteAPI
)
from .Metrics.Rastrigin import Rastrigin as Metrics_Rastrigin
from .Metrics.Experimental.Rastrigin2 import (
    Rastrigin2 as Metrics_Experimental_Rastrigin2,
)
from .Metrics.Yan24.Geometry import Geometry as Metrics_Yan24_Geometry
from .Metrics.Yan24.Energy import Energy as Metrics_Yan24_Energy
from .Metrics.Yan24.Descriptors import Descriptors as Metrics_Yan24_Descriptors
from .Fitness.WeightedMean import WeightedMean as Fitness_WeightedMean
from .GA.Mutation.RandomResetting import RandomResetting as GA_Mutation_RandomResetting
from .GA.Mutation.Swap import Swap as GA_Mutation_Swap
from .GA.Recombination.MiddlePointCrossover import (
    MiddlePointCrossover as GA_Recombination_MiddlePointCrossover,
)
from .GA.Recombination.SinglePointCrossover import (
    SinglePointCrossover as GA_Recombination_SinglePointCrossover,
)
from .GA.Recombination.TwoPointsCrossover import (
    TwoPointsCrossover as GA_Recombination_TwoPointsCrossover,
)
from .GA.Recombination.UniformCrossover import (
    UniformCrossover as GA_Recombination_UniformCrossover,
)
from .GA.Selection.Tournament import Tournament as GA_Selection_Tournament
from .GA.Selection.Uniform import Uniform as GA_Selection_Uniform
from .Prediction.AlphaFold import AlphaFold as Prediction_AlphaFold
from .Prediction.ESMFold import ESMFold as Prediction_ESMFold
from .Prediction.ESMFoldColab import ESMFoldColab as Prediction_ESMFoldColab
from .Prediction.ESMFoldRemoteApi import ESMFoldRemoteApi as Prediction_ESMFoldRemoteApi
from .Prediction.Null import Null as Prediction_Null


def load(settings_path: str):
    with open(settings_path, "rt", encoding="utf-8") as json_file:
        settings = json.load(json_file)
    return parse(settings)


def parse(settings: dict):
    class_name = list(settings.keys())[0]
    params = copy.deepcopy(settings[class_name])
    actual_class = getattr(sys.modules[__name__], class_name.replace(".", "_"))
    for key, item in params.items():
        if type(item) == dict:
            params[key] = parse(item)
        if type(item) == list and type(item[0]) == dict:
            params[key] = [parse(s) for s in item] if type(item[0]) == dict else item
    return actual_class(**params)
