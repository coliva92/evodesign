from abc import ABC, abstractmethod
from ..RetrievableSettings import RetrievableSettings
from ..Chemistry.Chain import Chain
from ..Callbacks.StorageManager import StorageManager
from ..Prediction.Predictor import Predictor
from ..Prediction.DirectoryManager import DirectoryManager
from pymoo.core.algorithm import Algorithm as PyMOOAlgorithm
from pymoo.core.problem import Problem as PyMOOProblem
from ..Callbacks.CallbackCollection import CallbackCollection
from pymoo.optimize import minimize
from pymoo.core.result import Result
from ..Problems.IntegerProfileSampling import IntegerProfileSampling
import numpy as np
import numpy.typing as npt
from typing import Optional


class AlgorithmFactory(RetrievableSettings, ABC):

    def __init__(
        self,
        max_generations: int,
        population_size: int,
        predictor: Predictor,
    ):
        super().__init__()
        self.max_generations = max_generations
        self.population_size = population_size
        self.predictor = predictor
        self._sampling = IntegerProfileSampling()
        return

    @abstractmethod
    def create_algorithm(self) -> PyMOOAlgorithm:
        raise NotImplementedError

    @abstractmethod
    def create_problem(
        self,
        ref_chain: Chain,
        predictor_directory: DirectoryManager,
        aa_profile: Optional[npt.NDArray[np.float64]] = None,
    ) -> PyMOOProblem:
        raise NotImplementedError

    def create_callbacks(
        self,
        storage: StorageManager,
    ) -> CallbackCollection:
        return CallbackCollection([ storage ])
