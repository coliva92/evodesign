from abc import ABC, abstractmethod
from typing import Optional

import numpy as np
import numpy.typing as npt
from pymoo.core.algorithm import Algorithm as PyMOOAlgorithm
from pymoo.core.problem import Problem as PyMOOProblem

from .Callbacks.CallbackCollection import CallbackCollection
from .Callbacks.StorageManager import StorageManager
from .Chemistry.Chain import Chain
from .IntegerProfileSampling import IntegerProfileSampling
from .Prediction.DirectoryManager import DirectoryManager
from .Prediction.Predictor import Predictor
from .RetrievableSettings import RetrievableSettings




class AlgorithmFactory(RetrievableSettings, ABC):

    def __init__(
        self,
        max_generations: int,
        population_size: int,
        predictor: Predictor,
    ) -> None:
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
        return CallbackCollection([storage])
