from abc import ABC, abstractmethod
from ..RetrievableSettings import RetrievableSettings
from ..Utils.Chain import Chain
from ..Callbacks.StorageManager import StorageManager
from ..Prediction.Predictor import Predictor
from ..Prediction.DirectoryManager import DirectoryManager
from pymoo.core.algorithm import Algorithm as PyMOOAlgorithm
from pymoo.core.problem import Problem as PyMOOProblem
from pymoo.core.callback import Callback
from pymoo.optimize import minimize
from pymoo.operators.sampling.rnd import IntegerRandomSampling
import numpy as np
import numpy.typing as npt
from typing import Optional


class Algorithm(RetrievableSettings, ABC):

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
        self._sampling = IntegerRandomSampling()
        self._algorithm = None
        self._problem = None

    @abstractmethod
    def _create_algorithm(self) -> PyMOOAlgorithm:
        raise NotImplementedError

    @abstractmethod
    def _create_problem(
        self,
        ref_chain: Chain,
        predictor_directory: DirectoryManager,
        aa_profile: Optional[npt.NDArray[np.float64]] = None,
    ) -> PyMOOProblem:
        raise NotImplementedError

    @abstractmethod
    def num_terms(self) -> int:
        raise NotImplementedError

    def _callbacks_chain(
        self,
        storage: StorageManager,
    ) -> Callback:
        return storage

    def run(
        self,
        ref_chain: Chain,
        storage: StorageManager,
        aa_profile: Optional[npt.NDArray[np.float64]] = None,
        **kwargs,
    ):
        if self._algorithm is None:
            self._algorithm = self._create_algorithm()
        self._problem = self._create_problem(ref_chain, storage.predictor_directory, aa_profile)
        results = minimize(
            self._problem,
            self._algorithm,
            callback=self._callbacks_chain(storage),
            verbose=True,
            copy_algorithm=False,
            **kwargs,
        )
        return results
