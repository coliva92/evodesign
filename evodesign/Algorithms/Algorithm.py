from abc import ABC, abstractmethod
from ..RetrievableSettings import RetrievableSettings
from ..Utils.Chain import Chain
from ..Utils.SavingManager import SavingManager
from ..Prediction.Predictor import Predictor
from pymoo.core.algorithm import Algorithm as PyMOOAlgorithm
from pymoo.core.problem import Problem as PyMOOProblem
from pymoo.optimize import minimize
from pymoo.operators.sampling.rnd import IntegerRandomSampling
from typing import Optional


class Algorithm(RetrievableSettings, ABC):

    def __init__(
        self, max_generations: int, population_size: int, predictor: Predictor
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
        prediction_pdb_path: str,
    ) -> PyMOOProblem:
        raise NotImplementedError
    
    @abstractmethod
    def num_terms(self) -> int:
        raise NotImplementedError

    def run(self, ref_chain: Chain, saving: Optional[SavingManager] = None, **kwargs):
        if self._algorithm is None:
            self._algorithm = self._create_algorithm()
        prediction_pdb_path = "prediction.pdb.tmp"
        if saving is not None:
            prediction_pdb_path = saving.working_folder.prediction_pdb_path
        self._problem = self._create_problem(ref_chain, prediction_pdb_path)
        results = minimize(
            self._problem,
            self._algorithm,
            callback=saving,
            verbose=True,
            copy_algorithm=False,
            **kwargs
        )
        return results
