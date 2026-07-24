from .Generational import Generational
from ....Problems.MonoObjective.PenalizedMonoCPD import PenalizedMonoCPD
from ....Chemistry.Chain import Chain
from ....Prediction.DirectoryManager import DirectoryManager
from pymoo.algorithms.soo.nonconvex.ga import GA
from ....GA.Termination.MaximumDiversityLoss import MaximumDiversityLoss
from ....GA.Replacement.PyMOO.Generational import (
    Generational as GenerationalReplacement,
)
from ....GA.Selection.StochasticRankingBinaryTournament import StochasticRankingBinaryTournament
from typing import Optional
import numpy as np
import numpy.typing as npt


class GenerationalPenalized(Generational):

    def __init__(
        self,
        identity_tol: float = 0.3,
        **kwargs
    ):
        super().__init__(**kwargs)
        assert isinstance(kwargs["selection"], StochasticRankingBinaryTournament)
        self.identity_tol = identity_tol
        return
    
    def create_problem(
        self,
        ref_chain: Chain,
        predictor_directory: DirectoryManager,
        aa_profile: Optional[npt.NDArray[np.float64]] = None,
    ) -> PenalizedMonoCPD:
        return PenalizedMonoCPD(
            identity_tol=self.identity_tol,
            ref_chain=ref_chain,
            fitness_fn=self.fitness_fn,
            predictor=self.predictor,
            predictor_directory=predictor_directory,
            aa_profile=aa_profile
        )
    
    def create_algorithm(self) -> GA:
        replacement = GenerationalReplacement()
        termination = MaximumDiversityLoss(
            self.max_generations, 
            self.diversity_loss_tol, 
            self.sample_size
        )
        algorithm = GA(
            pop_size=self.population_size,
            n_offsprings=self.population_size,
            sampling=self._sampling,
            selection=self.selection._pymoo_selection,
            crossover=self.crossover._pymoo_crossover,
            mutation=self.mutation._pymoo_mutation,
            survival=replacement,
            eliminate_duplicates=False,
        )
        algorithm.termination = termination
        return algorithm
