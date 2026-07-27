from .Generational import Generational
from ....Chemistry.Chain import Chain
from ....Prediction.DirectoryManager import DirectoryManager
from ...Problems.AlternatingCPD import AlternatingCPD
from ....Callbacks.StorageManager import StorageManager
from ....Callbacks.MonoFitnessFnAlternator import MonoFitnessFnAlternator
from ....Fitness.FitnessFunction import FitnessFunction
from ....Callbacks.CallbackCollection import CallbackCollection
from typing import Optional
import numpy as np
import numpy.typing as npt


class GenerationalAlternatingFitness(Generational):

    def __init__(
        self,
        alt_fitness_fn: FitnessFunction,
        alt_fitness_fn_every_nth_gen: int = 10,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.alt_fitness_fn = alt_fitness_fn
        self.alt_fitness_fn_every_nth_gen = alt_fitness_fn_every_nth_gen
        return

    def num_terms(self):
        return self.fitness_fn.num_terms() + self.alt_fitness_fn.num_terms()

    def create_problem(
        self,
        ref_chain: Chain,
        predictor_directory: DirectoryManager,
        aa_profile: Optional[npt.NDArray[np.float64]] = None,
    ) -> AlternatingCPD:
        return AlternatingCPD(
            ref_chain=ref_chain,
            fitness_fn=self.fitness_fn,
            alt_fitness_fn=self.alt_fitness_fn,
            predictor=self.predictor,
            predictor_directory=predictor_directory,
            aa_profile=aa_profile,
        )

    def create_callbacks(
        self,
        storage: StorageManager,
    ) -> CallbackCollection:
        callbacks = super().create_callbacks(storage)
        new_callback = MonoFitnessFnAlternator(
            self.alt_fitness_fn_every_nth_gen,
            storage,
        )
        callbacks.append(new_callback)
        return callbacks
