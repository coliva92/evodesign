from .GASimple import GASimple
from ..Fitness.FitnessFunction import FitnessFunction
from ..Metrics.Metric import Metric
from ..Prediction.Predictor import Predictor
import pandas as pd
from typing import List, Optional





class RandomSearch(GASimple):

    def __init__(self, 
                 fitness_fn: FitnessFunction,
                 max_generations: int, 
                 population_size: int, 
                 predictor: Predictor,
                 metrics: List[Metric], 
                 sort_columns: Optional[List[str]] = None,
                 sort_ascending: Optional[List[bool]] = None
                 ) -> None:
        super().__init__(fitness_fn,
                         max_generations, 
                         population_size, 
                         predictor, 
                         None, 
                         None, 
                         None, 
                         metrics, 
                         sort_columns, 
                         sort_ascending)



    def next_population(self, population: pd.DataFrame) -> pd.DataFrame:
        """
        Apply the selection, recombination and mutation operators to the given 
        population to create a new population for the next generation of the
        evolutionary algorithm.

        Parameters
        ----------
        population : pandas.DataFrame
            The population from which the new population will be created.

        Returns
        -------
        pandas.DataFrame
            The generated population.
        """
        return self.initial_population()
