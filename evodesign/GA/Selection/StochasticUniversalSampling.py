from .Selection import Selection
from ...Context import Context
import numpy as np
import numpy.typing as npt
import pandas as pd





class StochasticUniversalSampling(Selection):

    def _params(self) -> dict:
        params = super()._params()
        params["scaling_factor"] = self._scaling_factor
        return params



    def __init__(self,
                 scaling_factor: float = 1.2
                 ) -> None:
        super().__init__()
        self._scaling_factor = scaling_factor
        self._prob_distribution = None



    def select_parents(self, 
                       population: pd.DataFrame, 
                       context: Context
                       ) -> pd.DataFrame:
        n = len(population)
        if self._prob_distribution is None:
            self._prob_distribution = self._cumm_prob_distribution(n)
        population.sort_values(by=context.sort_columns,
                               ascending=context.sort_ascending,
                               ignore_index=True,
                               inplace=True)
        selected_parents = pd.DataFrame(columns=population.columns)
        inc = 1 / n
        r = context.rng.uniform(0.0, inc)
        selections_count = 0
        i = n - 1
        while selections_count < n:
            while r <= self._prob_distribution[i]:
                row = population.iloc[i]
                selected_parents = pd.concat([ selected_parents, row.to_frame().T ],
                                             ignore_index=True)
                selections_count += 1
                r += inc
            i -= 1
        return selected_parents
    


    def _cumm_prob_distribution(self, 
                                population_size: int
                                ) -> npt.NDArray[np.float64]:
        # assuming population is sorted in descending order by fitness value
        probabilities = [
            self._linear_ranking_prob(i, population_size)
            for i in range(population_size)
        ]
        for i in reversed(range(population_size - 1)):
            probabilities[i] += probabilities[i + 1]
        return np.array(probabilities)



    def _linear_ranking_prob(self, 
                             idx: int,
                             population_size: int
                             ) -> float:
        s = self._scaling_factor
        i = population_size - idx # rank = # worse solutions after individual i
        n = population_size
        return (2 - s) / n + (2 * i * (s - 1)) / (n * (n - 1))
    


    def _rank(self, 
              idx: int, 
              population_size: int) -> int:
        # TODO ¿qué sucede con el rank cuando hay empates en la población?
        # TODO por ahora el rank se asigna suponiendo que no existen empates
        return population_size - idx
