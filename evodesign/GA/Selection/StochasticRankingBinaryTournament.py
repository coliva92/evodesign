import numpy as np
import numpy.typing as npt
from pymoo.core.population import Population
from pymoo.operators.selection.tournament import TournamentSelection

from .Selection import Selection


class StochasticRankingBinaryTournament(Selection):

    def __init__(
        self,
        rank_by_fitness_prob: float = 0.45,
    ) -> None:
        super().__init__(TournamentSelection(self.compare_by_rank, pressure=2))
        self.rank_by_fitness_prob = rank_by_fitness_prob
        return

    def _stochastic_ranking(
        self, pop: Population, *args, **kwargs
    ) -> npt.NDArray[np.int64]:
        population_size = pop.shape[0]
        N_sweeps = population_size
        sorted_indices = np.arange(population_size, dtype=np.int64)

        for _ in range(N_sweeps):
            swap_done = False
            for j in range(population_size - 1):
                idx1 = sorted_indices[j]
                idx2 = sorted_indices[j + 1]
                coin_toss = np.random.uniform(0, 1)

                cv1 = pop[idx1].CV[0] if pop[idx1].CV is not None else 0.0
                cv2 = pop[idx2].CV[0] if pop[idx2].CV is not None else 0.0
                f1 = pop[idx1].F[0]
                f2 = pop[idx2].F[0]

                # If both solutions are feasible (G <= 0) OR the coin toss is heads,
                # compare based on objective function (minimization)
                comp_by_fitness = (cv1 <= 0 and cv2 <= 0) or (
                    coin_toss < self.rank_by_fitness_prob
                )

                if comp_by_fitness:
                    if f1 > f2:
                        sorted_indices[j], sorted_indices[j + 1] = (
                            sorted_indices[j + 1],
                            sorted_indices[j],
                        )
                        swap_done = True
                else:
                    # Otherwise compare based on constraint violation (lesser violation wins)
                    if cv1 > cv2:
                        sorted_indices[j], sorted_indices[j + 1] = (
                            sorted_indices[j + 1],
                            sorted_indices[j],
                        )
                        swap_done = True
            if not swap_done:
                break
        return sorted_indices

    def compare_by_rank(
        self,
        pop: Population,
        P: npt.NDArray[np.int64],
        **kwargs,
    ) -> npt.NDArray[np.int64]:
        selection_size, tournament_size = P.shape
        assert tournament_size == 2
        sorted_indices = self._stochastic_ranking(pop)
        ranks = np.zeros(len(pop), dtype=int)
        for rank, original_idx in enumerate(sorted_indices):
            ranks[original_idx] = rank
        pop.set("rank", ranks)
        selected_parents = np.empty(selection_size, dtype=np.int64)
        for i in range(selection_size):
            idx_a, idx_b = P[i, 0], P[i, 1]
            rank_a = pop[idx_a].get("rank")
            rank_b = pop[idx_b].get("rank")
            if rank_a < rank_b:
                selected_parents[i] = idx_a
            else:
                selected_parents[i] = idx_b
        return selected_parents
