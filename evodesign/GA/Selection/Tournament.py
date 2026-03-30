from .Selection import Selection
from pymoo.operators.selection.tournament import TournamentSelection
import numpy as np
import numpy.typing as npt


class Tournament(Selection):

    def __init__(
        self,
        tournament_size: int = 2,
        win_probability: float = 1.0,
    ) -> None:
        super().__init__(TournamentSelection(self.tournament, tournament_size))
        self.tournament_size = tournament_size
        self.win_probability = win_probability
        return

    def tournament(
        self,
        pop,
        P: npt.NDArray[np.int64],
        **kwargs,
    ) -> npt.NDArray[np.int64]:
        selection_size, tournament_size = P.shape
        selected_parents = np.full(selection_size, -1, dtype=np.int64)
        if self.win_probability >= 1.0:
            for i in range(selection_size):
                fitness = np.array([pop[j].F[0] for j in P[i]])
                selected_parents[i] = P[i][fitness.argmin()]
            return selected_parents
        assert(tournament_size == 2)
        for i in range(selection_size):
            idx_a, idx_b = P[i]
            fit_a = pop[idx_a].F[0]
            fit_b = pop[idx_b].F[0]
            if fit_a < fit_b:
                best_idx, worst_idx = idx_a, idx_b
            else:
                best_idx, worst_idx = idx_b, idx_a
            if np.random.random() < self.win_probability:
                selected_parents[i] = best_idx
            else:
                selected_parents[i] = worst_idx
        return selected_parents
