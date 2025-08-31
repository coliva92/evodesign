from .Selection import Selection
from pymoo.operators.selection.tournament import TournamentSelection
import numpy as np
import numpy.typing as npt


class Tournament(Selection):

    def __init__(self, tournament_size: int = 2) -> None:
        super().__init__(TournamentSelection(self.tournament, tournament_size))
        self.tournament_size = tournament_size

    @classmethod
    def tournament(
        cls, pop, P: npt.NDArray[np.int64], **kwargs
    ) -> npt.NDArray[np.int64]:
        selection_size, tournament_size = P.shape
        selected_parents = np.full(selection_size, -1, dtype=np.int64)
        for i in range(selection_size):
            fitness = np.array([pop[j].F[0] for j in P[i]])
            selected_parents[i] = P[i][fitness.argmin()]
        return selected_parents
