from abc import ABC
from ...RetrievableSettings import RetrievableSettings
from pymoo.core.mutation import Mutation as PyMOOMutation


class Mutation(RetrievableSettings, ABC):

    def __init__(
        self,
        sequence_mutation_prob: float,
        residue_mutation_prob: float,
        pymoo_mutation: PyMOOMutation,
    ) -> None:
        super().__init__()
        self.sequence_mutation_prob = sequence_mutation_prob
        self.residue_mutation_prob = residue_mutation_prob
        self._pymoo_mutation = pymoo_mutation

    # def do(
    #     self,
    #     rng: np.random.Generator,
    #     children: npt.NDArray[np.int64],
    #     constraints: Optional[npt.NDArray[np.float64]] = None,
    # ) -> npt.NDArray[np.int64]:
    #     raise NotImplementedError
