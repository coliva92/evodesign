from .Mutation import Mutation
from evodesign.Sequence import AMINO_ACID_INT_VALUES
import numpy as np
import numpy.typing as npt
from typing import Optional


class Swap(Mutation):

    def __init__(self, mutation_probability: float = 1.0, num_swaps: int = 1) -> None:
        super().__init__(mutation_probability)
        self.num_swaps = num_swaps

    def do(
        self,
        rng: np.random.Generator,
        children: npt.NDArray[np.int64],
        constraints: Optional[npt.NDArray[np.float64]] = None,
    ) -> npt.NDArray[np.int64]:
        rows_mask = rng.random(children.shape[0]) <= self.mutation_probability
        swap_indices = rng.integers(
            0, children.shape[1], (children.shape[0], self.num_swaps)
        ).T
        if constraints is None:
            mutations = rng.choice(AMINO_ACID_INT_VALUES[1:], size=children.shape)
            children[rows_mask, swap_indices[rows_mask]] = (
                children[rows_mask, swap_indices[rows_mask]]
                + mutations[rows_mask, swap_indices[rows_mask]]
            ) % AMINO_ACID_INT_VALUES.shape[0]
        else:
            # constraints.shape = sequence_length x len(AMINO_ACID_INT_VALUES)
            mutations = np.vstack(
                [
                    rng.choice(
                        AMINO_ACID_INT_VALUES, children.shape[0], p=constraints[i]
                    )
                    for i in range(constraints.shape[0])
                ]
            ).T
            children[rows_mask, swap_indices[rows_mask]] = mutations[
                rows_mask, swap_indices[rows_mask]
            ]
        return children
