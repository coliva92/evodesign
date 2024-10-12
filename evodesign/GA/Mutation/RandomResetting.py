from .Mutation import Mutation
from evodesign.Sequence import AMINO_ACID_INT_VALUES
import numpy as np
import numpy.typing as npt
from typing import Optional


class RandomResetting(Mutation):

    def __init__(
        self,
        mutation_probability: float = 1.0,
        residue_mutation_probability: float = 0.1,
    ) -> None:
        super().__init__(mutation_probability)
        self.residue_mutation_probability = residue_mutation_probability

    def do(
        self,
        rng: np.random.Generator,
        children: npt.NDArray[np.int64],
        constraints: Optional[npt.NDArray[np.float64]] = None,
    ) -> npt.NDArray[np.int64]:
        rows_mask = rng.random(children.shape[0]) <= self.mutation_probability
        mask = rng.random(children.shape) <= self.residue_mutation_probability
        mask[~rows_mask] = False
        if constraints is None:
            mutations = rng.choice(AMINO_ACID_INT_VALUES[1:], size=children.shape)
            children[mask] = (
                children[mask] + mutations[mask]
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
            children[mask] = mutations[mask]
        return children
