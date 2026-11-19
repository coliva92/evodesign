from abc import ABC  # noqa: N999

import numpy as np
import numpy.typing as npt
from pymoo.core.problem import Problem

from .Prediction.DirectoryManager import DirectoryManager
from .Prediction.Predictor import Predictor
from .Protein.Chain import Chain
from .Protein.Sequences import AMINO_ACIDS


class CPD(Problem, ABC):  # Computational Protein Design
    def __init__(
        self,
        ref_chain: Chain,
        predictor: Predictor,
        predictor_directory: DirectoryManager,
        aa_profile: npt.NDArray[np.float64] | None = None,
        **kwargs,
    ):
        if "n_ieq_constr" not in kwargs:
            kwargs["n_ieq_constr"] = 0
        super().__init__(
            n_var=len(ref_chain.sequence),
            n_eq_constr=0,
            xl=0,
            xu=len(AMINO_ACIDS) - 1,
            vtype=np.int64,
            **kwargs,
        )
        self.ref_chain = ref_chain
        self.predictor = predictor
        self.predictor_directory = predictor_directory
        self.term_values = None
        self.aa_profile = aa_profile
