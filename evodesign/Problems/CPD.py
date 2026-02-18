from abc import ABC
from pymoo.core.problem import Problem
from ..Prediction.Predictor import Predictor
from ..Chemistry.Chain import Chain
from ..Prediction.DirectoryManager import DirectoryManager
import numpy as np
import numpy.typing as npt
from typing import Optional


class CPD(Problem, ABC):

    def __init__(
        self,
        ref_chain: Chain,
        predictor: Predictor,
        predictor_directory: DirectoryManager,
        aa_profile: Optional[npt.NDArray[np.float64]] = None,
        num_objectives: int = 1, # mono objective
    ):
        super().__init__(
            n_var=len(ref_chain.sequence),
            n_obj=num_objectives,
            n_eq_constr=0,
            n_ieq_constr=0,
            xl=0,  # working on a 20-letters amino acid alphabet
            xu=19,  # represented by integers from 0 to 19
            vtype=np.int64,
        )
        self.ref_chain = ref_chain
        self.predictor = predictor
        self.predictor_directory = predictor_directory
        self.term_values = None
        self.aa_profile = aa_profile
        return
