from abc import ABC
from pymoo.core.problem import Problem
from ..Prediction.Predictor import Predictor
from ..Utils.Chain import Chain
from ..Prediction.DirectoryManager import DirectoryManager
import numpy as np


class CPD(Problem, ABC):

    def __init__(
        self,
        predictor: Predictor,
        ref_chain: Chain,
        predictor_directory: DirectoryManager,
    ):
        super().__init__(
            n_var=len(ref_chain.sequence),
            n_obj=1,  # mono objective
            n_eq_constr=0,
            n_ieq_constr=0,
            xl=0,  # working on a 20-letters amino acid alphabet
            xu=19,  # represented by integers from 0 to 19
            vtype=np.int64,
        )
        self.predictor = predictor
        self.ref_chain = ref_chain
        self.term_values = None
        self.predictor_directory = predictor_directory
