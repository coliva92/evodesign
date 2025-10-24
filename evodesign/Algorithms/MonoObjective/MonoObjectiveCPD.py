from pymoo.core.problem import Problem
from ...Fitness.FitnessFunction import FitnessFunction
from ...Prediction.Predictor import Predictor
from ...Prediction.Null import Null
from ...Utils.Chain import Chain, ChainFactory, numpy_sequence_to_str
from ...Prediction.DirectoryManager import DirectoryManager
import numpy as np
import numpy.typing as npt
import os


class MonoObjectiveCPD(Problem):

    def __init__(
        self,
        predictor: Predictor,
        fitness_fn: FitnessFunction,
        ref_chain: Chain,
        predictor_directory: DirectoryManager,
    ):
        self.fitness_fn = fitness_fn
        self.predictor = predictor
        self.ref_chain = ref_chain
        self.term_values = None
        self.predictor_directory = predictor_directory
        super().__init__(
            n_var=len(self.ref_chain.sequence),
            n_obj=1,  # mono objective
            n_eq_constr=0,
            n_ieq_constr=0,
            xl=0,  # working on a 20-letters amino acid alphabet
            xu=19,  # represented by integers from 0 to 19
            vtype=np.int64,
        )

    def compute_fitness(
        self,
        sequence_idx: int,
    ) -> npt.NDArray[np.float64]:
        pdb_path = os.path.join(
            self.predictor_directory.prediction_pdbs_dir,
            f"{self.predictor_directory.prefix}_{sequence_idx}.pdb",
        )
        model_chain = ChainFactory.create(pdb_path)
        return self.fitness_fn.do(model_chain, self.ref_chain)

    def _evaluate(
        self,
        x,
        out,
        *args,
        **kwargs,
    ) -> None:
        # using a list comprehension in this case is faster than vectorizing
        sequences = [numpy_sequence_to_str(seq_numpy) for seq_numpy in x]
        self.predictor.do(sequences, self.predictor_directory)
        # Note: x.shape = population_size x sequence_length
        tmp_terms = np.array([self.compute_fitness(i) for i in range(x.shape[0])])
        self.term_values = tmp_terms[:, 1:]
        out["F"] = -1.0 * tmp_terms[:, 0]
