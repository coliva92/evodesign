from pymoo.core.problem import Problem
from ...Fitness.FitnessFunction import FitnessFunction
from ...Prediction.Predictor import Predictor
from ...Utils.Chain import Chain
from ...Utils.ChainFactory import ChainFactory
from ...Prediction.DirectoryManager import DirectoryManager
import numpy as np
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

    def _evaluate(
        self,
        x,
        out,
        *args,
        **kwargs,
    ) -> None:
        # using a list comprehension in this case is faster than vectorizing
        sequences = [ChainFactory.sequence_numpy_to_str(seq_numpy) for seq_numpy in x]
        # Note: x.shape = population_size x sequence_length
        if self.predictor.is_activated():
            self.predictor.do(sequences, self.predictor_directory)
            base_dir = self.predictor_directory.prediction_pdbs_dir
            prefix = self.predictor_directory.prefix
            model_chains = [
                ChainFactory.create_from_pdb(
                    os.path.join(base_dir, f"{prefix}_{i}.pdb")
                )
                for i in range(x.shape[0])
            ]
        else:
            model_chains = [
                ChainFactory.create_from_sequence(
                    ChainFactory.sequence_numpy_to_str(sequence_numpy)
                )
                for sequence_numpy in x
            ]
        tmp_terms = np.array(
            [
                self.fitness_fn.do(model_chain, self.ref_chain)
                for model_chain in model_chains
            ]
        )
        self.term_values = tmp_terms[:, 1:]
        out["F"] = -1.0 * tmp_terms[:, 0]
