from ..CPD import CPD
from ...Fitness.FitnessFunction import FitnessFunction
from ...Prediction.Predictor import Predictor
from ...Utils.Chain import Chain
from ...Utils.ChainFactory import ChainFactory
from ...Prediction.DirectoryManager import DirectoryManager
import numpy as np
import os


class MonoCPD(CPD):

    def __init__(
        self,
        fitness_fn: FitnessFunction,
        predictor: Predictor,
        ref_chain: Chain,
        predictor_directory: DirectoryManager,
    ):
        super().__init__(predictor, ref_chain, predictor_directory)
        self.fitness_fn = fitness_fn

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
                ChainFactory.create_from_numpy(sequence_numpy)
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
