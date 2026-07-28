import os
from typing import List, Tuple

import numpy as np
import numpy.typing as npt

from ...Chemistry.Chain import Chain
from ...Chemistry.ChainFactory import ChainFactory
from ...CPD import CPD
from ...Fitness.FitnessFunction import FitnessFunction




class SingleObjectiveCPD(CPD):

    def __init__(self, fitness_fn: FitnessFunction, **kwargs):
        super().__init__(n_obj=1, **kwargs)
        self.fitness_fn = fitness_fn
        self.archive = {}
        return



    def _predict_structures(self, sequences: List[str]) -> List[Chain]:
        self.predictor.do(sequences, self.predictor_directory)
        base_dir = self.predictor_directory.prediction_pdbs_dir
        prefix = self.predictor_directory.prefix
        model_chains = [
            ChainFactory.create_from_pdb(os.path.join(base_dir, f"{prefix}_{i}.pdb"))
            for i in range(len(sequences))
        ]
        return model_chains



    def _compute_term_values(
        self,
        model_chains: List[Chain],
    ) -> npt.NDArray[np.float64]:
        return np.array(
            [
                self.fitness_fn.do(model_chain, self.ref_chain)
                for model_chain in model_chains
            ]
        )



    def _find_sequences_in_archive(
        self, sequences: List[str]
    ) -> Tuple[List[int], List[int]]:
        idx_archived, idx_not_archived = [], []
        for k, seq in enumerate(sequences):
            if seq in self.archive:
                idx_archived.append(k)
                continue
            idx_not_archived.append(k)
        return idx_archived, idx_not_archived



    def _update_archive(
        self,
        sequences: List[str],
        term_values: npt.NDArray[np.float64],
        indices: List[int],
    ) -> None:
        for k, i in enumerate(indices):
            self.archive[sequences[i]] = term_values[k]
        return



    def _evaluate(
        self,
        x,
        out,
        *args,
        **kwargs,
    ) -> None:
        # Note: x.shape = population_size x sequence_length
        sequences = [ChainFactory.sequence_numpy_to_str(seq_numpy) for seq_numpy in x]
        idx_archived, idx_not_archived = self._find_sequences_in_archive(sequences)
        if len(idx_not_archived) == 0:
            tmp_terms = np.array([self.archive[seq] for seq in sequences])
            self.term_values = tmp_terms[:, 1:]
            out["F"] = -1.0 * tmp_terms[:, 0]
            return
        if self.fitness_fn.requires_structure_predictor():
            model_chains = self._predict_structures(
                [sequences[k] for k in idx_not_archived]
            )
        else:
            model_chains = [
                ChainFactory.create_from_numpy(sequence_numpy)
                for sequence_numpy in x[idx_not_archived, :]
            ]
        new_terms = self._compute_term_values(model_chains)
        self._update_archive(sequences, new_terms, idx_not_archived)
        tmp_terms = np.zeros((len(sequences), new_terms.shape[1]))
        tmp_terms[idx_not_archived, :] = new_terms
        if len(idx_archived) > 0:
            tmp_terms[idx_archived, :] = np.array(
                [self.archive[sequences[k]] for k in idx_archived]
            )
        self.term_values = tmp_terms[:, 1:]  # from 1 to N - 1 are the term values
        out["F"] = -1.0 * tmp_terms[:, 0]  # element 0 is the fitness value
        return
