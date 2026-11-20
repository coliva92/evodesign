from __future__ import annotations  # noqa: N999

from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

from ...Files.Graph import create_graph_file
from ...Protein.Geometry import compute_contacts_map
from ...System.Subprocess import stream_command_line
from ..ContextInterface import ContextInterface
from ..ESM2.ModelInterface import ModelInterface as ESM2Interface
from ..metrics import euclidean_distance, norm_cosine_similarity, rmse
from ..Normalization.Normalization import Normalization
from ..SequenceBasedEvaluator import SequenceBasedEvaluator

RCC_LENGTH = 26
RCC_CLASSES = {
    # Length = 3
    (3,): 2,
    (1, 2): 1,
    (1, 1, 1): 0,
    # Length = 4
    (4,): 7,
    (1, 3): 6,
    (2, 2): 5,
    (1, 1, 2): 4,
    (1, 1, 1, 1): 3,
    # Length = 5
    (5,): 14,
    (2, 3): 13,
    (1, 2, 2): 12,
    (1, 4): 11,
    (1, 1, 3): 10,
    (1, 1, 1, 2): 9,
    (1, 1, 1, 1, 1): 8,
    # Length = 6
    (6,): 25,
    (3, 3): 24,
    (2, 4): 23,
    (2, 2, 2): 22,
    (1, 5): 21,
    (1, 2, 3): 20,
    (1, 1, 4): 19,
    (1, 1, 2, 2): 18,
    (1, 1, 1, 3): 17,
    (1, 1, 1, 1, 2): 16,
    (1, 1, 1, 1, 1, 1): 15,
}


class RCC(SequenceBasedEvaluator):  # Residue Cluster Class
    def __init__(
        self,
        path_to_quick_cliques: str,
        esm_model: ESM2Interface | None = None,  # other models could be supported later
        distance_threshold: float = 8.0,
        min_separation: int = 6,
        rmse_normalization: Normalization | None = None,
        distance_normalization: Normalization | None = None,
    ) -> None:
        super().__init__()
        self.path_to_quick_cliques = path_to_quick_cliques
        self.esm_model = esm_model
        self.distance_threshold = distance_threshold
        self.min_separation = min_separation
        self.rmse_normalization = rmse_normalization
        self.distance_normalization = distance_normalization

    @classmethod
    def get_rcc_class(cls, nodes: Sequence[int]) -> int:
        """
        Computes runs of contiguous residue indices along the sequence, sorts the run
        lengths in ascending order, and classifies cliques of size 3 to 6 into classes 0-25.
        Returns -1 if the clique size is outside 3-6 or unrecognized.
        """
        if not (3 <= len(nodes) <= 6):
            return -1

        # 1. Identify lengths of contiguous segments
        contiguous = []
        run_len = 1

        for i in range(1, len(nodes)):
            if nodes[i] == nodes[i - 1] + 1:
                run_len += 1
            else:
                contiguous.append(run_len)
                run_len = 1

        contiguous.append(run_len)
        contiguous.sort()

        # 2. Look up the partition directly
        return RCC_CLASSES.get(tuple(contiguous), -1)

    def compute_rcc_vector(
        self,
        contact_map: npt.NDArray[np.float64],
        filename: str,
    ) -> npt.NDArray[np.int64]:
        # create the graph file that is the input to quick-cliques
        create_graph_file(contact_map, filename)

        # initialize the final RCC vector
        rcc = np.zeros(RCC_LENGTH, dtype=np.int64)

        # run quick-cliques and construct the RCC at the same time
        cmd = [
            self.path_to_quick_cliques,
            f"--input-file={filename}",
            "--algorithm=tomita",
        ]
        for raw_line in stream_command_line(cmd):
            line = raw_line.strip()
            if not line or line.startswith(("NOTE:", "Reading", "tomita:")):
                continue

            # Process clique directly in memory
            nodes = sorted([int(token) for token in line.split()])
            rcc_class = self.get_rcc_class(nodes)
            if rcc_class > -1:
                rcc[rcc_class] += 1

        return rcc

    def do(
        self,
        model_rcc_vector: npt.NDArray[npt.float64],
        reference_rcc_vector: npt.NDArray[npt.float64],
        **kwargs,
    ) -> tuple[float, ...]:
        rmse_val = rmse(model_rcc_vector, reference_rcc_vector)
        norm_rmse = rmse_val
        if self.rmse_normalization is not None:
            norm_rmse = self.rmse_normalization.do(rmse_val)
        norm_cos = norm_cosine_similarity(model_rcc_vector, reference_rcc_vector)
        dist = euclidean_distance(model_rcc_vector, reference_rcc_vector)
        norm_dist = dist
        if self.distance_normalization is not None:
            norm_dist = self.distance_normalization.do(dist)
        return rmse_val, norm_rmse, norm_cos, dist, norm_dist

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:

        if self.esm_model is None:
            from ..ESM2.ESM2LocalAllLayers import ESM2LocalAllLayers

            # ESM-2 model is instanced only when `query_model` is called
            self.esm_model = ESM2LocalAllLayers()

        # generate a temporary filename if none is provided
        if self._tmp_filename is None:
            import datetime

            self._tmp_filename = (
                "Graph_" + datetime.datetime.now().strftime("%Y%m%d%H%M%S") + ".txt"
            )

        ref_rcc_vector = context.get_persistent_param_value("ref_rcc_vector")
        if ref_rcc_vector is None:
            ref_contact_map = context.get_persistent_param_value("ref_contact_map")
            if ref_contact_map is None:
                ref_atoms = context.get_reference_chain().ca_atoms
                ref_contact_map = compute_contacts_map(
                    ref_atoms,
                    self.min_separation,
                )
                context.set_persistent_param_value("ref_contact_map", ref_contact_map)
            ref_rcc_vector = self.compute_rcc_vector(
                ref_contact_map,
                self._tmp_filename,
            )
            context.set_persistent_param_value("ref_rcc_vector", ref_rcc_vector)
        predicted_contacts = context.get_extra_param_value("esm2_predicted_contacts")
        if predicted_contacts is None:
            model_sequence = context.get_model_chain().sequence
            model_desc_matrix, predicted_contacts = self.esm_model.query_model(
                model_sequence
            )
            context.set_extra_param_value("esm2_model_desc_matrix", model_desc_matrix)
            context.set_extra_param_value("esm2_predicted_contacts", predicted_contacts)
        predicted_contacts = self.esm_model.contacts_map_from_prediction(
            predicted_contacts,
            self.min_separation,
        )
        model_rcc_vector = self.compute_rcc_vector(
            predicted_contacts,
            self._tmp_filename,
        )
        rmse_val, norm_rmse, norm_cos, dist, norm_dist = self.do(
            model_rcc_vector,
            ref_rcc_vector,
        )
        return {
            "rmse": rmse_val,
            "norm_rmse": norm_rmse,
            "norm_cos_similarity": norm_cos,
            "distance": dist,
            "norm_distance": norm_dist,
        }
