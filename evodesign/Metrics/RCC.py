from typing import Dict, Optional, Sequence

import numpy as np
import numpy.typing as npt

from ..Chemistry.Geometry import compute_contacts_map
from ..Files.Graph import create_graph_file
from ..System.Subprocess import stream_command_line
from .ContextInterface import ContextInterface
from .ESM2 import ESM2
from .NonStructuralMetric import NonStructuralMetric

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


class RCC(NonStructuralMetric):

    def __init__(
        self,
        path_to_quick_cliques: str,
        esm_model: ESM2 = ESM2(),
        distance_threshold: float = 8.0,
        min_separation: int = 6,
    ) -> None:
        super().__init__()
        self.path_to_quick_cliques = path_to_quick_cliques
        self.esm_model = esm_model
        self.distance_threshold = distance_threshold
        self.min_separation = min_separation
        return

    def _rmse(
        self,
        model_rcc_vector: npt.NDArray[np.int64],
        ref_rcc_vector: npt.NDArray[np.int64],
    ) -> float:
        return np.sqrt(np.mean((model_rcc_vector - ref_rcc_vector) ** 2))

    def _cos_similarity(
        self,
        model_rcc_vector: npt.NDArray[np.int64],
        ref_rcc_vector: npt.NDArray[np.int64],
    ) -> float:
        n1 = np.linalg.norm(model_rcc_vector)
        n2 = np.linalg.norm(ref_rcc_vector)
        return np.dot(model_rcc_vector, ref_rcc_vector) / (n1 * n2)

    def _euclidean_distance(
        self,
        model_rcc_vector: npt.NDArray[np.int64],
        ref_rcc_vector: npt.NDArray[np.int64],
    ) -> float:
        return np.linalg.norm(model_rcc_vector - ref_rcc_vector)

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
        filename: Optional[str] = None,
    ) -> npt.NDArray[np.int64]:
        # generate a temporary filename if none is provided
        if filename is None:
            import datetime

            filename = datetime.datetime.now().strftime("%Y%m%d%H%M%S") + ".txt"

        # create the graph file that is the input to quick-cliques
        create_graph_file(contact_map, filename)

        # initialize the final RCC vector
        rcc = np.zeros(RCC_LENGTH, dtype=np.int64)

        # run quick-cliques and construct the RCC at the same time
        cmd = [
            self.path_to_quick_cliques,
            f"--input-file={filename}",
            f"--algorithm=tomita",
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

    def do(self, model_rcc_vector, reference_rcc_vector):
        rmse = self._rmse(model_rcc_vector, reference_rcc_vector)
        cos_sim = self._cos_similarity(model_rcc_vector, reference_rcc_vector)
        dist = self._euclidean_distance(model_rcc_vector, reference_rcc_vector)
        return rmse, cos_sim, dist

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        ref_rcc_vector = context.get_extra_param_value("ref_rcc_vector")
        if ref_rcc_vector is None:
            ref_contact_map = context.get_extra_param_value("ref_contact_map")
            if ref_contact_map is None:
                ref_atoms = context.get_reference_chain().ca_atoms
                ref_contact_map = compute_contacts_map(
                    ref_atoms,
                    self.distance_threshold,
                    self.min_separation,
                )
                context.set_extra_param_value("ref_contact_map", ref_contact_map)
            ref_rcc_vector = self.compute_rcc_vector(ref_contact_map)
            context.set_extra_param_value("ref_rcc_vector", ref_rcc_vector)
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
        model_rcc_vector = self.compute_rcc_vector(predicted_contacts)
        return self.do(model_rcc_vector, ref_rcc_vector)
