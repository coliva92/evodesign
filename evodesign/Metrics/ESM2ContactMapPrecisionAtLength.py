from typing import Dict

import numpy as np
import numpy.typing as npt

from ..Chemistry.Geometry import compute_contacts_map
from .ContextInterface import ContextInterface
from .ESM2 import ESM2
from .NonStructuralMetric import NonStructuralMetric


class ESM2ContactMapPrecisionAtLength(NonStructuralMetric):

    def __init__(
        self,
        esm_model: ESM2 = ESM2(),
        distance_threshold: float = 8.0,
        min_separation: int = 6,
    ) -> None:
        super().__init__()
        self.esm_model = esm_model
        self.distance_threshold = distance_threshold
        self.min_separation = min_separation
        return

    def do(
        self,
        predicted_map: npt.NDArray[np.int64],
        reference_map: npt.NDArray[np.int64],
        num_contacts: int,
    ) -> float:
        assert predicted_map.shape[0] == predicted_map.shape[1]
        assert predicted_map.shape == reference_map.shape
        L = reference_map.shape[0]
        num_hits = 0
        for i in range(L):
            for j in range(i + 1, L):
                if predicted_map[i, j] == 1 and reference_map[i, j] == 1:
                    num_hits += 1
        return num_hits / num_contacts

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        ref_contact_map = context.get_extra_param_value("ref_contact_map")
        if ref_contact_map is None:
            ref_atoms = context.get_reference_chain().ca_atoms
            ref_contact_map = compute_contacts_map(
                ref_atoms,
                self.distance_threshold,
                self.min_separation,
            )
            context.set_extra_param_value("ref_contact_map", ref_contact_map)
        num_contacts = context.get_extra_param_value("num_ref_contacts")
        if num_contacts is None:
            num_contacts = np.triu(ref_contact_map, k=1).sum()
            context.set_extra_param_value("num_ref_contacts")
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
            num_contacts,
        )
        P_at_L = self.do(predicted_contacts, ref_contact_map, num_contacts)
        return P_at_L
