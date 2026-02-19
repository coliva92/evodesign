from .NonStructuralMetric import NonStructuralMetric
from .ContextInterface import ContextInterface
from .ESM2 import ESM2
from .Normalization.Formulas import cos_similarity
import numpy as np
import numpy.typing as npt
from typing import Optional, Dict, Tuple, List


class ESM2ContactMapPredictedRef(NonStructuralMetric):

    def __init__(
        self,
        esm_model: ESM2 = ESM2(),
        submap_indices: Optional[List[int]] = None,
    ) -> None:
        super().__init__()
        self.esm_model = esm_model
        self.submap_indices = submap_indices
        return

    def do(
        self,
        predicted_contacts: npt.NDArray[np.float64],
        ref_contact_map: npt.NDArray[np.float64],
        **kwargs,
    ) -> Tuple[float, Optional[float]]:
        return cos_similarity(ref_contact_map.flatten(), predicted_contacts.flatten())

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        ref_contact_map = context.get_extra_param_value("esm2_ref_contact_map")
        scaling_factor = context.get_extra_param_value(
            "esm2_contact_map_scaling_factor"
        )
        if ref_contact_map is None or scaling_factor is None:
            ref_sequence = context.get_reference_chain().sequence
            _, ref_contact_map = self.esm_model.query_model(
                ref_sequence, submap_indices=self.submap_indices
            )
            a, b = np.min(ref_contact_map), np.max(ref_contact_map)
            scaling_factor = 1 / (b - a)
            context.set_extra_param_value(
                "esm2_ref_contact_map", ref_contact_map / scaling_factor
            )
            context.set_extra_param_value(
                "esm2_contact_map_scaling_factor", scaling_factor
            )
        predicted_contacts = context.get_extra_param_value("esm2_predicted_contacts")
        if predicted_contacts is None:
            model_sequence = context.get_model_chain().sequence
            model_desc_matrix, predicted_contacts = self.esm_model.query_model(
                model_sequence, submap_indices=self.submap_indices
            )
            context.set_extra_param_value("esm2_model_desc_matrix", model_desc_matrix)
            context.set_extra_param_value(
                "esm2_predicted_contacts", predicted_contacts / scaling_factor
            )
        sim = self.do(predicted_contacts, ref_contact_map)
        return {"cos_similarity": sim}
