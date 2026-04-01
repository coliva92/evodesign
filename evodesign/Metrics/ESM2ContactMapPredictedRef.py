from .NonStructuralMetric import NonStructuralMetric
from .ContextInterface import ContextInterface
from .ESM2 import ESM2
from .Normalization.Formulas import cos_similarity
import numpy as np
import numpy.typing as npt
from typing import Optional, Dict, Tuple, List
from .Normalization.Normalization import Normalization
from scipy.stats import entropy
from scipy.spatial.distance import jensenshannon


class ESM2ContactMapPredictedRef(NonStructuralMetric):

    def __init__(
        self,
        esm_model: ESM2 = ESM2(),
        submap_indices: Optional[List[int]] = None,
        regularization: Optional[Normalization] = None,
    ) -> None:
        super().__init__()
        self.esm_model = esm_model
        self.submap_indices = submap_indices
        self.regularization = regularization
        return

    def do(
        self,
        predicted_contacts: npt.NDArray[np.float64],
        ref_contact_map: npt.NDArray[np.float64],
        **kwargs,
    ) -> Tuple[float, float]:
        assert(predicted_contacts.shape[0] == ref_contact_map.shape[0])
        mean_distance = np.mean([
            jensenshannon(predicted_contacts[i, :], ref_contact_map[i, :])
            for i in range(predicted_contacts.shape[0])
        ])
        return mean_distance

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        ref_contact_map = context.get_extra_param_value("esm2_ref_contact_map")
        if ref_contact_map is None:
            ref_sequence = context.get_reference_chain().sequence
            _, ref_contact_map = self.esm_model.query_model(
                ref_sequence, submap_indices=self.submap_indices
            )
            a, b = np.min(ref_contact_map), np.max(ref_contact_map)
            context.set_extra_param_value(
                "esm2_ref_contact_map", ref_contact_map
            )
        predicted_contacts = context.get_extra_param_value("esm2_predicted_contacts")
        if predicted_contacts is None:
            model_sequence = context.get_model_chain().sequence
            model_desc_matrix, predicted_contacts = self.esm_model.query_model(
                model_sequence, submap_indices=self.submap_indices
            )
            context.set_extra_param_value("esm2_model_desc_matrix", model_desc_matrix)
            context.set_extra_param_value(
                "esm2_predicted_contacts", predicted_contacts
            )
        distance = self.do(predicted_contacts, ref_contact_map)
        return {"distance": distance}
