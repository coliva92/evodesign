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
        rmse_norm_const: float = 1.0,
        dist_norm_const: float = 1.0,
    ) -> None:
        super().__init__()
        self.esm_model = esm_model
        self.submap_indices = submap_indices
        self.rmse_norm_const = rmse_norm_const
        self.dist_norm_const = dist_norm_const
        return
    
    def _rmse(
        self,
        predicted_contacts: npt.NDArray[np.float64],
        ref_contact_map: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        return np.sqrt(np.mean((predicted_contacts - ref_contact_map)**2))
    
    def _euclidean_distance(
        self,
        predicted_contacts: npt.NDArray[np.float64],
        ref_contact_map: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        return np.linalg.norm(predicted_contacts - ref_contact_map)
    
    def _jensen_shannon(
        self,
        predicted_contacts: npt.NDArray[np.float64],
        ref_contact_map: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        return np.mean([
            jensenshannon(predicted_contacts[i, :], ref_contact_map[i, :])
            for i in range(predicted_contacts.shape[0])
        ])

    def do(
        self,
        predicted_contacts: npt.NDArray[np.float64],
        ref_contact_map: npt.NDArray[np.float64],
        **kwargs,
    ) -> Tuple[float]:
        assert(predicted_contacts.shape[0] == ref_contact_map.shape[0])
        row_idx, col_idx = np.triu_indices_from(predicted_contacts)
        u = predicted_contacts[row_idx, col_idx]
        row_idx, col_idx = np.triu_indices_from(ref_contact_map)
        v = ref_contact_map[row_idx, col_idx]
        rmse = self._rmse(u, v)
        norm_rmse = 1 - (rmse / self.rmse_norm_const)
        dist = self._euclidean_distance(u, v)
        norm_dist = 1 - (dist / self.dist_norm_const)
        predicted_contacts += 0.0001
        predicted_contacts /= np.sum(predicted_contacts, axis=1)[:, np.newaxis]
        ref_contact_map += 0.0001
        ref_contact_map /= np.sum(ref_contact_map, axis=1)[:, np.newaxis]
        js_dist = self._jensen_shannon(predicted_contacts, ref_contact_map)
        return rmse, norm_rmse, dist, norm_dist, js_dist

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
        rmse, norm_rmse, dist, norm_dist, js_dist = self.do(predicted_contacts, ref_contact_map)
        return {
            "rmse": rmse,
            "norm_rmse": norm_rmse,
            "distance": dist,
            "norm_distance": norm_dist,
            "js_distance": js_dist
        }
