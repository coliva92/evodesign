import numpy as np  # noqa: N999
import numpy.typing as npt

from ..ContextInterface import ContextInterface
from ..metrics import euclidean_distance, jensen_shannon_distance, rmse
from ..Normalization.Normalization import Normalization
from ..SequenceBasedEvaluator import SequenceBasedEvaluator
from .ModelInterface import ModelInterface


class ContactMapJensenShannon(SequenceBasedEvaluator):
    def __init__(
        self,
        model: ModelInterface | None = None,
        submap_indices: list[int] | None = None,
        rmse_normalization: Normalization | None = None,
        distance_normalization: Normalization | None = None,
    ) -> None:
        super().__init__()
        self.model = model
        self.submap_indices = submap_indices
        self.rmse_normalization = rmse_normalization
        self.distance_normalization = distance_normalization

    def _mean_jensen_shannon(
        self,
        predicted_map: npt.NDArray[np.float64],
        ref_contact_map: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        return np.mean(
            [
                jensen_shannon_distance(predicted_map[i, :], ref_contact_map[i, :])
                for i in range(predicted_map.shape[0])
            ]
        )

    def do(
        self,
        predicted_map: npt.NDArray[np.float64],
        ref_contact_map: npt.NDArray[np.float64],
        **kwargs,
    ) -> tuple[float, ...]:
        assert predicted_map.shape[0] == ref_contact_map.shape[0]
        row_idx, col_idx = np.triu_indices_from(predicted_map)
        u = predicted_map[row_idx, col_idx]
        row_idx, col_idx = np.triu_indices_from(ref_contact_map)
        v = ref_contact_map[row_idx, col_idx]
        rmse_val = rmse(u, v)
        norm_rmse = rmse_val
        if self.rmse_normalization is not None:
            norm_rmse = self.rmse_normalization.do(rmse_val)
        dist = euclidean_distance(u, v)
        norm_dist = dist
        if self.distance_normalization is not None:
            norm_dist = self.distance_normalization.do(dist)
        predicted_map += 0.0001
        predicted_map /= np.sum(predicted_map, axis=1)[:, np.newaxis]
        ref_contact_map += 0.0001
        ref_contact_map /= np.sum(ref_contact_map, axis=1)[:, np.newaxis]
        js_dist = self._mean_jensen_shannon(predicted_map, ref_contact_map)
        neg_js_dist = 1 - js_dist
        return rmse_val, norm_rmse, dist, norm_dist, neg_js_dist

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        ref_contact_map = context.get_extra_param_value("esm2_ref_contact_map")
        if ref_contact_map is None:
            if self.model is None:
                from .ESM2LocalAllLayers import ESM2LocalAllLayers

                self.model = ESM2LocalAllLayers()
            ref_sequence = context.get_reference_chain().sequence
            _, ref_contact_map = self.model.query_model(
                ref_sequence,
                submap_indices=self.submap_indices,
            )
            context.set_extra_param_value("esm2_ref_contact_map", ref_contact_map)
        predicted_map = context.get_extra_param_value("esm2_predicted_contacts")
        if predicted_map is None:
            if self.model is None:
                from .ESM2LocalAllLayers import ESM2LocalAllLayers

                self.model = ESM2LocalAllLayers()
            model_sequence = context.get_model_chain().sequence
            model_desc_matrix, predicted_map = self.model.query_model(
                model_sequence, submap_indices=self.submap_indices
            )
            context.set_extra_param_value("esm2_model_desc_matrix", model_desc_matrix)
            context.set_extra_param_value("esm2_predicted_contacts", predicted_map)
        rmse_val, norm_rmse, dist, norm_dist, js_dist = self.do(
            predicted_map,
            ref_contact_map,
        )
        return {
            "rmse": rmse_val,
            "norm_rmse": norm_rmse,
            "distance": dist,
            "norm_distance": norm_dist,
            "js_distance": js_dist,
        }
