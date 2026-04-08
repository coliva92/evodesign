from .NonStructuralMetric import NonStructuralMetric
from .ContextInterface import ContextInterface
from .ESM2 import ESM2
import numpy as np
import numpy.typing as npt
from typing import Optional, Dict, List, Tuple


class ESM2Descriptors(NonStructuralMetric):

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
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        u = model_desc_matrix.flatten()
        v = ref_desc_matrix.flatten()
        return np.sqrt(np.mean((u - v)**2))

    def _cos_similarity(
        self,
        u: npt.NDArray[np.float64],
        v: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        n1 = np.linalg.norm(u)
        n2 = np.linalg.norm(v)
        return np.dot(u, v) / (n1*n2)

    def _mean_cos_similarity(
        self,
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        values = []
        for i in range(model_desc_matrix.shape[0]):
            u = model_desc_matrix[i, :]
            v = ref_desc_matrix[i, :]
            values.append((self._cos_similarity(u, v) + 1) / 2)
        return np.mean(values)
    
    def _mean_euclidean_distance(
        self,
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        values = []
        for i in range(model_desc_matrix.shape[0]):
            u = model_desc_matrix[i, :]
            v = ref_desc_matrix[i, :]
            values.append(np.linalg.norm(u - v))
        return np.mean(values)

    def do(
        self,
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
        **kwargs,
    ) -> Tuple[float]:
        rmse = self._rmse(model_desc_matrix, ref_desc_matrix)
        norm_rmse = 1 - (rmse / self.rmse_norm_const)
        mean_cos = self._mean_cos_similarity(model_desc_matrix, ref_desc_matrix)
        distance = self._mean_euclidean_distance(model_desc_matrix, ref_desc_matrix)
        norm_distance = 1 - (distance / self.dist_norm_const)
        return rmse, norm_rmse, mean_cos, distance, norm_distance

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        ref_desc_matrix = context.get_extra_param_value("esm2_ref_desc_matrix")
        if ref_desc_matrix is None:
            ref_sequence = context.get_reference_chain().sequence
            ref_desc_matrix, _ = self.esm_model.query_model(
                ref_sequence, submap_indices=self.submap_indices
            )
            context.set_extra_param_value("esm2_ref_desc_matrix", ref_desc_matrix)
        model_desc_matrix = context.get_extra_param_value("esm2_model_desc_matrix")
        if model_desc_matrix is None:
            model_sequence = context.get_model_chain().sequence
            model_desc_matrix, model_contact_map = self.esm_model.query_model(
                model_sequence, submap_indices=self.submap_indices
            )
            context.set_extra_param_value("esm2_model_desc_matrix", model_desc_matrix)
            context.set_extra_param_value("esm2_predicted_contacts", model_contact_map)
        rmse, norm_rmse, mean_norm_cos, distance, norm_distance = self.do(model_desc_matrix, ref_desc_matrix)
        return {
            "rmse": rmse,
            "norm_rmse": norm_rmse, 
            "mean_norm_cos": mean_norm_cos,
            "distance": distance,
            "norm_distance": norm_distance,
        }
