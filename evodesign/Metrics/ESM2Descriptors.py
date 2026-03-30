from .NonStructuralMetric import NonStructuralMetric
from .ContextInterface import ContextInterface
from .ESM2 import ESM2
import numpy as np
import numpy.typing as npt
from typing import Optional, Dict, List, Tuple
from .Normalization.Normalization import Normalization
from .Normalization.Reciprocal import Reciprocal
from scipy.stats import entropy


class ESM2Descriptors(NonStructuralMetric):

    def __init__(
        self,
        esm_model: ESM2 = ESM2(),
        submap_indices: Optional[List[int]] = None,
        normalization: Optional[Normalization] = None,
        upper_bound: float = 5488.444555777253
    ) -> None:
        super().__init__()
        self.esm_model = esm_model
        self.submap_indices = submap_indices
        self.normalization = normalization
        self.upper_bound = upper_bound
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
            values.append(self._cos_similarity(u, v))
        return np.mean(values)
    
    def _mean_distance(
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

    def _kullback_leibler(
        self,
        u: npt.NDArray[np.float64],
        v: npt.NDArray[np.float64],
        pseudo_count: float = 0.5
    ) -> float:
        min_value = min([ u.min(), v.min() ])
        max_value = max([ u.max(), v.max() ])
        dist_range = (min_value, max_value)
        a, _ = np.histogram(u, 
                            bins=20, 
                            range=dist_range, 
                            density=False)
        b, _ = np.histogram(v, 
                            bins=20, 
                            range=dist_range,
                            density=False)
        a += pseudo_count
        b += pseudo_count
        divergence = entropy(a / a.sum(), b / b.sum())
        return divergence

    def _symmetric_kullback_leibler(
        self,
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        u = model_desc_matrix.flatten()
        v = ref_desc_matrix.flatten()
        a = self._kullback_leibler(u, v)
        b = self._kullback_leibler(v, u)
        return a + b

    def do(
        self,
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
        **kwargs,
    ) -> Tuple[float]:
        # rms = self._rmse(model_desc_matrix, ref_desc_matrix)
        distance = self._mean_distance(model_desc_matrix, ref_desc_matrix)
        # mean_cos = self._mean_cos_similarity(model_desc_matrix, ref_desc_matrix)
        # divergence = self._symmetric_kullback_leibler(model_desc_matrix, ref_desc_matrix)
        # return rms, mean_cos, divergence, distance
        norm_distance = 1 - (np.exp(distance) / self.upper_bound)
        return distance, norm_distance

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
        distance, norm_distance = self.do(model_desc_matrix, ref_desc_matrix)
        return {
            "distance": distance,
            "norm_distance": norm_distance,
        }
