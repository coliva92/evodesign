from .NonStructuralMetric import NonStructuralMetric
from .ContextInterface import ContextInterface
from .ESM2 import ESM2
import numpy as np
import numpy.typing as npt
from typing import Optional, Dict, List, Tuple
from .Normalization.Normalization import Normalization
from .Normalization.Reciprocal import Reciprocal


class ESM2Descriptors(NonStructuralMetric):

    def __init__(
        self,
        esm_model: ESM2 = ESM2(),
        submap_indices: Optional[List[int]] = None,
        normalization: Optional[Normalization] = Reciprocal(),
    ) -> None:
        super().__init__()
        self.esm_model = esm_model
        self.submap_indices = submap_indices
        self.normalization = normalization
        return

    def do(
        self,
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
        **kwargs,
    ) -> Tuple[float, float, float]:
        distance = np.mean([
            np.linalg.norm(model_desc_matrix[i, :] - ref_desc_matrix[i, :])
            for i in range(model_desc_matrix.shape[0])    
        ])
        if self.normalization is not None:
            norm_distance = self.normalization.do(distance)
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
