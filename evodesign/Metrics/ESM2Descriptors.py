from .NonStructuralMetric import NonStructuralMetric
from .ContextInterface import ContextInterface
from .ESM2 import ESM2
import numpy as np
import numpy.typing as npt
from typing import Optional, Dict, Tuple, List
from .Normalization.Normalization import Normalization
from .Normalization.Reciprocal import Reciprocal


class ESM2Descriptors(NonStructuralMetric):

    def __init__(
        self,
        esm_model: ESM2 = ESM2(),
        normalization: Optional[Normalization] = Reciprocal(),
        submap_indices: Optional[List[int]] = None,
    ) -> None:
        super().__init__()
        self.esm_model = esm_model
        self.normalization = normalization
        self.submap_indices = submap_indices
        return

    def do(
        self,
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
        **kwargs,
    ) -> Tuple[float, Optional[float]]:
        rmse = np.sqrt(
            np.mean((ref_desc_matrix.flatten() - model_desc_matrix.flatten()) ** 2)
        )
        norm = None
        if self.normalization is not None:
            norm = self.normalization.do(rmse)
        return rmse, norm

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
        rmse, norm = self.do(model_desc_matrix, ref_desc_matrix)
        return {
            "rmse": rmse,
            "norm_rmse": norm,
        }
