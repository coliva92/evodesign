from .Metric import Metric
from .ContextInterface import ContextInterface
from .ESM2 import ESM2
import numpy as np
import numpy.typing as npt
from typing import Optional, Dict, Tuple
from .Normalization.Normalization import Normalization
from .Normalization.Reciprocal import Reciprocal


class ESM2Descriptors(Metric):

    def __init__(
        self,
        esm_model: ESM2 = ESM2(),
        normalization: Optional[Normalization] = Reciprocal(),
    ) -> None:
        super().__init__()
        self.esm_model = esm_model
        self.normalization = normalization

    def uses_predictor(self) -> bool:
        return False

    def do(
        self,
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
        **kwargs,
    ) -> Tuple[float, Optional[float]]:
        # compute the squared error of each residue
        error = np.array(
            [
                np.sqrt(np.mean((ref_desc_matrix[i] - model_desc_matrix[i]) ** 2))
                for i in range(ref_desc_matrix.shape[0])
            ]
        )

        # then compute average over all residues; this is the same as computing
        # the weighted mean since all residues have the same number of descriptors
        rmse = np.mean(error)

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
            ref_desc_matrix, _ = self.esm_model.query_model(ref_sequence)
            context.set_extra_param_value("esm2_ref_desc_matrix", ref_desc_matrix)
        model_desc_matrix = context.get_extra_param_value("esm2_model_desc_matrix")
        if model_desc_matrix is None:
            model_sequence = context.get_model_chain().sequence
            model_desc_matrix, model_contact_map = self.esm_model.query_model(
                model_sequence
            )
            context.set_extra_param_value("esm2_model_desc_matrix", model_desc_matrix)
            context.set_extra_param_value("esm2_predicted_contacts", model_contact_map)
        rmse, norm = self.do(model_desc_matrix, ref_desc_matrix)
        return {
            "rmse": rmse,
            "norm_rmse": norm,
        }
