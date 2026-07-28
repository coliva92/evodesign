from typing import Dict, Tuple

import numpy as np
import numpy.typing as npt

from .ContextInterface import ContextInterface
from .ESM2 import ESM2_3
from .NonStructuralMetric import NonStructuralMetric




class ESMSinglePairRepresentation(NonStructuralMetric):

    def __init__(
        self,
        rmse_norm_const: float = 1.0,
        esm_model=None,
    ) -> None:
        super().__init__()
        self.esm_model = ESM2_3()  # hardcoded for now
        self.rmse_norm_const = rmse_norm_const
        return



    def _rmse(
        self,
        a: npt.NDArray[np.float64],
        b: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        u = a.flatten()
        v = b.flatten()
        return np.sqrt(np.mean((u - v) ** 2))



    def do(
        self,
        model_single_rep: npt.NDArray[np.float64],
        model_pair_rep: npt.NDArray[np.float64],
        ref_single_rep: npt.NDArray[np.float64],
        ref_pair_rep: npt.NDArray[np.float64],
        **kwargs,
    ) -> Tuple[float, float]:
        single_rmse = self._rmse(model_single_rep, ref_single_rep)
        pair_rmse = self._rmse(model_pair_rep, ref_pair_rep)
        rmse = (single_rmse + pair_rmse) / 2
        norm_rmse = 1 - (rmse / self.rmse_norm_const)
        return rmse, norm_rmse



    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        ref_single_rep = context.get_extra_param_value("esmfold_ref_single_rep")
        ref_pair_rep = context.get_extra_param_value("esmfold_ref_pair_rep")
        if ref_single_rep is None or ref_pair_rep is None:
            ref_sequence = context.get_reference_chain().sequence
            ref_single_rep, ref_pair_rep = self.esm_model.query_model(ref_sequence)
            context.set_extra_param_value("esmfold_ref_single_rep", ref_single_rep)
            context.set_extra_param_value("esmfold_ref_pair_rep", ref_pair_rep)
        model_single_rep = context.get_extra_param_value("esmfold_model_single_rep")
        model_pair_rep = context.get_extra_param_value("esmfold_model_pair_rep")
        if model_single_rep is None or model_pair_rep is None:
            model_sequence = context.get_model_chain().sequence
            model_single_rep, model_pair_rep = self.esm_model.query_model(
                model_sequence
            )
            context.set_extra_param_value("esmfold_model_single_rep", model_single_rep)
            context.set_extra_param_value("esmfold_model_pair_rep", model_pair_rep)
        rmse, norm_rmse = self.do(
            model_single_rep, model_pair_rep, ref_single_rep, ref_pair_rep
        )
        return {
            "rmse": rmse,
            "norm_rmse": norm_rmse,
        }
