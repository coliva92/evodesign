import numpy as np  # noqa: N999
import numpy.typing as npt

from ..ContextInterface import ContextInterface
from ..metrics import rmse
from ..Normalization.Normalization import Normalization
from ..SequenceBasedEvaluator import SequenceBasedEvaluator
from .ModelInterface import ModelInterface


class SinglePairRepresentationSimilarity(SequenceBasedEvaluator):
    def __init__(
        self,
        model: ModelInterface | None = None,
        rmse_normalization: Normalization | None = None,
    ) -> None:
        super().__init__()
        self.model = model
        self.rmse_normalization = rmse_normalization

    def do(
        self,
        model_single_rep: npt.NDArray[np.float64],
        model_pair_rep: npt.NDArray[np.float64],
        ref_single_rep: npt.NDArray[np.float64],
        ref_pair_rep: npt.NDArray[np.float64],
        **kwargs,
    ) -> tuple[float, ...]:
        single_rmse = rmse(model_single_rep.flatten(), ref_single_rep.flatten())
        pair_rmse = rmse(model_pair_rep.flatten(), ref_pair_rep.flatten())
        rmse_val = (single_rmse + pair_rmse) / 2
        norm_rmse = rmse_val
        if self.rmse_normalization is not None:
            norm_rmse = self.rmse_normalization.do(rmse_val)
        return rmse, norm_rmse

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:

        if self.model is None:
            from .ESM2LocalSinglePairRepresentation import (
                ESM2LocalSinglePairRepresentation,
            )

            # ESM-2 model is instanced only when `query_model` is called
            self.model = ESM2LocalSinglePairRepresentation()

        ref_single_rep = context.get_extra_param_value("esmfold_ref_single_rep")
        ref_pair_rep = context.get_extra_param_value("esmfold_ref_pair_rep")
        if ref_single_rep is None or ref_pair_rep is None:
            ref_sequence = context.get_reference_chain().sequence
            ref_single_rep, ref_pair_rep = self.model.query_model(ref_sequence)
            context.set_extra_param_value("esmfold_ref_single_rep", ref_single_rep)
            context.set_extra_param_value("esmfold_ref_pair_rep", ref_pair_rep)
        model_single_rep = context.get_extra_param_value("esmfold_model_single_rep")
        model_pair_rep = context.get_extra_param_value("esmfold_model_pair_rep")
        if model_single_rep is None or model_pair_rep is None:
            if self.model is None:
                from .ESM2LocalSinglePairRepresentation import (
                    ESM2LocalSinglePairRepresentation,
                )

                self.model = ESM2LocalSinglePairRepresentation()
            model_sequence = context.get_model_chain().sequence
            model_single_rep, model_pair_rep = self.model.query_model(model_sequence)
            context.set_extra_param_value("esmfold_model_single_rep", model_single_rep)
            context.set_extra_param_value("esmfold_model_pair_rep", model_pair_rep)
        rmse, norm_rmse = self.do(
            model_single_rep, model_pair_rep, ref_single_rep, ref_pair_rep
        )
        return {
            "rmse": rmse,
            "norm_rmse": norm_rmse,
        }
