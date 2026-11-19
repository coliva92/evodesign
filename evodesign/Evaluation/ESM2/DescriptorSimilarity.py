import numpy as np  # noqa: N999
import numpy.typing as npt

from ..ContextInterface import ContextInterface
from ..metrics import euclidean_distance, norm_cosine_similarity, rmse
from ..Normalization.Normalization import Normalization
from ..SequenceBasedEvaluator import SequenceBasedEvaluator
from .ModelInterface import ModelInterface


class DescriptorSimilarity(SequenceBasedEvaluator):
    def __init__(
        self,
        model: ModelInterface | None,
        submap_indices: list[int] | None = None,
        rmse_normalization: Normalization | None = None,
        distance_normalization: Normalization | None = None,
    ) -> None:
        super().__init__()
        self.model = model
        self.submap_indices = submap_indices
        self.rmse_normalization = rmse_normalization
        self.distance_normalization = distance_normalization

    def _mean_cos_similarity(
        self,
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
    ) -> float:
        values = []
        for i in range(model_desc_matrix.shape[0]):
            values.append(
                norm_cosine_similarity(
                    model_desc_matrix[i, :],
                    ref_desc_matrix[i, :],
                )
            )
        return np.mean(values)

    def _mean_euclidean_distance(
        self,
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
    ) -> float:
        values = []
        for i in range(model_desc_matrix.shape[0]):
            values.append(
                euclidean_distance(
                    model_desc_matrix[i, :],
                    ref_desc_matrix[i, :],
                )
            )
        return np.mean(values)

    def do(
        self,
        model_desc_matrix: npt.NDArray[np.float64],
        ref_desc_matrix: npt.NDArray[np.float64],
        **kwargs,
    ) -> tuple[float, ...]:
        rmse_val = rmse(model_desc_matrix.flatten(), ref_desc_matrix.flatten())
        norm_rmse = rmse_val
        if self.rmse_normalization is not None:
            norm_rmse = self.rmse_normalization.do(rmse_val)
        mean_cos = self._mean_cos_similarity(model_desc_matrix, ref_desc_matrix)
        distance = self._mean_euclidean_distance(model_desc_matrix, ref_desc_matrix)
        norm_distance = distance
        if self.distance_normalization is not None:
            norm_distance = self.distance_normalization.do(distance)
        return rmse_val, norm_rmse, mean_cos, distance, norm_distance

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:

        if self.model is None:
            from .ESM2LocalLastLayer import ESM2LocalLastLayer

            # ESM-2 model is instanced only when `query_model` is called
            self.model = ESM2LocalLastLayer()

        ref_desc_matrix = context.get_persistent_param_value("esm2_ref_desc_matrix")
        if ref_desc_matrix is None:
            ref_sequence = context.get_reference_chain().sequence
            ref_desc_matrix, _ = self.model.query_model(ref_sequence)
            context.set_persistent_param_value("esm2_ref_desc_matrix", ref_desc_matrix)
        model_desc_matrix = context.get_extra_param_value("esm2_model_desc_matrix")
        if model_desc_matrix is None:
            model_sequence = context.get_model_chain().sequence
            model_desc_matrix, model_contact_map = self.model.query_model(
                model_sequence
            )
            context.set_extra_param_value("esm2_model_desc_matrix", model_desc_matrix)
            context.set_extra_param_value("esm2_predicted_contacts", model_contact_map)

        # get the requested subset of descriptors
        if self.submap_indices is not None:
            ref_desc_matrix = self.model.get_descriptors_submap(
                ref_desc_matrix,
                self.submap_indices,
            )
            model_desc_matrix = self.model.get_descriptors_submap(
                model_desc_matrix,
                self.submap_indices,
            )

        rmse_val, norm_rmse, mean_norm_cos, distance, norm_distance = self.do(
            model_desc_matrix,
            ref_desc_matrix,
        )
        return {
            "rmse": rmse_val,
            "norm_rmse": norm_rmse,
            "mean_norm_cos": mean_norm_cos,
            "distance": distance,
            "norm_distance": norm_distance,
        }
