import numpy as np  # noqa: N999
import numpy.typing as npt

from ...Chemistry.Geometry import compute_contacts_map
from ..ContextInterface import ContextInterface
from ..metrics import euclidean_distance
from ..NonStructuralMetric import NonStructuralMetric
from ..Normalization.Normalization import Normalization
from .ModelInterface import ModelInterface


class ContactMapSpectralDistance(NonStructuralMetric):
    def __init__(
        self,
        model: ModelInterface | None = None,
        distance_threshold: float = 8.0,
        min_separation: int = 6,
        distance_normalization: Normalization | None = None,
    ) -> None:
        super().__init__()
        self.model = model
        self.distance_threshold = distance_threshold
        self.min_separation = min_separation
        self.distance_normalization = distance_normalization

    def normalized_laplacian(
        self,
        adj_matrix: npt.NDArray[np.int64],
    ) -> npt.NDArray[np.float64]:
        degrees = np.sum(adj_matrix, axis=1)
        # Compute D^(-1/2)
        # We use np.where to safely handle isolated nodes (degree = 0)
        # to avoid dividing by zero.
        with np.errstate(divide="ignore"):
            d_inv_sqrt = np.where(degrees > 0, 1.0 / np.sqrt(degrees), 0.0)
        # Create D^(-1/2) as a diagonal matrix
        D_inv_sqrt_mat = np.diag(d_inv_sqrt)
        # Create the Identity matrix (I)
        I = np.eye(adj_matrix.shape[0])
        # Compute L_norm = I - D^(-1/2) * A * D^(-1/2)
        L_norm = I - np.dot(D_inv_sqrt_mat, np.dot(adj_matrix, D_inv_sqrt_mat))
        return L_norm

    def spectral_distance(
        self,
        predicted_laplacian: npt.NDArray[np.int64],
        ref_laplacian: npt.NDArray[np.int64],
    ) -> float:
        assert predicted_laplacian.shape == ref_laplacian.shape
        eig_predicted = np.linalg.eigvalsh(predicted_laplacian)
        eig_ref = np.linalg.eigvalsh(ref_laplacian)
        eig_predicted /= 2.0
        eig_ref /= 2.0
        # Compute Euclidean distance between the scaled vectors
        distance = euclidean_distance(eig_predicted, eig_ref)
        return distance

    def do(
        self,
        predicted_laplacian: npt.NDArray[np.int64],
        ref_laplacian: npt.NDArray[np.int64],
        **kwargs,
    ) -> tuple[float, ...]:
        distance = self.spectral_distance(predicted_laplacian, ref_laplacian)
        norm_distance = distance
        if self.distance_normalization is not None:
            norm_distance = self.distance_normalization.do(distance)
        return distance, norm_distance

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        ref_laplacian = context.get_extra_param_value("ref_norm_laplacian")
        if ref_laplacian is None:
            ref_contact_map = context.get_extra_param_value("ref_contact_map")
            if ref_contact_map is None:
                ref_atoms = context.get_reference_chain().ca_atoms
                ref_contact_map = compute_contacts_map(
                    ref_atoms,
                    self.distance_threshold,
                    self.min_separation,
                )
                context.set_extra_param_value("ref_contact_map", ref_contact_map)
            ref_laplacian = self.normalized_laplacian(ref_contact_map)
            context.set_extra_param_value("ref_norm_laplacian", ref_laplacian)
        predicted_laplacian = context.get_extra_param_value(
            "esm2_predicted_norm_laplacian"
        )
        if predicted_laplacian is None:
            predicted_map = context.get_extra_param_value("esm2_predicted_contacts")
            if predicted_map is None:
                if self.model is None:
                    from .ESM2LocalAllLayers import ESM2LocalAllLayers

                    self.model = ESM2LocalAllLayers()
                model_sequence = context.get_model_chain().sequence
                model_desc_matrix, predicted_map = self.model.query_model(
                    model_sequence
                )
                context.set_extra_param_value(
                    "esm2_model_desc_matrix", model_desc_matrix
                )
                context.set_extra_param_value("esm2_predicted_contacts", predicted_map)
            predicted_map = self.model.contacts_map_from_prediction(
                predicted_map,
                self.min_separation,
            )
            predicted_laplacian = self.normalized_laplacian(predicted_map)
            context.set_extra_param_value(
                "esm2_predicted_norm_laplacian", predicted_laplacian
            )
        dist, norm_dist = self.do(predicted_laplacian, ref_laplacian)
        return {
            "distance": dist,
            "norm_distance": norm_dist,
        }
