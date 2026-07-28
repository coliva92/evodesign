from typing import Dict, Tuple

import numpy as np
import numpy.typing as npt

from ..Chemistry.Geometry import compute_contacts_map
from .ContextInterface import ContextInterface
from .ESM2 import ESM2
from .NonStructuralMetric import NonStructuralMetric




class ESM2ContactMapRealRef(NonStructuralMetric):

    def __init__(
        self,
        esm_model: ESM2 = ESM2(),
        distance_threshold: float = 8.0,
        min_separation: int = 6,
        dist_norm_const: float = 1.0,
    ) -> None:
        super().__init__()
        self.esm_model = esm_model
        self.distance_threshold = distance_threshold
        self.min_separation = min_separation
        self.dist_norm_const = dist_norm_const
        return



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
        distance = np.linalg.norm(eig_predicted - eig_ref)
        return distance



    def do(
        self,
        predicted_laplacian: npt.NDArray[np.int64],
        ref_laplacian: npt.NDArray[np.int64],
        **kwargs,
    ) -> Tuple[float, float]:
        dist = self.spectral_distance(predicted_laplacian, ref_laplacian)
        norm_dist = 1 - (dist / self.dist_norm_const)
        return dist, norm_dist



    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        ref_laplacian = context.get_extra_param_value("ref_norm_laplacian")
        if ref_laplacian is None:
            ref_contact_map = context.get_extra_param_value("ref_contact_map")
            if ref_contact_map is None:
                ref_atoms = context.get_reference_chain().ca_atoms
                ref_contact_map = compute_contacts_map(
                    ref_atoms, self.distance_threshold, self.min_separation
                )
                context.set_extra_param_value("ref_contact_map", ref_contact_map)
            ref_laplacian = self.normalized_laplacian(ref_contact_map)
            context.set_extra_param_value("ref_norm_laplacian", ref_laplacian)
        predicted_laplacian = context.get_extra_param_value(
            "esm2_predicted_norm_laplacian"
        )
        if predicted_laplacian is None:
            predicted_contacts = context.get_extra_param_value(
                "esm2_predicted_contacts"
            )
            if predicted_contacts is None:
                model_sequence = context.get_model_chain().sequence
                model_desc_matrix, predicted_contacts = self.esm_model.query_model(
                    model_sequence
                )
                context.set_extra_param_value(
                    "esm2_model_desc_matrix", model_desc_matrix
                )
                context.set_extra_param_value(
                    "esm2_predicted_contacts", predicted_contacts
                )
            predicted_contacts = self.esm_model.contacts_map_from_prediction(
                predicted_contacts, self.min_separation
            )
            predicted_laplacian = self.normalized_laplacian(predicted_contacts)
            context.set_extra_param_value(
                "esm2_predicted_norm_laplacian", predicted_laplacian
            )
        dist, norm_dist = self.do(predicted_laplacian, ref_laplacian)
        return {
            "distance": dist,
            "norm_distance": norm_dist,
        }
