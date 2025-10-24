from .Metric import Metric
from .ContextInterface import ContextInterface
from .ESM2 import ESM2
from .Normalization.Normalization import Normalization
from .Normalization.Reciprocal import Reciprocal
from .DistanceMap import DistanceMap
import numpy as np
import numpy.typing as npt
from Bio.PDB.Atom import Atom
from typing import Optional, Dict, Tuple, List


class ESM2ContactMap(Metric):

    def __init__(
        self,
        esm_model: ESM2 = ESM2(),
        distance_threshold: float = 8.0,
        normalization: Optional[Normalization] = Reciprocal(),
    ) -> None:
        super().__init__()
        self.esm_model = esm_model
        self.distance_threshold = distance_threshold
        self.normalization = normalization

    def do(
        self,
        predicted_contacts: npt.NDArray[np.float64],
        ref_contact_map: npt.NDArray[np.float64],
        **kwargs,
    ) -> Tuple[float, Optional[float]]:
        error = np.array(
            [
                np.sqrt(np.mean((ref_contact_map[i] - predicted_contacts[i]) ** 2))
                for i in range(ref_contact_map.shape[0])
            ]
        )
        rmse = np.mean(error)
        norm = None
        if self.normalization is not None:
            norm = self.normalization.do(rmse)
        return rmse, norm

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        ref_contact_map = context.get_extra_param_value("esm2_ref_contact_map")
        if ref_contact_map is None:
            ref_atoms = context.get_reference_chain().ca_atoms
            ref_distance_map = DistanceMap.compute_map(ref_atoms)
            ref_contact_map = (
                (ref_distance_map <= self.distance_threshold)
                .astype(int)
                .astype(np.float64)
            )
            context.set_extra_param_value("esm2_ref_contact_map", ref_contact_map)
        predicted_contacts = context.get_extra_param_value("esm2_predicted_contacts")
        if predicted_contacts is None:
            model_sequence = context.get_model_chain().sequence
            model_desc_matrix, predicted_contacts = self.esm_model.query_model(
                model_sequence
            )
            context.set_extra_param_value("esm2_model_desc_matrix", model_desc_matrix)
            context.set_extra_param_value("esm2_predicted_contacts", predicted_contacts)
        rmse, norm = self.do(predicted_contacts, ref_contact_map)
        return {"rmse": rmse, "norm": norm}
