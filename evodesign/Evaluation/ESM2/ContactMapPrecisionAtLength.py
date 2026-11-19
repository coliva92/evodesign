import numpy as np  # noqa: N999
import numpy.typing as npt

from ...Protein.Geometry import compute_contacts_map
from ..ContextInterface import ContextInterface
from ..SequenceBasedEvaluator import SequenceBasedEvaluator
from .ModelInterface import ModelInterface


class ContactMapPrecisionAtLength(SequenceBasedEvaluator):
    def __init__(
        self,
        model: ModelInterface | None = None,
        distance_threshold: float = 8.0,
        min_separation: int = 6,
    ) -> None:
        super().__init__()
        self.model = model
        self.distance_threshold = distance_threshold
        self.min_separation = min_separation

    def do(
        self,
        predicted_map: npt.NDArray[np.int64],
        reference_map: npt.NDArray[np.int64],
        total_contacts: int,
    ) -> tuple[float, ...]:
        assert predicted_map.shape[0] == predicted_map.shape[1]
        assert predicted_map.shape == reference_map.shape
        L = reference_map.shape[0]
        num_hits = 0
        for i in range(L):
            for j in range(i + 1, L):
                if predicted_map[i, j] == 1 and reference_map[i, j] == 1:
                    num_hits += 1
        return num_hits / total_contacts

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        ref_contact_map = context.get_extra_param_value("ref_contact_map")
        if ref_contact_map is None:
            ref_atoms = context.get_reference_chain().ca_atoms
            ref_contact_map = compute_contacts_map(
                ref_atoms,
                self.distance_threshold,
                self.min_separation,
            )
            context.set_extra_param_value("ref_contact_map", ref_contact_map)
        total_contacts = context.get_extra_param_value("num_ref_contacts")
        if total_contacts is None:
            total_contacts = np.triu(ref_contact_map, k=1).sum()
            context.set_extra_param_value("num_ref_contacts", total_contacts)
        predicted_map = context.get_extra_param_value("esm2_predicted_contacts")
        if predicted_map is None:
            if self.model is None:
                from .ESM2LocalAllLayers import ESM2LocalAllLayers

                self.model = ESM2LocalAllLayers()
            model_sequence = context.get_model_chain().sequence
            model_desc_matrix, predicted_map = self.model.query_model(model_sequence)
            context.set_extra_param_value("esm2_model_desc_matrix", model_desc_matrix)
            context.set_extra_param_value("esm2_predicted_contacts", predicted_map)
        predicted_map = self.model.contacts_map_from_prediction(
            predicted_map,
            self.min_separation,
            total_contacts,
        )
        P_at_L = self.do(predicted_map, ref_contact_map, total_contacts)
        return {"precision": P_at_L}
