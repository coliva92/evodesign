from __future__ import annotations  # noqa: N999

import numpy as np
import numpy.typing as npt

from .ModelContainer import ModelContainer


class ModelInterface:
    def __init__(
        self,
        model: ModelContainer | None = None,
        **kwargs,
    ) -> None:
        self.model = model

    def query_model(
        self, sequence: str, sequence_name: str
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        raise NotImplementedError

    def predict_pdb(self, sequence: str) -> str:
        if self.model is None:
            self.model = ModelContainer()
        return self.model.esmfold_model.infer_pdb(sequence)

    @classmethod
    def contacts_map_from_prediction(  # TODO: rename to `contact_map_from_prediction`
        cls,
        predicted_contact_map: npt.NDArray[np.float64],
        min_separation: int = 6,
        total_contacts: int | None = None,
    ) -> npt.NDArray[np.int64]:
        L = predicted_contact_map.shape[0]

        if total_contacts is None:
            total_contacts = L

        row_indices, column_indices = np.triu_indices(L, k=min_separation)
        pred_upper = predicted_contact_map[row_indices, column_indices]
        # argsort can't sort in descending order;
        # highest probabilities must come first
        top_indices = np.argsort(-pred_upper)[:total_contacts]
        binary_contact_map = np.zeros((L, L), dtype=np.int64)
        top_x = row_indices[top_indices]
        top_y = column_indices[top_indices]
        binary_contact_map[top_x, top_y] = 1
        binary_contact_map[top_y, top_x] = 1
        return binary_contact_map

    @classmethod
    def get_descriptors_submap(
        cls,
        descriptors_matrix: npt.NDArray[np.float64],
        submap_indices: list[int],
    ) -> npt.NDArray[np.float64]:
        return descriptors_matrix[submap_indices]

    @classmethod
    def get_contact_map_submap(
        cls,
        contact_map: npt.NDArray[np.float64],
        submap_indices: list[int],
    ) -> npt.NDArray[np.float64]:
        return contact_map[np.ix_(submap_indices, submap_indices)]
