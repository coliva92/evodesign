from typing import List, Optional, Tuple

import numpy as np
import numpy.typing as npt
import torch

from ..ESM2Interface import ESM2Interface
from ..ESM2ModelContainer import ESM2ModelContainer


class ESM2_v1(ESM2Interface):

    _model = None

    def __init__(
        self,
        gpu_device: Optional[str] = "cuda:0",
    ) -> None:
        self.gpu_device = gpu_device
        if ESM2_v1._model is None:
            ESM2_v1._model = ESM2ModelContainer(self.gpu_device)
        return

    def contacts_map_from_prediction(
        self,
        predicted_contacts_matrix: npt.NDArray[np.float64],
        min_separation: int = 6,
    ) -> npt.NDArray[np.int64]:
        L = predicted_contacts_matrix.shape[0]
        row_indices, column_indices = np.triu_indices(L, k=min_separation)
        pred_upper = predicted_contacts_matrix[row_indices, column_indices]
        top_L = min(L, pred_upper.size)
        # argsort can't sort in descending order;
        # highest probabilities must come first
        top_indices = np.argsort(-pred_upper)[:top_L]
        binary_map = np.zeros((L, L), dtype=np.int64)
        top_x = row_indices[top_indices]
        top_y = column_indices[top_indices]
        binary_map[top_x, top_y] = 1
        binary_map[top_y, top_x] = 1
        return binary_map

    def query_model(
        self,
        sequence: str,
        sequence_name: str = "tmp_protein",
        submap_indices: Optional[List[int]] = None,
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        data = [(sequence_name, sequence)]
        seq_ids, seqs, tokens = self._model.batch_converter(data)
        if torch.cuda.is_available() and self._model.gpu_device is not None:
            tokens = tokens.to(device=self.gpu_device, non_blocking=True)
        with torch.no_grad():
            result = self._model.esm_model(
                tokens,
                repr_layers=[self._model.esm_model.num_layers],
                return_contacts=True,
            )

        # `result['representations']` contains the weights of each layer in the
        # neural net; we only want the weights of the last layer
        last_layer = result["representations"][self._model.esm_model.num_layers]

        # the last layer contains a certain number of weights per token; a token is
        # an integer representation of each AA in the input sequence, however,
        # additional tokens are appended at the beginning and at the end of said
        # sequence; we only want to retrieve the weights corresponding to the AA
        # in the sequence
        desc_matrix = last_layer[0][1 : len(seqs[0]) + 1].cpu().numpy()

        predicted_contacts_matrix = result["contacts"][0].cpu().numpy()

        if submap_indices is not None:
            desc_matrix = desc_matrix[submap_indices]
            predicted_contacts_matrix = predicted_contacts_matrix[
                np.ix_(submap_indices, submap_indices)
            ]

        # predicted_contacts_matrix += 0.0001
        # np.fill_diagonal(predicted_contacts_matrix, 0)
        # predicted_contacts_matrix /= np.sum(predicted_contacts_matrix, axis=1)[:, np.newaxis]

        # row_idx, col_idx = np.triu_indices_from(predicted_contacts_matrix)
        # predicted_contacts = predicted_contacts_matrix[row_idx, col_idx]

        # free GPU memory
        del tokens
        del result

        return desc_matrix.astype(np.float64), predicted_contacts_matrix.astype(
            np.float64
        )
