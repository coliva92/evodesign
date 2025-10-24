from ..RetrievableSettings import RetrievableSettings
from typing import Optional, Tuple
import numpy as np
import numpy.typing as npt


class ESM2(RetrievableSettings):

    _model = None
    _batch_converter = None

    def __init__(
        self,
        gpu_device: Optional[str] = "cuda:0",
    ) -> None:
        self.gpu_device = gpu_device

    def query_model(
        self,
        sequence: str,
        sequence_name: str = "tmp_protein",
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        # initialize the model if not yet initialized
        import torch

        if self._model is None or self._batch_converter is None:
            import esm

            self._model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
            self._batch_converter = alphabet.get_batch_converter()
            self._model.eval()
            if torch.cuda.is_available() and self.gpu_device is not None:
                device = torch.device(self.gpu_device)
                self._model = self._model.to(device)
        data = [(sequence_name, sequence)]
        seq_ids, seqs, tokens = self._batch_converter(data)
        if torch.cuda.is_available() and self.gpu_device is not None:
            tokens = tokens.to(device=self.gpu_device, non_blocking=True)
        with torch.no_grad():
            result = self._model(
                tokens, repr_layers=[self._model.num_layers], return_contacts=True
            )

        # `result['representations']` contains the weights of each layer in the
        # neural net; we only want the weights of the last layer
        last_layer = result["representations"][self._model.num_layers]

        # the last layer contains a certain number of weights per token; a token is
        # an integer representation of each AA in the input sequence, however,
        # additional tokens are appended at the beginning and at the end of said
        # sequence; we only want to retrieve the weights corresponding to the AA
        # in the sequence
        desc_matrix = last_layer[0][1 : len(seqs[0]) + 1].cpu().numpy()

        predicted_contacts_matrix = result["contacts"][0].cpu().numpy()
        np.fill_diagonal(predicted_contacts_matrix, 1.0)
        row_idx, col_idx = np.triu_indices_from(predicted_contacts_matrix, k=0)
        predicted_contacts = predicted_contacts_matrix[row_idx, col_idx]

        # free GPU memory
        del tokens
        del result

        return desc_matrix, predicted_contacts
