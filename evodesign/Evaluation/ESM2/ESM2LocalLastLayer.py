import numpy as np  # noqa: N999
import numpy.typing as npt
import torch

from .ModelContainer import ModelContainer
from .ModelInterface import ModelInterface


class ESM2LocalLastLayer(ModelInterface):
    def query_model(
        self,
        sequence: str,
        sequence_name: str = "tmp_protein",
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:

        if self.model is None:
            self.model = ModelContainer()

        # prepare the model's input
        data = [(sequence_name, sequence)]
        _, seqs, tokens = self.model.batch_converter(data)

        # migrate the tokens to GPU memory if we're using the GPU
        if torch.cuda.is_available() and self.model.gpu_device is not None:
            tokens = tokens.to(device=self.model.gpu_device, non_blocking=True)

        with torch.no_grad():  # start inference
            result = self.model.esm_model(
                tokens,
                repr_layers=[self.model.esm_model.num_layers],  # last layer only
                return_contacts=True,  # also predict the contact map
            )

        # `result['representations']` contains the weights of each layer in the
        # neural net; we only want the weights of the last layer
        last_layer = result["representations"][self.model.esm_model.num_layers]

        # the last layer contains a certain number of weights per token; a token is
        # an integer representation of each AA in the input sequence, however,
        # additional tokens are appended at the beginning and at the end of said
        # sequence; we only want to retrieve the weights corresponding to the AA
        # in the sequence
        molecular_descriptors = last_layer[0][1 : len(seqs[0]) + 1].cpu().numpy()
        predicted_contact_map = result["contacts"][0].cpu().numpy()

        # free GPU memory
        del tokens
        del result

        return molecular_descriptors.astype(np.float64), predicted_contact_map.astype(
            np.float64
        )
