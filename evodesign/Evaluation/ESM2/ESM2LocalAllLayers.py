import numpy as np  # noqa: N999
import numpy.typing as npt
import torch

from .ModelContainer import ModelContainer
from .ModelInterface import ModelInterface


class ESM2LocalAllLayers(ModelInterface):  # chronollogically, this is ESM2Local_v2
    def query_model(
        self,
        sequence: str,
        sequence_name: str = "tmp_protein",
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:

        if self.model is None:
            self.model = ModelContainer()

        # prepare the model's input
        data = [(sequence_name, sequence)]
        _, _, tokens = self.model.batch_converter(data)

        # migrate the tokens to GPU memory if we're using the GPU
        if torch.cuda.is_available() and self.model.gpu_device is not None:
            tokens = tokens.to(device=self.model.gpu_device, non_blocking=True)

        with torch.no_grad():  # start inference
            result = self.model.esm_model(
                tokens,
                repr_layers=range(self.model.esm_model.num_layers + 1),  # all layers
                return_contacts=True,  # also predict the contact map
            )

            # reproduce the Weighted Learned Sum in ESMFold's architecture
            esm_s = torch.stack(
                [v for _, v in sorted(result["representations"].items())],
                dim=2,
            )
            esm_s = esm_s[:, 1:-1]
            esm_s = esm_s.to(self.model.esmfold_model.esm_s_combine.dtype)
            esm_s = esm_s.detach()
            esm_s = (
                self.model.esmfold_model.esm_s_combine.softmax(0).unsqueeze(0) @ esm_s
            ).squeeze(2)

        # move the queried descriptors and contact map back to CPU memory
        esm_s = esm_s.to("cpu")
        molecular_descriptors = esm_s[0].numpy()
        predicted_contact_map = result["contacts"][0].cpu().numpy()

        # free GPU memory
        del tokens
        del result

        return molecular_descriptors.astype(np.float64), predicted_contact_map.astype(
            np.float64
        )
