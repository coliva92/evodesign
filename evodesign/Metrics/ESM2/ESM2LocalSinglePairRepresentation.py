import numpy as np  # noqa: N999
import numpy.typing as npt
import torch

from ...Chemistry.Sequences import OPENFOLD_MAP_AMINO_ACID_TO_INT
from .ModelContainer import ModelContainer
from .ModelInterface import ModelInterface


class ESM2LocalSinglePairRepresentation(ModelInterface):
    def query_model(
        self,
        sequence: str,
        sequence_name: str = "tmp_protein",
        **kwargs,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:

        if self.model is None:
            self.model = ModelContainer()

        # initialize the residue index representation
        L = len(sequence)
        residx = torch.from_numpy(np.array([list(range(L))]))
        aa = torch.from_numpy(
            np.array([OPENFOLD_MAP_AMINO_ACID_TO_INT[i] for i in sequence])
        )

        # prepare the model's input
        data = [(sequence_name, sequence)]
        _, _, tokens = self.model.batch_converter(data)

        # migrate the input data to GPU memory if we're using the GPU
        if torch.cuda.is_available() and self.model.gpu_device is not None:
            tokens = tokens.to(device=self.model.gpu_device, non_blocking=True)
            residx = residx.to(device=self.model.gpu_device, non_blocking=True)
            aa = aa.to(device=self.model.gpu_device, non_blocking=True)

        with torch.no_grad():  # start inference
            result = self.model.esm_model(
                tokens,
                repr_layers=range(self.model.esm_model.num_layers + 1),  # all layers
                return_contacts=False,  # no contact map prediction
            )

            # reproduce the Simplified Evoformer in ESMFold's architecture
            esm_s = torch.stack(
                [v for _, v in sorted(result["representations"].items())], dim=2
            )
            esm_s = esm_s[:, 1:-1]
            esm_s = esm_s.to(self.model.esmfold_model.esm_s_combine.dtype)
            esm_s = esm_s.detach()
            esm_s = (
                self.model.esmfold_model.esm_s_combine.softmax(0).unsqueeze(0) @ esm_s
            ).squeeze(2)
            s_s_0 = self.model.esmfold_model.esm_s_mlp(esm_s)
            s_s_0 += self.model.esmfold_model.embedding(aa)

            s_z_0 = s_s_0.new_zeros(
                1, L, L, self.model.esmfold_model.cfg.trunk.pairwise_state_dim
            )
            recycle_z = torch.zeros_like(s_z_0)
            recycle_s = torch.zeros_like(s_s_0)
            recycle_bins = torch.zeros(
                *s_z_0.shape[:-1], device=self.model.gpu_device, dtype=torch.int64
            )
            recycle_s = self.model.esmfold_model.trunk.recycle_s_norm(
                recycle_s.detach()
            )
            recycle_z = self.model.esmfold_model.trunk.recycle_z_norm(
                recycle_z.detach()
            )
            recycle_z += self.model.esmfold_model.trunk.recycle_disto(
                recycle_bins.detach()
            )
            s_s_0 += recycle_s
            z = s_z_0 + recycle_z
            s_z_0 = z + self.model.esmfold_model.trunk.pairwise_positional_embedding(
                residx, mask=None
            )

        # migrate the single and pair representations to CPU memory
        single_rep = s_s_0.to("cpu")[0].numpy()
        pair_rep = s_z_0.to("cpu")[0].numpy()

        # free GPU memory
        del tokens
        del result
        del residx
        del aa

        return single_rep.astype(np.float64), pair_rep.astype(np.float64)
