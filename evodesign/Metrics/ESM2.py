from .ESM2Interface import ESM2Interface
from .ESM2ModelContainer import ESM2ModelContainer
from typing import Optional, Tuple
import torch
import numpy as np
import numpy.typing as npt
from ..Chemistry.Sequences import OPENFOLD_MAP_AMINO_ACID_TO_INT


class ESM2(ESM2Interface):

    _model = None

    def __init__(
        self,
        gpu_device: Optional[str] = "cuda:0",
    ) -> None:
        self.gpu_device = gpu_device
        if ESM2._model is None:
            ESM2._model = ESM2ModelContainer(self.gpu_device)
        return


    def contacts_map_from_prediction(self, 
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
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        L = len(sequence)
        residx = torch.from_numpy(np.array([ list(range(L)) ]))
        aa = torch.from_numpy(np.array([ OPENFOLD_MAP_AMINO_ACID_TO_INT[i] for i in sequence ]))
        data = [(sequence_name, sequence)]
        seq_ids, seqs, tokens = self._model.batch_converter(data)
        if torch.cuda.is_available() and self._model.gpu_device is not None:
            tokens = tokens.to(device=self.gpu_device, non_blocking=True)
            residx = residx.to(device=self.gpu_device, non_blocking=True)
            aa = aa.to(device=self.gpu_device, non_blocking=True)

        with torch.no_grad():
            result = self._model.esm_model(
                tokens,
                repr_layers=range(self._model.esm_model.num_layers + 1),
                return_contacts=False,
            )
            esm_s = torch.stack(
                [v for _, v in sorted(result["representations"].items())], dim=2
            )
            esm_s = esm_s[:, 1:-1]
            esm_s = esm_s.to(self._model.esmfold_model.esm_s_combine.dtype)
            esm_s = esm_s.detach()
            esm_s = (self._model.esmfold_model.esm_s_combine.softmax(0).unsqueeze(0) @ esm_s).squeeze(2)
            s_s_0 = self._model.esmfold_model.esm_s_mlp(esm_s)
            s_s_0 += self._model.esmfold_model.embedding(aa)

            s_z_0 = s_s_0.new_zeros(1, L, L, self._model.esmfold_model.cfg.trunk.pairwise_state_dim)
            recycle_z = torch.zeros_like(s_z_0)
            recycle_s = torch.zeros_like(s_s_0)
            recycle_bins = torch.zeros(*s_z_0.shape[:-1], device=self.gpu_device, dtype=torch.int64)
            recycle_s = self._model.esmfold_model.trunk.recycle_s_norm(recycle_s.detach())
            recycle_z = self._model.esmfold_model.trunk.recycle_z_norm(recycle_z.detach())
            recycle_z += self._model.esmfold_model.trunk.recycle_disto(recycle_bins.detach())
            s_s_0 += recycle_s
            z = s_z_0 + recycle_z
            s_z_0 = z + self._model.esmfold_model.trunk.pairwise_positional_embedding(residx, mask=None)

        single_rep = s_s_0.to("cpu")[0].numpy()
        pair_rep = s_z_0.to("cpu")[0].numpy()

        # free GPU memory
        del tokens
        del result
        del residx
        del aa

        return single_rep.astype(np.float64), pair_rep.astype(np.float64)
