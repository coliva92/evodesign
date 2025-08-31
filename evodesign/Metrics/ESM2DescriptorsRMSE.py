from .Metric import Metric
from .ContextInterface import ContextInterface
import numpy as np
import numpy.typing as npt
from typing import Optional, Dict


class ESM2DescriptorsRMSE(Metric):

    _model = None
    _batch_converter = None

    def __init__(self, gpu_device: Optional[str] = "cuda:0") -> None:
        super().__init__()
        self.gpu_device = gpu_device

    def do(
        self, model_sequence: str, ref_desc_matrix: npt.NDArray[np.float64], **kwargs
    ) -> float:
        # compute the RMS of each residue
        model_desc_matrix = self.compute_descriptors_matrix(model_sequence)
        rmse = np.array(
            [
                np.sqrt(np.mean((ref_desc_matrix[i] - model_desc_matrix[i]) ** 2))
                for i in range(ref_desc_matrix.shape[0])
            ]
        )

        # then compute the average of the RMS of all residues; this is the same as computing
        # the weighted mean since all residues have the same number of descriptors
        return np.mean(rmse)

    def do_for_fitness_fn(self, context: ContextInterface) -> Dict[str, float]:
        model_sequence = context.get_model_chain().sequence
        ref_desc_matrix = context.get_extra_param_value("reference_esm2_descriptors")
        if ref_desc_matrix is None:
            ref_sequence = context.get_reference_chain().sequence
            ref_desc_matrix = self.compute_descriptors_matrix(ref_sequence)
            context.set_extra_param_value("reference_esm2_descriptors", ref_desc_matrix)
        rmse = self.do(model_sequence, ref_desc_matrix)
        return {"rmse": rmse}

    def compute_descriptors_matrix(
        self, sequence: str, sequence_name: str = "temp_protein"
    ) -> npt.NDArray[np.float64]:
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
                tokens, repr_layers=[self._model.num_layers], return_contacts=False
            )

        # `result['representations']` contains the weights of each layer in the
        # neural net; we only want the weights of the last layer
        last_layer = result["representations"][self._model.num_layers]

        # the last layer contains a certain number of weights per token; a token is
        # an integer representation of each AA in the input sequence, however,
        # additional tokens are appended at the beginning and at the end of said
        # sequence; we only want to retrieve the weights corresponding to the AA
        # in the sequence
        matrix = last_layer[0][1 : len(seqs[0]) + 1].cpu().numpy()

        # free GPU memory
        del tokens
        del result

        return matrix
