from .ESMFoldInterface import ESMFoldInterface
from ..Metrics.ESM2ModelContainer import ESM2ModelContainer
from typing import Optional


class ESMFold(ESMFoldInterface):

    _model = None

    def __init__(
        self,
        gpu_device: Optional[str] = "cuda:0",
    ) -> None:
        super().__init__()
        self.gpu_device = gpu_device
        if ESMFold._model is None:
            ESMFold._model = ESM2ModelContainer(self.gpu_device)
        return

    def predict_single_pdb_str(
        self,
        sequence: str,
    ) -> str:
        import torch

        with torch.no_grad():
            prediction = self._model.esmfold_model.infer_pdb(sequence)
        return prediction
