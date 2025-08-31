from .Predictor import Predictor
from typing import Optional


class ESMFold(Predictor):

    def __init__(self, gpu_device: Optional[str] = "cuda:0") -> None:
        super().__init__()
        self.gpu_device = gpu_device
        self._model = None

    def predict_pdb_str(self, sequence: str) -> str:
        import torch

        if self._model is None:
            import esm

            self._model = esm.pretrained.esmfold_v1()
            self._model.eval()
            if torch.cuda.is_available() and self.gpu_device is not None:
                device = torch.device(self.gpu_device)
                self._model = self._model.to(device)
                self._model.set_chunk_size(128)
        with torch.no_grad():
            prediction = self._model.infer_pdb(sequence)
        return prediction
