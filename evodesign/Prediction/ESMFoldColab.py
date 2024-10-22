from .Predictor import Predictor


class ESMFoldColab(Predictor):

    def __init__(self) -> None:
        super().__init__()
        self._model = None

    def predict_pdb_str(self, sequence: str) -> str:
        import torch

        if self._model is None:
            self._model = torch.load("esmfold.model")
            self._model = self._model.eval().cuda()
            self._model.set_chunk_size(128)
        with torch.no_grad():
            prediction = self._model.infer_pdb(sequence)
        return prediction
