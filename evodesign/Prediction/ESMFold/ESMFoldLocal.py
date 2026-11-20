from __future__ import annotations  # noqa: N999

from ...Evaluation.ESM2.ModelContainer import ModelContainer as ESM2ModelContainer
from .ESMFoldInterface import ESMFoldInterface


class ESMFoldLocal(ESMFoldInterface):
    _model = None

    def __init__(
        self,
        model: ESM2ModelContainer | None = None,
    ) -> None:
        super().__init__()
        self.model = model

    def predict_single_pdb_str(
        self,
        sequence: str,
    ) -> str:
        if self.model is None:
            from ...Evaluation.ESM2.ESM2LocalAllLayers import ESM2LocalAllLayers

            self.model = ESM2LocalAllLayers()
        import torch

        with torch.no_grad():
            prediction = self._model.esmfold_model.infer_pdb(sequence)
        return prediction
