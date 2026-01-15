from typing import Optional


class ESM2ModelContainer:

    esmfold_model = None
    esm_model = None
    batch_converter = None

    def __init__(
        self,
        gpu_device: Optional[str] = "cuda:0",
    ) -> None:
        self.gpu_device = gpu_device
        if self.esmfold_model is not None and self.esm_model is not None:
            return

        import torch
        import esm

        self.esmfold_model = esm.pretrained.esmfold_v1()
        self.esmfold_model.eval()
        if torch.cuda.is_available() and self.gpu_device is not None:
            device = torch.device(self.gpu_device)
            self.esmfold_model = self.esmfold_model.to(device)
            self.esmfold_model.set_chunk_size(128)
        self.esm_model = self.esmfold_model.esm
        self.batch_converter = self.esmfold_model.esm_dict.get_batch_converter()
        return
