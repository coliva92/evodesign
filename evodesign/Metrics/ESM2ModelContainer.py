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
        if (
            ESM2ModelContainer.esmfold_model is not None
            and ESM2ModelContainer.esm_model is not None
        ):
            return

        import torch
        import esm

        ESM2ModelContainer.esmfold_model = esm.pretrained.esmfold_v1()
        ESM2ModelContainer.esmfold_model.eval()
        if torch.cuda.is_available() and self.gpu_device is not None:
            device = torch.device(self.gpu_device)
            ESM2ModelContainer.esmfold_model = ESM2ModelContainer.esmfold_model.to(
                device
            )
            ESM2ModelContainer.esmfold_model.set_chunk_size(128)
        ESM2ModelContainer.esm_model = ESM2ModelContainer.esmfold_model.esm
        ESM2ModelContainer.batch_converter = (
            ESM2ModelContainer.esmfold_model.esm_dict.get_batch_converter()
        )
        return
