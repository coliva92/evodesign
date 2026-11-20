from __future__ import annotations  # noqa: N999


class ModelContainer:
    """
    Singleton container that loads and caches the pretrained ESMFold v1 model and its
    underlying ESM-2 language model in GPU memory.

    Because ESMFold bundles an ESM-2 model internally, storing these components as
    class attributes ensures they are initialized and moved to the GPU only once,
    allowing structure predictors and other interfaces share the same model instance.

    Attributes:
        esmfold_model: The pretrained ESMFold v1 model instance.
        esm_model: The internal ESM-2 language model extracted from `esmfold_model`.
        batch_converter: Tokenizer utility for converting sequences into model tokens.
        gpu_device (str | None): Target CUDA device identifier (e.g., `"cuda:0"`), or
            `None` to run on CPU.
    """

    esmfold_model = None
    esm_model = None
    batch_converter = None
    gpu_device = None

    def __init__(
        self,
        gpu_device: str | None = "cuda:0",
    ) -> None:
        """
        Initializes the ESMFold and ESM-2 models on the specified device if they have
        not already been loaded into the class attributes.
        """
        ModelContainer.gpu_device = gpu_device
        if (
            ModelContainer.esmfold_model is not None
            and ModelContainer.esm_model is not None
        ):
            return

        import esm
        import torch

        ModelContainer.esmfold_model = esm.pretrained.esmfold_v1()
        ModelContainer.esmfold_model.eval()
        if torch.cuda.is_available() and ModelContainer.gpu_device is not None:
            device = torch.device(ModelContainer.gpu_device)
            ModelContainer.esmfold_model = ModelContainer.esmfold_model.to(device)
            ModelContainer.esmfold_model.set_chunk_size(128)
        ModelContainer.esm_model = ModelContainer.esmfold_model.esm
        ModelContainer.batch_converter = (
            ModelContainer.esmfold_model.esm_dict.get_batch_converter()
        )
        return
