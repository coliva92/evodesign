from pathlib import Path

from evodesign.Prediction.AlphaFold.AlphaFold3Docker import AlphaFold3Docker
from evodesign.Prediction.DirectoryManager import DirectoryManager

MODEL_DIR = "/home/biocomp/Documents/af3_params"
MEMORY_FRACTION = 0.2


pred_path = Path.cwd() / "predictions"
dirman = DirectoryManager(pred_path, pred_path, pred_path)
af3 = AlphaFold3Docker(
    model_dir=MODEL_DIR,
    jax_preallocate_memory=MEMORY_FRACTION,
    using_rootless_docker=False,
)
af3.do(["GREETINGSTRAVELER"], dirman)
