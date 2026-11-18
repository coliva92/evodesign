import os
from pathlib import Path

from evodesign.Prediction.AlphaFold.AlphaFold3Docker import AlphaFold3Docker
from evodesign.Prediction.DirectoryManager import DirectoryManager

HOME = os.environ["HOME"]
MODEL_DIR = os.path.join(HOME, "af3_params")
MEMORY_FRACTION = 0.2


pred_path = Path.cwd() / "predictions"
dirman = DirectoryManager(pred_path, pred_path, pred_path)
af3 = AlphaFold3Docker(
    model_dir=MODEL_DIR,
    jax_preallocate_memory=MEMORY_FRACTION,
    using_rootless_docker=False,
)
af3.do(["GREETINGSTRAVELER"], dirman)
