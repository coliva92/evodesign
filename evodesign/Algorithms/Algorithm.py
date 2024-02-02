from abc import ABC, abstractmethod
from ..SettingsRetrievable import SettingsRetrievable
from ..Workspace import Workspace
from ..Random import Random





class Algorithm(SettingsRetrievable, ABC):

  def setup(self,
            targetPdb: str,
            workspaceDir: str):
    """
    Initializes the workspace and the RNG, as well as the prerequisite data
    before running the evolutinary algorithm.

    Parameters
    ----------
    targetPdb : str
        The path to the PDB file containing the target protein backbone.
    workspaceDir : str
        The folder where all the output files generated by the evolutionary
        algorithm will be stored.
    """
    # initialize the workspace
    self.workspace = Workspace(workspaceDir, targetPdb)

    # initialize the RNG
    rng = Random.generator()
    initial = self.workspace.load_rng_state(loadCheckpoint=False)
    checkpoint = self.workspace.load_rng_state()
    state = checkpoint if checkpoint else initial
    if state:
      rng.bit_generator.state = state
    else:
      self.workspace.save_rng_state(rng.bit_generator.state, checkpoint=False)
    self.workspace.save_settings(self.settings())

  

  @abstractmethod
  def __call__(self, **kwargs) -> None:
    raise NotImplementedError
