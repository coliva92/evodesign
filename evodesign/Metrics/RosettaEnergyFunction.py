from .Metric import Metric
import pyrosetta


class RosettaEnergyFunction(Metric):

    _score_fn = None

    def do(self, pdb_path: str, **kwargs) -> float:
        if self._score_fn is None:
            pyrosetta.init()
            self._score_fn = pyrosetta.get_score_function(True)
        pose = pyrosetta.pose_from_pdb(pdb_path)
        return self._score_fn(pose)
