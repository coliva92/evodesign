from __future__ import annotations  # noqa: N999

from ..ContextInterface import ContextInterface
from ..StructureBasedEvaluator import StructureBasedEvaluator


class RosettaEnergyFunction(StructureBasedEvaluator):
    _score_fn = None

    def do(
        self,
        pdb_path: str,
        **kwargs,
    ) -> float:
        if self._score_fn is None:
            import pyrosetta

            pyrosetta.init()
            self._score_fn = pyrosetta.get_score_function(True)
        pose = pyrosetta.pose_from_pdb(pdb_path)
        return self._score_fn(pose)

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        pdb_path = context.get_model_chain().pdb_path
        energy_score = self.do(pdb_path)
        return {"energy_score": energy_score}
