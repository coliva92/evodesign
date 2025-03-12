from .Metric import Metric
from Bio.PDB import Superimposer
from typing import List
from Bio.PDB.Atom import Atom


class RMSD(Metric):

    _superimposer = None

    def do(
        self,
        model_backbone: List[Atom],
        ref_backbone: List[Atom],
        **kwargs
    ) -> float:
        if not self._superimposer:
            self._superimposer = Superimposer()
        self._superimposer.set_atoms(ref_backbone, model_backbone)
        self._superimposer.apply(model_backbone)
        return self._superimposer.rms
