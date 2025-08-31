from .Metric import Metric
from Bio.PDB import Superimposer
from typing import List, Dict, Optional
from Bio.PDB.Atom import Atom
from .ContextInterface import ContextInterface


class RMSD(Metric):

    # singleton superimposer
    _default_superimposer = None

    def do(
        self,
        model_backbone: List[Atom],
        ref_backbone: List[Atom],
        superimposer: Optional[Superimposer] = None,
        **kwargs
    ) -> float:
        if superimposer is None:
            if self._default_superimposer is None:
                self._default_superimposer = Superimposer()
            superimposer = self._default_superimposer
        superimposer.set_atoms(ref_backbone, model_backbone)
        superimposer.apply(model_backbone)
        return superimposer.rms

    def do_for_fitness_fn(self, context: ContextInterface) -> Dict[str, float]:
        model_backbone = context.get_model_chain().backbone_atoms
        ref_backbone = context.get_reference_chain().backbone_atoms
        superimposer = context.get_extra_param_value("superimposer")
        rmsd = self.do(model_backbone, ref_backbone, superimposer)
        if superimposer is None:
            context.set_extra_param_value("superimposer", self._default_superimposer)
        return {"rmsd": rmsd}
