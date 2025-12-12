from .StructuralMetric import StructuralMetric
from Bio.PDB import Superimposer
from typing import List, Dict, Optional, Tuple
from Bio.PDB.Atom import Atom
from .ContextInterface import ContextInterface
from .Normalization.Normalization import Normalization
from .Normalization.Reciprocal import Reciprocal


class RMSD(StructuralMetric):

    # singleton superimposer
    _default_superimposer = None

    def __init__(
        self,
        normalization: Optional[Normalization] = Reciprocal(),
    ):
        super().__init__()
        self.normalization = normalization
        return

    def do(
        self,
        model_backbone: List[Atom],
        ref_backbone: List[Atom],
        superimposer: Optional[Superimposer] = None,
        **kwargs,
    ) -> Tuple[float, Optional[float]]:
        if superimposer is None:
            if self._default_superimposer is None:
                self._default_superimposer = Superimposer()
            superimposer = self._default_superimposer
        superimposer.set_atoms(ref_backbone, model_backbone)
        superimposer.apply(model_backbone)
        rmsd = superimposer.rms
        norm = None
        if self.normalization is not None:
            norm = self.normalization.do(rmsd)
        return rmsd, norm

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        model_backbone = context.get_model_chain().backbone_atoms
        ref_backbone = context.get_reference_chain().backbone_atoms
        superimposer = context.get_extra_param_value("superimposer")
        rmsd, norm = self.do(model_backbone, ref_backbone, superimposer)
        if superimposer is None:
            context.set_extra_param_value("superimposer", self._default_superimposer)
        return {
            "rmsd": rmsd,
            "norm_rmsd": norm,
        }
