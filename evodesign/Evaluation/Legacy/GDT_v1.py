from warnings import deprecated  # noqa: N999

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Superimposer import Superimposer

from ..ContextInterface import ContextInterface
from ..StructureBasedEvaluator import StructureBasedEvaluator


@deprecated("Use `StructureSimilarity.GDT` instead")
class GDT_v1(StructureBasedEvaluator):
    """
    This version of the GDT class only considers a single structure superposition,
    despite the fact that Zemla's original LGA program actually calculates the GDT from
    multiple superpositions.
    """

    def __init__(
        self,
        cutoffs: list[float] | None = None,
    ) -> None:
        super().__init__()
        self.cutoffs = cutoffs
        if self.cutoffs is None:
            self.cutoffs = [1.0, 2.0, 4.0, 8.0]

    def do(
        self,
        model_backbone: list[Atom],
        ref_backbone: list[Atom],
        superimposer: Superimposer | None = None,
        **kwargs,
    ) -> float:
        # assume the backbones are superimposed if no superimposer is provided
        if superimposer is not None:
            superimposer.set_atoms(ref_backbone, model_backbone)
            superimposer.apply(model_backbone)
        distances = np.array([a - b for a, b in zip(model_backbone, ref_backbone)])
        gdt = np.mean([np.mean([d <= c for d in distances]) for c in self.cutoffs])
        return gdt

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        model_backbone = context.get_model_chain().backbone_atoms
        ref_backbone = context.get_reference_chain().backbone_atoms
        superimposer = context.get_persistent_param_value("superimposer")
        if superimposer is None:
            superimposer = Superimposer()
            context.set_persistent_param_value("superimposer", superimposer)
        gdt = self.do(model_backbone, ref_backbone, superimposer)
        return {"gdt": gdt}
