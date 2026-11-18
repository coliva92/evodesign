from Bio.PDB import Superimposer  # noqa: N999
from Bio.PDB.Atom import Atom

from ..ContextInterface import ContextInterface
from ..Normalization.Normalization import Normalization
from ..StructuralMetric import StructuralMetric


class RMSD(StructuralMetric):
    def __init__(
        self,
        normalization: Normalization | None = None,
    ):
        super().__init__()
        self.normalization = normalization

    def do(
        self,
        model_backbone: list[Atom],
        ref_backbone: list[Atom],
        superimposer: Superimposer | None = None,
        **kwargs,
    ) -> tuple[float, ...]:

        if superimposer is None:
            superimposer = Superimposer()

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
    ) -> dict[str, float]:
        model_backbone = context.get_model_chain().backbone_atoms
        ref_backbone = context.get_reference_chain().backbone_atoms
        superimposer = context.get_extra_param_value("superimposer")
        rmsd, norm = self.do(model_backbone, ref_backbone, superimposer)
        if superimposer is None:
            superimposer = Superimposer()
            context.set_extra_param_value("superimposer", superimposer)
        return {"rmsd": rmsd, "norm_rmsd": norm}
