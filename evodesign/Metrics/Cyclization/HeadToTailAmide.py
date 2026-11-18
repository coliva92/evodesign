from typing import ClassVar  # noqa: N999

from Bio.PDB.Atom import Atom

from ..ContextInterface import ContextInterface
from ..Normalization.Reciprocal import Reciprocal
from ..Normalization.ZScore import ZScore
from ..StructuralMetric import StructuralMetric


class HeadToTailAmide(StructuralMetric):
    # statistical values obtained from measuring peptide bond lengths
    # in a collection of cyclic peptide PDBs
    _z_calc: ClassVar[ZScore] = ZScore(mean_value=1.3248119, stdev=0.10498072)
    # scaling factor is to make a z-score less than 1 be normalized to larger than 0.95
    _reciprocal: ClassVar[Reciprocal] = Reciprocal(scaling_factor=0.05)

    def do(
        self,
        model_backbone: list[Atom],
        **kwargs,
    ) -> tuple[float, ...]:
        # assuming the backbone consists of atoms N-CA-C-O
        distance = model_backbone[-2] - model_backbone[0]
        z_score = self._z_calc.do(distance)
        norm_z_score = self._reciprocal.do(abs(z_score))
        return distance, z_score, norm_z_score

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        model_backbone = context.get_model_chain().backbone_atoms
        distance, z_score, norm_z_score = self.do(model_backbone)
        return {
            "cyclization": distance,
            "z_score": z_score,
            "norm_z_score": norm_z_score,
        }
