from warnings import deprecated  # noqa: N999

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Superimposer import Superimposer

from ..ContextInterface import ContextInterface
from ..Normalization.Normalization import Normalization
from ..StructureBasedEvaluator import StructureBasedEvaluator


@deprecated("Use `StructureSimilarity.TMScore` instead")
class TMScore_v1(StructureBasedEvaluator):
    """
    This version of the TMScore class only considers a single structure superposition,
    despite the fact that Zhang's original TM-align program actually calculates the
    TM-score from multiple superpositions.
    """

    def __init__(
        self,
        residue_weights: list[float] | None = None,
        regularization: Normalization | None = None,
    ) -> None:
        super().__init__()
        self.regularization = regularization
        # this version of the TM-score supports evaluating the similarity of a
        # residue subset instead of the whole structure
        self.residue_weights = residue_weights
        self._backbone_weights = None

    def normalizing_constant(
        self,
        num_residues: int,
    ) -> float:
        if num_residues > 21:
            return 1.24 * np.cbrt(num_residues - 15) - 1.8
        return 0.5

    def do(
        self,
        model_backbone: list[Atom],
        ref_backbone: list[Atom],
        superimposer: Superimposer | None = None,
        **kwargs,
    ) -> tuple[float, ...]:
        if self.residue_weights is not None and self._backbone_weights is None:
            self._backbone_weights = []
            for i in range(len(self.residue_weights)):
                w = self.residue_weights[i] / 4.0
                self._backbone_weights.extend(4 * [w])
        # assume the backbones are superimposed if no superimposer is provided
        if superimposer is not None:
            superimposer.set_atoms(ref_backbone, model_backbone)
            superimposer.apply(model_backbone)
        d0 = self.normalizing_constant(len(ref_backbone))
        distances = np.array([a - b for a, b in zip(model_backbone, ref_backbone)])
        tmp = 1 / (1 + (distances / d0) ** 2)
        tm_score = np.average(tmp, weights=self._backbone_weights)
        reg_tm_score = tm_score
        if self.regularization is not None:
            reg_tm_score = self.regularization.do(tm_score)
        return tm_score, reg_tm_score

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        model_backbone = context.get_model_chain().ca_atoms
        ref_backbone = context.get_reference_chain().ca_atoms
        superimposer = context.get_persistent_param_value("superimposer")
        if superimposer is None:
            superimposer = Superimposer()
            context.set_persistent_param_value("superimposer", superimposer)
        tm_score, reg_tm_score = self.do(model_backbone, ref_backbone, superimposer)
        return {"tm_score": tm_score, "reg_tm_score": reg_tm_score}
