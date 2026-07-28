from typing import Dict, List, Optional, Tuple

import numpy as np
import numpy.typing as npt

from Bio.PDB.Atom import Atom
from Bio.PDB.Superimposer import Superimposer

from .ContextInterface import ContextInterface
from .StructuralMetric import StructuralMetric




class GDT(StructuralMetric):

    def __init__(
        self,
        cutoffs: List[float] = [1.0, 2.0, 4.0, 8.0],
    ) -> None:
        super().__init__()
        self.cutoffs = cutoffs



    def do(
        self,
        model_backbone: List[Atom],
        ref_backbone: List[Atom],
        superimposer: Optional[Superimposer] = None,
        **kwargs,
    ) -> Tuple[float, float, npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        self.seed_sizes = [3, 5, 7]

        if superimposer is None:
            superimposer = Superimposer()

        total_atoms = len(ref_backbone)
        if total_atoms < 3:
            return 0.0  # SVD requires at least 3 points

        # Save original coordinates to reset between iterations
        original_coords = [atom.coord.copy() for atom in model_backbone]

        # Pre-calculate all initial seed indices (sliding windows)
        all_seeds = []
        for size in self.seed_sizes:
            for i in range(total_atoms - size + 1):
                all_seeds.append(list(range(i, i + size)))

        # Fallback if sequence is too short for the seeds
        if not all_seeds:
            all_seeds = [list(range(total_atoms))]

        scores = []
        rotation = None
        translation = None

        for cutoff in self.cutoffs:
            max_active_for_cutoff = 0

            # Test every seed to find the maximum possible set for this cutoff
            for initial_seed in all_seeds:

                # 1. Reset model coordinates for this specific seed trial
                for atom, coord in zip(model_backbone, original_coords):
                    atom.coord = coord.copy()

                active_indices = initial_seed

                while True:
                    if len(active_indices) < 3:
                        break

                    ref_subset = [ref_backbone[i] for i in active_indices]
                    model_subset = [model_backbone[i] for i in active_indices]

                    # 2. Obtain transform based on the current subset
                    superimposer.set_atoms(ref_subset, model_subset)

                    # 3. Apply to the entire model
                    superimposer.apply(model_backbone)

                    # Calculate distances
                    distances = np.array(
                        [a - b for a, b in zip(model_backbone, ref_backbone)]
                    )

                    # 4. Identify atom pairs under threshold
                    new_active_indices = [
                        i for i, d in enumerate(distances) if d <= cutoff
                    ]

                    # 5. Check for stabilization
                    if new_active_indices == active_indices:
                        break

                    active_indices = new_active_indices

                # Update the maximum found so far
                if len(active_indices) > max_active_for_cutoff:
                    max_active_for_cutoff = len(active_indices)
                    rotation = superimposer.rotran[0].copy()
                    translation = superimposer.rotran[1].copy()

            # Store the highest percentage achieved for this cutoff
            scores.append(max_active_for_cutoff / total_atoms)

        # GDT is the average of the maximums across all cutoffs
        gdt = float(np.mean(scores))
        for atom, coord in zip(model_backbone, original_coords):
            atom.coord = coord.copy()
        superimposer.rotran = (rotation, translation)
        superimposer.apply(model_backbone)
        rmsd = np.sqrt(
            np.mean(
                np.array([(a - b) ** 2 for a, b in zip(model_backbone, ref_backbone)])
            )
        )
        return gdt, rmsd, rotation, translation



    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        model_backbone = context.get_model_chain().backbone_atoms
        ref_backbone = context.get_reference_chain().backbone_atoms
        superimposer = context.get_extra_param_value("superimposer")
        if superimposer is None:
            superimposer = Superimposer()
            context.set_extra_param_value("superimposer", superimposer)
        gdt = self.do(model_backbone, ref_backbone, superimposer)
        return {"gdt": gdt}
