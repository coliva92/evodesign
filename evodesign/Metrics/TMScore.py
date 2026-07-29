from typing import Dict, List, Optional, Tuple

import numpy as np
import numpy.typing as npt
from Bio.PDB.Atom import Atom
from Bio.PDB.Superimposer import Superimposer

from .ContextInterface import ContextInterface
from .StructuralMetric import StructuralMetric


class TMScore(StructuralMetric):

    def normalizing_constant(
        self,
        num_residues: int,
    ) -> float:
        """Calculates the d_0 normalization distance scale."""
        if num_residues > 21:
            return 1.24 * np.cbrt(num_residues - 15) - 1.8
        return 0.5

    def _calculate_total_tm_score(
        self, model_atoms: List[Atom], ref_atoms: List[Atom], d_0: float, L_T: int
    ) -> float:
        """Calculates the overall TM-score for the superimposed structure."""
        tm_sum = 0.0
        for m, r in zip(model_atoms, ref_atoms):
            d = m - r
            tm_sum += 1.0 / (1.0 + (d / d_0) ** 2)
        return tm_sum / L_T

    def do(
        self,
        model_ca_atoms: List[Atom],
        ref_ca_atoms: List[Atom],
        superimposer: Optional[Superimposer] = None,
    ) -> Tuple[float, npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        if superimposer is None:
            superimposer = Superimposer()

        L_T = len(ref_ca_atoms)  # Length of the target protein sequence
        if L_T < 3:
            return 0.0, None, None  # Superimposer requires at least 3 points

        d_0 = self.normalizing_constant(L_T)

        # Save original coordinates to reset between iterations
        original_coords = [atom.coord.copy() for atom in model_ca_atoms]

        # Initialization
        best_tm_score = 0.0
        best_rotran = None

        # Define the different initial fragment lengths to test
        # Ensure lengths do not exceed L_T and are uniquely evaluated
        raw_lengths = [L_T, L_T // 2, L_T // 4, 4]
        fragment_lengths = sorted(
            list(set([l for l in raw_lengths if l >= 3])), reverse=True
        )

        for L_int in fragment_lengths:
            # Determine the starting positions for the fragment
            if L_int == L_T:
                start_positions = [0]
            else:
                start_positions = list(range(0, L_T - L_int + 1))

            for start in start_positions:
                # 1. Initial Superposition
                # Reset model coordinates for this specific seed trial
                for atom, coord in zip(model_ca_atoms, original_coords):
                    atom.coord = coord.copy()

                active_indices = list(range(start, start + L_int))

                # 2. Refining by Distance (Iterative Loop)
                while True:
                    if len(active_indices) < 3:
                        break  # Cannot calculate Kabsch matrix with < 3 atoms

                    ref_subset = [ref_ca_atoms[i] for i in active_indices]
                    model_subset = [model_ca_atoms[i] for i in active_indices]

                    # Calculate the Kabsch rotation matrix for this subset
                    superimposer.set_atoms(ref_subset, model_subset)
                    current_rotran = superimposer.rotran

                    # Apply the current rotation matrix to superimpose the full structures
                    superimposer.apply(model_ca_atoms)

                    # Scan and collect ALL aligned residues closer than the d_0 threshold
                    distances = np.array(
                        [
                            np.linalg.norm(a.coord - b.coord)
                            for a, b in zip(model_ca_atoms, ref_ca_atoms)
                        ]
                    )
                    new_active_indices = [i for i, d in enumerate(distances) if d < d_0]

                    # Check for convergence
                    # Note: Checking subset indices is numerically safer than comparing
                    # floating-point rotation matrices for equality. If the subsets match,
                    # the resulting matrix will be identical.
                    if new_active_indices == active_indices:
                        break

                    active_indices = new_active_indices

                # 3. Evaluate the Converged Matrix
                # The model_ca_atoms is currently superimposed using the converged matrix
                current_tm_score = self._calculate_total_tm_score(
                    model_ca_atoms, ref_ca_atoms, d_0, L_T
                )

                # 4. Update Best Matrix
                if current_tm_score > best_tm_score:
                    best_tm_score = current_tm_score
                    best_rotran = (current_rotran[0].copy(), current_rotran[1].copy())

        # Reset model backbones to their true optimal positions based on the best search result
        for atom, coord in zip(model_ca_atoms, original_coords):
            atom.coord = coord.copy()

        if best_rotran is not None:
            superimposer.rotran = best_rotran
            superimposer.apply(model_ca_atoms)
            rotation, translation = best_rotran
        else:
            rotation, translation = None, None

        return best_tm_score, rotation, translation

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        model_backbone = context.get_model_chain().ca_atoms
        ref_backbone = context.get_reference_chain().ca_atoms
        superimposer = context.get_extra_param_value("superimposer")

        if superimposer is None:
            superimposer = Superimposer()
            context.set_extra_param_value("superimposer", superimposer)

        tm_score, _, _ = self.do(model_backbone, ref_backbone, superimposer)
        return {"tm_score": tm_score}
