import numpy as np  # noqa: N999
from Bio.PDB.Atom import Atom
from Bio.PDB.Superimposer import Superimposer

from ..ContextInterface import ContextInterface
from ..StructuralMetric import StructuralMetric


class GDT(StructuralMetric):
    def __init__(
        self,
        cutoffs: list[float] | None = None,
        sliding_window_sizes: list[int] | None = None,
    ) -> None:
        super().__init__()
        self.cutoffs = cutoffs
        if self.cutoffs is None:
            self.cutoffs = [1.0, 2.0, 4.0, 8.0]
        self.sliding_window_sizes = sliding_window_sizes
        if self.sliding_window_sizes is None:
            self.sliding_window_sizes = [3, 5, 7]
        self._rotation = None
        self._translation = None

    def do(
        self,
        model_ca_atoms: list[Atom],
        ref_ca_atoms: list[Atom],
        superimposer: Superimposer | None,
        **kwargs,
    ) -> tuple[float, ...]:

        if superimposer is None:
            superimposer = Superimposer()

        total_atoms = len(ref_ca_atoms)
        if total_atoms < 3:
            return 0.0  # SVD requires at least 3 points

        # Save original coordinates to reset between iterations
        original_coords = [atom.coord.copy() for atom in model_ca_atoms]

        # Pre-calculate all initial sliding window indices (sliding windows)
        all_windows = []
        for size in self.sliding_window_sizes:
            for i in range(total_atoms - size + 1):
                all_windows.append(list(range(i, i + size)))

        # Fallback if sequence is too short for the windows
        if not all_windows:
            all_windows = [list(range(total_atoms))]

        scores = []
        rotation = None
        translation = None

        for cutoff in self.cutoffs:
            max_active_for_cutoff = 0

            # Test every window to find the maximum possible set for this cutoff
            for initial_window in all_windows:
                # 1. Reset model coordinates for this specific seed trial
                for atom, coord in zip(model_ca_atoms, original_coords):
                    atom.coord = coord.copy()

                active_indices = initial_window

                while True:
                    if len(active_indices) < 3:
                        break

                    ref_subset = [ref_ca_atoms[i] for i in active_indices]
                    model_subset = [model_ca_atoms[i] for i in active_indices]

                    # 2. Obtain transform based on the current subset
                    superimposer.set_atoms(ref_subset, model_subset)

                    # 3. Apply to the entire model
                    superimposer.apply(model_ca_atoms)

                    # Calculate distances
                    distances = np.array(
                        [a - b for a, b in zip(model_ca_atoms, ref_ca_atoms)]
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
        for atom, coord in zip(model_ca_atoms, original_coords):
            atom.coord = coord.copy()
        superimposer.rotran = (rotation, translation)
        superimposer.apply(model_ca_atoms)
        rmsd = np.sqrt(
            np.mean(
                np.array([(a - b) ** 2 for a, b in zip(model_ca_atoms, ref_ca_atoms)])
            )
        )
        self._rotation = rotation
        self._translation = translation
        return gdt, rmsd

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        model_ca_atoms = context.get_model_chain().ca_atoms
        ref_ca_atoms = context.get_reference_chain().ca_atoms
        superimposer = context.get_extra_param_value("superimposer")

        if superimposer is None:
            superimposer = Superimposer()
            context.set_extra_param_value("superimposer", superimposer)

        gdt, rmsd = self.do(model_ca_atoms, ref_ca_atoms, superimposer)
        return {"gdt": gdt, "rmsd": rmsd}
