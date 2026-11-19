# EvoDesign Refactoring Plan: Temporary Storage, Predictor Decoupling, and Evaluation Context

This document captures the proposed refactorings for:
1. Introducing a dedicated `TempStorage` capability to separate ephemeral temporary files from persistent job storage.
2. Removing `evodesign/Prediction/DirectoryManager.py` and decoupling structure predictors so they work cleanly both standalone and inside the optimization loop.
3. Refactoring `evodesign/Metrics/Context.py` and `ContextInterface.py` to support typed capabilities, composite cache keys, per-run reference caching, and shared domain cache helpers.
4. Small related cleanups and bug fixes identified during the analysis.

---

## 1. Ephemeral Files vs. Persistent Storage (`TempStorage`)

### Problem
Currently, there is no clean mechanism for components (such as `iLearnDescriptors` or external bioinformatics tools) to create temporary files without violating Separation of Concerns and the Interface Segregation Principle (ISP):
- Passing `PathsContainer` into `Context` exposes the entire persistent job layout (`results_npz_path`, `pymoo_algorithm_bin_path`, `settings_json_path`, etc.) to evaluators that have no business knowing about checkpoints or output paths.
- Passing `StorageManager` into `Context` couples evaluators to a PyMOO `Callback` responsible for job persistence and exposes mutating persistence/deletion methods to evaluation code.
- Constructor injection on evaluators is awkward because `evodesign/Settings.py` instantiates evaluators from JSON *before* `PathsContainer` exists, and `_SettingsEncoder` serializes public instance attributes back into `settings.json`.
- `StorageManager` currently owns `delete_file`, `delete_folder`, and `delete_non_essential_files_and_folders()` solely to sweep up temporary folders dumped directly into `PathsContainer.path`.

There are two fundamentally different categories of I/O in `evodesign`:
1. **Persistent job artifacts** (`results.npz`, `pymoo_algorithm.bin`, RNG states, `settings.json`, `version.txt`, `aa_profile.txt`, copied target PDB) — owned by `StorageManager`.
2. **Ephemeral scratch / IPC files** (temporary FASTA files for `iLearn`, AlphaFold input JSONs / raw CIFs / intermediate PDBs) — short-lived files that only exist because external tools require file paths on disk.

### Solution: Introduce `evodesign/System/TempStorage.py`

Create a small, focused capability class backed by Python's `tempfile` module.
- When given a `base_dir` (e.g., `<job_dir>/tmp` during an optimization run), all temporary files and subdirectories are created there (helpful on HPC clusters with limited `/tmp`, for Docker bind mounts, and for centralized crash cleanup).
- When `base_dir=None` (standalone library use or unit tests), it transparently falls back to the OS temporary directory (`$TMPDIR` / `/tmp`).
- Every temporary file/directory is managed by a context manager for deterministic cleanup and concurrency safety (no colliding hardcoded filenames like `tmp_sequence.fasta`).

```python
from contextlib import contextmanager
import os
import shutil
import tempfile
from typing import Iterator


class TempStorage:
    """Scoped provider for ephemeral files and directories."""

    def __init__(self, base_dir: str | None = None) -> None:
        self._base_dir = os.path.abspath(base_dir) if base_dir is not None else None

    def _ensure_base_dir(self) -> str | None:
        if self._base_dir is not None:
            os.makedirs(self._base_dir, exist_ok=True)
        return self._base_dir

    @contextmanager
    def temp_file(self, suffix: str = "", prefix: str = "tmp_") -> Iterator[str]:
        fd, path = tempfile.mkstemp(
            suffix=suffix,
            prefix=prefix,
            dir=self._ensure_base_dir(),
        )
        os.close(fd)
        try:
            yield path
        finally:
            if os.path.exists(path):
                os.remove(path)

    @contextmanager
    def temp_dir(self, prefix: str = "tmp_") -> Iterator[str]:
        with tempfile.TemporaryDirectory(
            prefix=prefix,
            dir=self._ensure_base_dir(),
        ) as dir_path:
            yield dir_path

    def cleanup(self) -> None:
        """Remove the entire managed temp directory (e.g., on job completion or error)."""
        if self._base_dir is not None and os.path.isdir(self._base_dir):
            shutil.rmtree(self._base_dir)
```

*(Note: For completely self-contained single-call temporary files that never need to live under the job directory, using `tempfile` directly is also valid; `TempStorage` unifies both cases because `TempStorage()` defaults to system `tempfile`.)*

### Updates to `PathsContainer`, `StorageManager`, and `__main__.py`

1. **`evodesign/System/PathsContainer.py`**:
   - Replace `prediction_pdbs_dir`, `predictor_input_dir`, and `predictor_output_dir` with a single `temp_dir: str` (e.g., `os.path.join(base_path, "tmp")`).
   - Remove `essential_files: frozenset[str]` — it is no longer needed once temporary files are isolated in `temp_dir`.
2. **`evodesign/Callbacks/StorageManager.py`**:
   - Remove `self.predictor_directory`.
   - Remove `delete_file`, `delete_folder`, and `delete_non_essential_files_and_folders`.
   - `StorageManager` becomes strictly responsible for saving and loading persistent run data and checkpoints.
   - *(Bug fix bonus: previously, `PathsContainer.essential_files` omitted the target PDB copied by `StorageManager.save_target_pdb()`, causing `delete_non_essential_files_and_folders()` to delete it at job completion.)*
3. **`evodesign/__main__.py` (Composition Root)**:
   - Construct `paths = PathsContainer.create(args.output_dir, args.jobname)` once.
   - Instantiate `storage = StorageManager(paths, ...)` and `temp_storage = TempStorage(paths.temp_dir)` side-by-side.
   - Pass `temp_storage` into `algorithm_factory.create_problem(ref_chain, temp_storage, aa_profile)`.
   - Replace `storage.delete_non_essential_files_and_folders()` calls on error/retry/completion with `temp_storage.cleanup()`.

---

## 2. Remove `Prediction/DirectoryManager.py` & Decouple `Predictor`

### Problem
`evodesign/Prediction/DirectoryManager.py` couples callers to internal predictor staging details:
- It mixes two distinct concerns: the caller's desired output location (`prediction_pdbs_dir`) and backend-specific staging folders (`model_input_dir`, `model_output_dir`).
- Backends like `ESMFoldInterface` and `Null` do not even use `model_input_dir` or `model_output_dir`.
- Standalone usage (e.g., importing `AlphaFold3Docker` outside the optimization loop) forces the caller to construct a `DirectoryManager` with three directories even when they only want a single PDB file or PDB string.
- Because `Predictor.do(...)` returns `None`, `SingleObjectiveCPD._predict_structures` has to know the internal filename convention (`f"{prefix}_{i}.pdb"`) to locate the generated files.
- Predictor-specific directory paths leak all the way up into `PathsContainer` and `StorageManager`.

### Solution: Output-Oriented `Predictor` API + Optional `TempStorage` Injection

Separate **what the caller wants (output)** from **how the predictor stages intermediates (private `TempStorage`)**, using **Optional Dependency Injection with a Sensible Default** (`temp_storage: TempStorage | None = None`).

#### 1. Updated `evodesign/Prediction/Predictor.py`

```python
import os
from abc import ABC, abstractmethod

from requests.exceptions import ConnectTimeout

from ..System.Exceptions import HttpGatewayTimeout
from ..System.TempStorage import TempStorage


class Predictor(ABC):
    @abstractmethod
    def predict_single_pdb_file(
        self,
        sequence: str,
        output_pdb_path: str,
        temp_storage: TempStorage | None = None,
    ) -> str:
        """Predict structure for `sequence`, write PDB to `output_pdb_path`, and return the path."""
        raise NotImplementedError

    def predict_single_pdb_str(
        self,
        sequence: str,
        temp_storage: TempStorage | None = None,
    ) -> str:
        """Standalone convenience: return raw PDB text without requiring an output path."""
        storage = temp_storage or TempStorage()
        with storage.temp_file(suffix=".pdb") as tmp_pdb:
            self.predict_single_pdb_file(sequence, tmp_pdb, temp_storage=storage)
            with open(tmp_pdb, "rt", encoding="utf-8") as f:
                return f.read()

    def do(
        self,
        sequences: list[str],
        output_dir: str,
        temp_storage: TempStorage | None = None,
    ) -> list[str]:
        """Predict a batch of sequences into `output_dir` and return the ordered list of PDB paths."""
        os.makedirs(output_dir, exist_ok=True)
        pdb_paths: list[str] = []
        for i, sequence in enumerate(sequences):
            out_path = os.path.join(output_dir, f"prediction_{i}.pdb")
            while True:
                try:
                    self.predict_single_pdb_file(
                        sequence, out_path, temp_storage=temp_storage
                    )
                    break
                except (HttpGatewayTimeout, ConnectTimeout):
                    continue
            pdb_paths.append(out_path)
        return pdb_paths
```

#### 2. Updated `AlphaFoldInterface`, `AF3Interface`, `ESMFoldInterface`, and `Null`

- **`evodesign/Prediction/AlphaFold/AlphaFoldInterface.py`**:
  `input_dir` and `output_dir` become private subdirectories inside a `with storage.temp_dir(...)` block, automatically cleaned up right after copying the final PDB to `output_pdb_path`:
  ```python
  def predict_single_pdb_file(
      self,
      sequence: str,
      output_pdb_path: str,
      temp_storage: TempStorage | None = None,
  ) -> str:
      storage = temp_storage or TempStorage()
      with storage.temp_dir(prefix="af_") as work_dir:
          input_dir = os.path.join(work_dir, "input")
          output_dir = os.path.join(work_dir, "output")
          os.makedirs(input_dir)
          os.makedirs(output_dir)

          protein_name = "protein_0"
          input_path = self._create_model_input(
              sequence, protein_name, input_dir, output_dir
          )
          self.run_inference(input_path, output_dir, do_batch_inference=False)
          prediction_pdb_path = self._prediction_pdb_path(protein_name, output_dir)

          os.makedirs(os.path.dirname(os.path.abspath(output_pdb_path)), exist_ok=True)
          shutil.copyfile(prediction_pdb_path, output_pdb_path)
      return output_pdb_path
  ```
- **`evodesign/Prediction/AlphaFold/AF3Interface.py` (batch `do`)**:
  Use a single `with storage.temp_dir(prefix="af3_batch_") as work_dir:` to stage all input JSONs in `work_dir/input`, run `self.run_inference(..., do_batch_inference=True)` into `work_dir/output`, copy each converted PDB into `output_dir/prediction_{i}.pdb`, and return the list of output PDB paths. All raw CIF/JSON/MSA artifacts are deleted automatically when the `with` block exits.
- **`evodesign/Prediction/ESMFold/ESMFoldInterface.py`**:
  Ignores `temp_storage` completely: `predict_single_pdb_str(sequence)` stays the abstract method, and `predict_single_pdb_file(sequence, output_pdb_path, temp_storage=None)` simply writes that string to `output_pdb_path` and returns `output_pdb_path`.
- **`evodesign/Prediction/Null.py`**:
  Writes directly to `output_pdb_path` and returns it (no `DirectoryManager` needed).

#### 3. Docker Volume Mount Note (`AlphaFold3Docker.py`)
`AlphaFold3Docker` bind-mounts `input_path` and `output_dir` via `--volume`.
- Inside an `evodesign` run, passing `TempStorage(paths.temp_dir)` ensures those temporary mount directories are created inside the job directory (which is already on a mountable host filesystem with sufficient disk space).
- In standalone mode, `TempStorage()` uses the system temp dir (respecting `$TMPDIR`, which can be overridden if `/tmp` is restricted or small).

#### 4. Standalone Usage Examples
```python
predictor = AlphaFold3Docker(model_dir="/path/to/models")

# Save directly to a chosen output path — no DirectoryManager needed
predictor.predict_single_pdb_file("MKWVTFISLL...", "my_result.pdb")

# Or skip managing output files altogether
pdb_text = predictor.predict_single_pdb_str("MKWVTFISLL...")
```

#### 5. Updated Optimization Loop (`CPD` & `SingleObjectiveCPD`)
Because structural evaluators like `evodesign/Metrics/Misc/RosettaEnergyFunction.py` read `context.get_model_chain().pdb_path` during fitness evaluation, predicted PDBs must stay on disk for the duration of `_compute_term_values(model_chains)`. Scoping a `temp_dir` around both prediction and evaluation handles this cleanly and deletes the generation's PDBs immediately afterward:

```python
# In SingleObjectiveCPD._evaluate:
if self.fitness_fn.requires_structure_predictor():
    with self.temp_storage.temp_dir(prefix="prediction_pdbs_") as pdbs_dir:
        pdb_paths = self.predictor.do(
            [sequences[k] for k in idx_not_archived],
            output_dir=pdbs_dir,
            temp_storage=self.temp_storage,
        )
        model_chains = [ChainFactory.create_from_pdb(p) for p in pdb_paths]
        new_terms = self._compute_term_values(model_chains)
else:
    model_chains = [
        ChainFactory.create_from_numpy(sequence_numpy)
        for sequence_numpy in x[idx_not_archived, :]
    ]
    new_terms = self._compute_term_values(model_chains)
```

In `_compute_term_values` (both `SingleObjectiveCPD` and `AlternatingCPD`), pass `temp_storage=self.temp_storage` into `self.fitness_fn.do(...)` so evaluators can access it through `Context`.

### Files Affected by `DirectoryManager` Removal
- **Delete**: `evodesign/Prediction/DirectoryManager.py`
- **Update**:
  - `evodesign/Prediction/Predictor.py`
  - `evodesign/Prediction/AlphaFold/AlphaFoldInterface.py`
  - `evodesign/Prediction/AlphaFold/AF3Interface.py`
  - `evodesign/Prediction/ESMFold/ESMFoldInterface.py`
  - `evodesign/Prediction/Null.py`
  - `evodesign/CPD.py`
  - `evodesign/SingleObjective/Problems/SingleObjectiveCPD.py`
  - `evodesign/SingleObjective/Problems/AlternatingCPD.py`
  - `evodesign/AlgorithmFactory.py`
  - `evodesign/SingleObjective/Algorithms/SingleObjectiveAlgorithm.py`
  - `evodesign/SingleObjective/Algorithms/GA/GenerationalPenalized.py`
  - `evodesign/SingleObjective/Algorithms/GA/GenerationalAlternatingFitness.py`
  - `evodesign/System/PathsContainer.py`
  - `evodesign/Callbacks/StorageManager.py`
  - `evodesign/__main__.py`

---

## 3. Refactor `Metrics/Context.py` & `ContextInterface.py`

### Problem
Using `Context` as a Blackboard / Evaluation Context to share intermediates across evaluators running in arbitrary order is a good pattern, and avoids creating a coupled `Manager` class. However, the current `get_extra_param_value` / `set_extra_param_value` implementation has four real issues:

1. **Silent parameter collisions and order-dependence (flat string keys)**:
   - `DescriptorSimilarity` defaults to `ESM2LocalLastLayer` with `submap_indices`, whereas `ContactMapSpectralDistance` defaults to `ESM2LocalAllLayers` with no `submap_indices`, yet both read/write `"esm2_model_desc_matrix"` and `"esm2_predicted_contacts"`. Whichever runs first silently dictates the tensor the other receives.
   - `"ref_contact_map"` is cached without `distance_threshold` or `min_separation` in both `ContactMapSpectralDistance` and `RCC`.
2. **Duplicated `if None: ... set(...)` boilerplate causing real bugs**:
   - `ContactMapPrecisionAtLength.py` (line 54): `context.set_extra_param_value("num_ref_contacts")` is missing the value argument.
   - `ContactMapSpectralDistance.py` (line 105), `ContactMapPrecisionAtLength.py` (line 65), and `RCC.py` (line 183): if another evaluator already cached `"esm2_predicted_contacts"`, the `if predicted_map is None:` block is skipped so `self.model` remains `None`, causing an `AttributeError` when calling `self.model.contacts_map_from_prediction(...)`.
   - `iLearnDescriptors.py` (lines 57–59): arguments to `context.set_extra_param_value(ref_descriptors, key)` are swapped relative to `(param_name, param_value)`.
3. **Lifetime mismatch (reference intermediates are never cached across candidates)**:
   - `FitnessFunction.do(model_chain, ref_chain)` creates a new `Context` for **every candidate sequence**.
   - As a result, reference intermediates (`esm2_ref_desc_matrix`, `ref_contact_map`, `ref_norm_laplacian`, `ref_rcc_vector`, `ref_distance_map`, `reference_ilearn_desc_*`) are recomputed for every candidate in every generation even though `ref_chain` is constant for the entire run.
4. **Mixing infrastructure with computation cache**:
   - System capabilities like `TempStorage` should be explicit, typed methods on `ContextInterface` rather than untyped keys in `extra_params`.

### Solution

#### 1. Updated `ContextInterface` & `Context`
Replace `get_extra_param_value` / `set_extra_param_value` with:
- An explicit typed accessor: `get_temp_storage() -> TempStorage`.
- Two generic memoization methods taking **hashable composite keys** and a **lazy factory**:
  - `get_or_compute(key, factory)`: memoized for the current candidate evaluation.
  - `get_or_compute_shared(key, factory)`: memoized across candidates for the same reference chain / run.

```python
from abc import ABC, abstractmethod
from typing import Any, Callable, Hashable, TypeVar

from ..Protein.Chain import Chain
from ..System.TempStorage import TempStorage

T = TypeVar("T")


class ContextInterface(ABC):
    @abstractmethod
    def get_metric_components(self, metric_name: str) -> dict[str, float]:
        raise NotImplementedError

    @abstractmethod
    def get_component_value(self, component_name: str) -> float:
        raise NotImplementedError

    @abstractmethod
    def get_model_chain(self) -> Chain:
        raise NotImplementedError

    @abstractmethod
    def get_reference_chain(self) -> Chain:
        raise NotImplementedError

    @abstractmethod
    def get_temp_storage(self) -> TempStorage:
        raise NotImplementedError

    @abstractmethod
    def get_or_compute(self, key: Hashable, factory: Callable[[], T]) -> T:
        raise NotImplementedError

    @abstractmethod
    def get_or_compute_shared(self, key: Hashable, factory: Callable[[], T]) -> T:
        raise NotImplementedError
```

```python
class Context(ContextInterface):
    def __init__(
        self,
        model_chain: Chain,
        ref_chain: Chain,
        metric_calculators: dict[str, Metric],
        temp_storage: TempStorage | None = None,
        shared_cache: dict[Hashable, Any] | None = None,
    ) -> None:
        self._model_chain = model_chain
        self._ref_chain = ref_chain
        self._metric_calculators = metric_calculators
        self._temp_storage = temp_storage or TempStorage()
        self._cached_results: dict[str, dict[str, float]] = {}
        self._candidate_cache: dict[Hashable, Any] = {}
        self._shared_cache: dict[Hashable, Any] = (
            shared_cache if shared_cache is not None else {}
        )

    def get_temp_storage(self) -> TempStorage:
        return self._temp_storage

    def get_or_compute(self, key: Hashable, factory: Callable[[], T]) -> T:
        """Memoized for the current candidate evaluation."""
        if key not in self._candidate_cache:
            self._candidate_cache[key] = factory()
        return self._candidate_cache[key]

    def get_or_compute_shared(self, key: Hashable, factory: Callable[[], T]) -> T:
        """Memoized across candidates for the same reference chain / run."""
        if key not in self._shared_cache:
            self._shared_cache[key] = factory()
        return self._shared_cache[key]
```

In `evodesign/Fitness/FitnessFunction.py`, hold a private `self._shared_cache: dict[Hashable, Any] = {}` (and `self._cached_ref_chain_id: int | None = None` so the cache resets automatically if `FitnessFunction.do(...)` is ever called with a different `ref_chain` instance) and pass both `temp_storage` and `shared_cache=self._shared_cache` into `Context(...)`. *(Because `_shared_cache` starts with `_`, `Settings.dump()` will ignore it.)*

#### 2. Shared Domain Cache Helpers (Purpose & Design)

**What a domain cache helper is and why it exists:**
- It is a tiny (~15-line) wrapper function that lives in the **domain module** that owns the computation (e.g., `evodesign/Metrics/ESM2/`), **not** in `Context`.
- Its job is to be the **single source of truth** for:
  1. The composite cache key (including all parameters that affect the output, preventing collisions).
  2. Choosing the right cache scope (`get_or_compute` for candidate vs. `get_or_compute_shared` for reference).
  3. Invoking the underlying computation inside `context.get_or_compute*`.
- Because the helper lives in the domain module, `Context` stays 100% generic and decoupled from ESM-2, contact maps, iLearn, etc.
- **When to write a helper**: only when **multiple evaluators share the same intermediate**. When an intermediate is private to a single evaluator (e.g., `DistanceMap` or `iLearnDescriptors`), call `context.get_or_compute_shared(...)` directly inside that evaluator.

In `evodesign`, only **two** shared helpers are needed:

##### Helper 1: ESM-2 query cache (`evodesign/Metrics/ESM2/cache.py`)
Shared by `ContactMapSpectralDistance`, `ContactMapPrecisionAtLength`, `ContactMapJensenShannon`, `DescriptorSimilarity`, and `RCC`:

```python
from typing import Literal

import numpy as np
import numpy.typing as npt

from ..ContextInterface import ContextInterface
from .ModelInterface import ModelInterface


def query_esm2_cached(
    context: ContextInterface,
    model: ModelInterface,
    target: Literal["model", "ref"],
    submap_indices: list[int] | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    submap_key = tuple(submap_indices) if submap_indices is not None else None
    key = ("esm2_query", type(model).__name__, target, submap_key)

    def _compute():
        chain = (
            context.get_model_chain()
            if target == "model"
            else context.get_reference_chain()
        )
        return model.query_model(chain.sequence, submap_indices=submap_indices)

    if target == "ref":
        return context.get_or_compute_shared(key, _compute)
    return context.get_or_compute(key, _compute)
```

##### Helper 2: Reference contact map cache (`evodesign/Metrics/cache_helpers.py`)
Shared by `ContactMapSpectralDistance`, `ContactMapPrecisionAtLength`, and `RCC`:

```python
import numpy as np
import numpy.typing as npt

from ..Protein.Geometry import compute_contacts_map
from .ContextInterface import ContextInterface


def get_ref_contact_map(
    context: ContextInterface,
    distance_threshold: float,
    min_separation: int,
) -> npt.NDArray[np.int64]:
    key = ("ref_contact_map", distance_threshold, min_separation)
    return context.get_or_compute_shared(
        key,
        lambda: compute_contacts_map(
            context.get_reference_chain().ca_atoms,
            distance_threshold,
            min_separation,
        ),
    )
```

*(Note: `SinglePairRepresentationSimilarity` can follow the same pattern if another evaluator ever shares ESMFold single/pair representations; otherwise it can call `context.get_or_compute` / `context.get_or_compute_shared` directly.)*

---

#### 3. How `do_for_fitness_fn` Changes Across Evaluators (Before & After)

##### Case A: Evaluators sharing intermediates via a helper

**`ContactMapSpectralDistance.do_for_fitness_fn` — Before:**
```python
def do_for_fitness_fn(self, context: ContextInterface) -> dict[str, float]:
    ref_laplacian = context.get_extra_param_value("ref_norm_laplacian")
    if ref_laplacian is None:
        ref_contact_map = context.get_extra_param_value("ref_contact_map")
        if ref_contact_map is None:
            ref_atoms = context.get_reference_chain().ca_atoms
            ref_contact_map = compute_contacts_map(
                ref_atoms, self.distance_threshold, self.min_separation
            )
            context.set_extra_param_value("ref_contact_map", ref_contact_map)
        ref_laplacian = self.normalized_laplacian(ref_contact_map)
        context.set_extra_param_value("ref_norm_laplacian", ref_laplacian)
    predicted_laplacian = context.get_extra_param_value("esm2_predicted_norm_laplacian")
    if predicted_laplacian is None:
        predicted_map = context.get_extra_param_value("esm2_predicted_contacts")
        if predicted_map is None:
            if self.model is None:
                from .ESM2LocalAllLayers import ESM2LocalAllLayers

                self.model = ESM2LocalAllLayers()
            model_sequence = context.get_model_chain().sequence
            model_desc_matrix, predicted_map = self.model.query_model(model_sequence)
            context.set_extra_param_value("esm2_model_desc_matrix", model_desc_matrix)
            context.set_extra_param_value("esm2_predicted_contacts", predicted_map)
        predicted_map = self.model.contacts_map_from_prediction(
            predicted_map, self.min_separation
        )
        predicted_laplacian = self.normalized_laplacian(predicted_map)
        context.set_extra_param_value(
            "esm2_predicted_norm_laplacian", predicted_laplacian
        )
    dist, norm_dist = self.do(predicted_laplacian, ref_laplacian)
    return {"distance": dist, "norm_distance": norm_dist}
```

**`ContactMapSpectralDistance.do_for_fitness_fn` — After:**
```python
def do_for_fitness_fn(self, context: ContextInterface) -> dict[str, float]:
    if self.model is None:
        from .ESM2LocalAllLayers import ESM2LocalAllLayers

        self.model = ESM2LocalAllLayers()

    ref_laplacian = context.get_or_compute_shared(
        ("ref_norm_laplacian", self.distance_threshold, self.min_separation),
        lambda: self.normalized_laplacian(
            get_ref_contact_map(context, self.distance_threshold, self.min_separation)
        ),
    )

    _, raw_contacts = query_esm2_cached(context, self.model, target="model")
    predicted_map = self.model.contacts_map_from_prediction(
        raw_contacts, self.min_separation
    )
    predicted_laplacian = self.normalized_laplacian(predicted_map)

    dist, norm_dist = self.do(predicted_laplacian, ref_laplacian)
    return {"distance": dist, "norm_distance": norm_dist}
```

**`DescriptorSimilarity.do_for_fitness_fn` — After:**
```python
def do_for_fitness_fn(self, context: ContextInterface) -> dict[str, float]:
    if self.model is None:
        from .ESM2LocalLastLayer import ESM2LocalLastLayer

        self.model = ESM2LocalLastLayer()

    ref_desc_matrix, _ = query_esm2_cached(
        context, self.model, target="ref", submap_indices=self.submap_indices
    )
    model_desc_matrix, _ = query_esm2_cached(
        context, self.model, target="model", submap_indices=self.submap_indices
    )
    rmse_val, norm_rmse, mean_norm_cos, distance, norm_distance = self.do(
        model_desc_matrix, ref_desc_matrix
    )
    return {
        "rmse": rmse_val,
        "norm_rmse": norm_rmse,
        "mean_norm_cos": mean_norm_cos,
        "distance": distance,
        "norm_distance": norm_distance,
    }
```

*(Apply the same pattern to `ContactMapPrecisionAtLength`, `ContactMapJensenShannon`, and `RCC`.)*

##### Case B: Evaluators with private reference caches (no helper needed)

**`DistanceMap.do_for_fitness_fn` — Before:**
```python
ref_distances = context.get_extra_param_value("ref_distance_map")
if ref_distances is None:
    ref_atoms = (
        context.get_reference_chain().ca_atoms
        if self.ca_atoms_only
        else context.get_reference_chain().backbone_atoms
    )
    ref_distances = self.compute_pairwise_distances(ref_atoms)
    context.set_extra_param_value("ref_distance_map", ref_distances)
```

**`DistanceMap.do_for_fitness_fn` — After:**
```python
ref_distances = context.get_or_compute_shared(
    ("ref_distance_map", self.ca_atoms_only),
    lambda: self.compute_pairwise_distances(
        context.get_reference_chain().ca_atoms
        if self.ca_atoms_only
        else context.get_reference_chain().backbone_atoms
    ),
)
```

**`iLearnDescriptors.do_for_fitness_fn` — After (combining `TempStorage` + `get_or_compute_shared`):**
```python
def do_for_fitness_fn(self, context: ContextInterface) -> dict[str, float]:
    temp_storage = context.get_temp_storage()

    def _compute_ref():
        with temp_storage.temp_file(suffix=".fasta") as fasta_path:
            return self.compute_descriptors_vector(
                fasta_path, context.get_reference_chain().sequence
            )

    ref_descriptors = context.get_or_compute_shared(
        ("ref_ilearn_desc", self.method), _compute_ref
    )

    with temp_storage.temp_file(suffix=".fasta") as fasta_path:
        model_descriptors = self.compute_descriptors_vector(
            fasta_path, context.get_model_chain().sequence
        )

    return {"rmse": self.do(model_descriptors, ref_descriptors)}
```

##### Case C: Evaluators with no shared/cached intermediates
No logic changes required for `pLDDT`, `RosettaEnergyFunction`, `HeadToTailAmide`, `HeadToTailDisulfide`, `Sum`, or `StructureSimilarity/{TMScore,RMSD,GDT,lDDT}`.

---

## 4. Additional Small Improvements & Bug Fixes

1. **`FitnessFunction.requires_structure_predictor()` (`evodesign/Fitness/FitnessFunction.py`)**:
   Currently checks every calculator in `self.term_calculators`, even if none of its outputs are listed in `self.terms`. Filter to only active calculators referenced in `self.terms`:
   ```python
   def requires_structure_predictor(self) -> bool:
       active_evaluators = {".".join(term.split(".")[:-1]) for term in self.terms}
       return any(
           self._term_calculators[name].requires_structure_predictor()
           for name in active_evaluators
       )
   ```
2. **`Metric.py` abstract method consistency**:
   `requires_structure_predictor` currently raises `RuntimeError`, whereas `do` and `do_for_fitness_fn` raise `NotImplementedError`. Standardize on `NotImplementedError`.
3. **Existing bugs to fix during this refactor**:
   - `ContactMapPrecisionAtLength.py`: missing value argument in `context.set_extra_param_value("num_ref_contacts")`.
   - `ContactMapSpectralDistance.py`, `ContactMapPrecisionAtLength.py`, `RCC.py`: `self.model` stays `None` if another evaluator populated the cache first.
   - `iLearnDescriptors.py`: swapped arguments in `context.set_extra_param_value(ref_descriptors, key)`.
   - `RCC.py`: `compute_contacts_map(ref_atoms, self.min_separation)` passes `min_separation` as `distance_threshold`.
   - `DistanceMap.py`: `self.compute_map(ref_atoms)` instead of `self.compute_pairwise_distances(ref_atoms)`.
   - `SinglePairRepresentationSimilarity.py` (line 35): `return rmse, norm_rmse` returns the imported `rmse` function instead of `rmse_val`.

---

## 5. Suggested Execution Order

To keep each commit small, testable, and easy to review:
1. **Step 1**: Add `evodesign/System/TempStorage.py`.
2. **Step 2**: Refactor `Prediction/` (`Predictor`, `AlphaFoldInterface`, `AF3Interface`, `ESMFoldInterface`, `Null`) to accept `output_pdb_path`/`output_dir` + optional `temp_storage`, and delete `DirectoryManager.py`.
3. **Step 3**: Update `PathsContainer`, `StorageManager`, `__main__.py`, `AlgorithmFactory`, `SingleObjectiveAlgorithm` subclasses, `CPD`, `SingleObjectiveCPD`, and `AlternatingCPD` to thread `TempStorage` instead of `DirectoryManager`.
4. **Step 4**: Refactor `ContextInterface`, `Context`, and `FitnessFunction` (`get_temp_storage`, `get_or_compute`, `get_or_compute_shared`, active-terms check in `requires_structure_predictor`).
5. **Step 5**: Add the two domain cache helpers (`query_esm2_cached`, `get_ref_contact_map`) and migrate the evaluators in `Metrics/ESM2/`, `Metrics/Misc/`, and `Metrics/Legacy/`.
