from .WorkingFolder import WorkingFolder
from pymoo.core.callback import Callback
from pymoo.core.algorithm import Algorithm as PyMOOAlgorithm
import numpy as np
import dill
import json
import os
import shutil


class SavingManager(Callback):

    def __init__(
        self,
        working_folder: WorkingFolder,
        num_generations: int,
        population_size: int,
        sequence_length: int,
        num_fitness_fn_terms: int,
        save_every_nth_generation: int = 10,
    ):
        super().__init__()
        self.working_folder = working_folder
        shape_3d = (num_generations, population_size, sequence_length)
        shape_2d = (num_generations, population_size)
        self.generations = np.full(shape_3d, -1, np.int64)
        self.fitness_values = np.zeros(shape_2d, np.float64)
        self.term_values = np.full(
            (num_generations, population_size, num_fitness_fn_terms), -1, np.float64
        )
        self.save_every_nth_generation = save_every_nth_generation

    def _extend_numpy_array(self, arr, new_size, fill_value=0):
        old_shape = arr.shape
        if old_shape[0] == new_size:
            return arr  # No resizing needed
        new_shape = (new_size,) + old_shape[1:]
        new_arr = np.full(new_shape, fill_value, dtype=arr.dtype)
        copy_len = min(old_shape[0], new_size)
        new_arr[:copy_len] = arr[:copy_len]
        return new_arr

    def extend_result_arrays(self, new_size: int) -> None:
        self.generations = self._extend_numpy_array(self.generations, new_size, -1)
        self.fitness_values = self._extend_numpy_array(self.fitness_values, new_size, 0)
        self.term_values = self._extend_numpy_array(self.term_values, new_size, -1)

    def notify(self, algorithm: PyMOOAlgorithm) -> None:
        generation_idx = algorithm.n_iter - 1
        for i, solution in enumerate(algorithm.pop.get("X")):
            for j, x in enumerate(solution):
                self.generations[generation_idx][i][j] = x
        for i, y in enumerate(algorithm.pop.get("F")):
            self.fitness_values[generation_idx][i] = -1.0 * y[0]
        for i, w in enumerate(algorithm.problem.term_values):
            self.term_values[generation_idx][i] = w
        if (
            generation_idx % self.save_every_nth_generation == 0
            or algorithm.termination.perc == 1.0
        ):
            self.save(algorithm)

    def save(self, algorithm: PyMOOAlgorithm) -> None:
        # save the results for offline analysis
        np.savez_compressed(
            self.working_folder.results_npz_path,
            generations=self.generations,
            fitness_values=self.fitness_values,
            term_values=self.term_values,
        )
        # save RNG for reproduction/resuming
        self.save_rng_state(
            np.random.get_state(), self.working_folder.last_rng_state_path
        )
        # save the settings for reproduction/resuming
        # save the algorithm for resuming later
        self.save_pymoo_algorithm(algorithm)

    def save_rng_state(self, state: tuple, file_path: str) -> None:
        os.makedirs(self.working_folder.path, exist_ok=True)
        with open(file_path, "wt", encoding="utf-8") as txt_file:
            txt_file.write(f"{state[0]}\n")
            values = ",".join([str(v) for v in state[1]])
            txt_file.write(f"{values}\n")
            for i in range(2, len(state)):
                txt_file.write(f"{state[i]}\n")

    def save_pymoo_algorithm(self, algorithm: PyMOOAlgorithm) -> None:
        os.makedirs(self.working_folder.path, exist_ok=True)
        file_path = self.working_folder.pymoo_algorithm_bin_path
        with open(file_path, "wb") as bin_file:
            tmp_problem = algorithm.problem
            tmp_callback = algorithm.callback
            algorithm.problem = algorithm.callback = None
            dill.dump(algorithm, bin_file)
            algorithm.problem = tmp_problem
            algorithm.callback = tmp_callback

    def save_git_commit_hash(self) -> None:
        os.makedirs(self.working_folder.path, exist_ok=True)
        file_path = self.working_folder.git_commit_hash_path
        with open(file_path, "wt", encoding="utf-8") as txt_file:
            commit_hash = self.working_folder.GIT_COMMIT_HASH
            txt_file.write(
                f"https://github.com/coliva92/evodesign/commit/{commit_hash}\n"
            )

    def save_settings(self, settings: dict) -> None:
        os.makedirs(self.working_folder.path, exist_ok=True)
        file_path = self.working_folder.settings_json_path
        with open(file_path, "wt", encoding="utf-8") as json_file:
            json_file.write(json.dumps(settings, indent=4) + "\n")

    def save_target_pdb(self, target_pdb_path: str) -> None:
        pdb_filename = os.path.basename(target_pdb_path)
        destination = os.path.join(self.working_folder.path, pdb_filename)
        shutil.copy(target_pdb_path, destination)

    def load_results_npz(self) -> None:
        data = np.load(self.working_folder.results_npz_path)
        self.generations = data["generations"]
        self.fitness_values = data["fitness_values"]
        self.term_values = data["term_values"]

    def load_rng_state(self, file_path: str) -> tuple:
        result = []
        i = 0
        for line in open(file_path, "rt", encoding="utf-8"):
            if i == 0:
                result.append(line.strip())
            elif i == 1:
                values = [int(float(x)) for x in line.strip().split(",")]
                np_values = np.zeros(len(values), dtype=np.uint32)
                for j in range(len(values)):
                    np_values[j] = values[j]
                result.append(np_values)
            elif i == 4:
                result.append(float(line.strip()))
            else:
                result.append(int(float(line.strip())))
            i += 1
        return tuple(result)

    def load_pymoo_algorithm(self) -> PyMOOAlgorithm:
        file_path = self.working_folder.pymoo_algorithm_bin_path
        with open(file_path, "rb") as bin_file:
            algorithm = dill.load(bin_file)
        return algorithm

    def delete_prediction_pdb(self) -> None:
        if os.path.isfile(self.working_folder.prediction_pdb_path):
            os.remove(self.working_folder.prediction_pdb_path)
