import os
import shutil
import argparse
import numpy as np


FILES_TO_COPY = [
    "initial_rng_state.txt",
    # "pymoo_algorithm.bin",
    "results.npz",
    "settings.json",
]


def copy_files(
    source_dir: str,
    dest_dir: str,
):
    if not os.path.exists(source_dir):
        raise RuntimeError
    os.makedirs(dest_dir, exist_ok=True)
    for filename in FILES_TO_COPY:
        source_file = os.path.join(source_dir, filename)
        dest_file = os.path.join(dest_dir, filename)
        if os.path.exists(source_file):
            # copy2 preserves metadata (timestamps, etc.)
            shutil.copy2(source_file, dest_file)
    return


def initialize_results_npz(destination_folder: str):
    npz_filename = os.path.join(destination_folder, "results.npz")
    data = np.load(npz_filename)
    generations = data["generations"]
    fitness_values = data["fitness_values"]
    term_values = data["term_values"]
    generations[1:] = -1
    fitness_values[1:] = 0.0
    term_values[1:] = -1.0
    fitness_values = np.zeros(fitness_values.shape, np.float64)
    term_values = np.full(term_values.shape, -1, np.float64)
    np.savez_compressed(
        npz_filename,
        generations=generations,
        fitness_values=fitness_values,
        term_values=term_values,
    )
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source_folder", type=str)
    parser.add_argument("destination_folder", type=str)
    args = parser.parse_args()
    copy_files(args.source_folder, args.destination_folder)
    initialize_results_npz(args.destination_folder)
