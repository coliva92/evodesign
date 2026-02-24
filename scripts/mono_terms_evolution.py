from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
from typing import List
from pathlib import Path
from evodesign.Settings import read_json
import os


def get_term_names(settings: dict) -> List[str]:
    algorithm = list(settings.keys())[0]
    fitness_fn = list(settings[algorithm]["fitness_fn"].keys())[0]
    return settings[algorithm]["fitness_fn"][fitness_fn]["terms"]


def terms_evolution_from_folder(
    indices: List[int],
    input_folder_path: str,
    output_folder_path: str,
):
    root_dir = Path(input_folder_path)
    for i, run_dir in enumerate(root_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        npz_path = run_dir / "results.npz"
        json_path = run_dir / "settings.json"
        if not npz_path.exists() or not json_path.exists():
            continue
        settings = read_json(json_path)
        tmp_names = get_term_names(settings)
        term_names = [tmp_names[j] for j in indices]
        prefix = os.path.join(output_folder_path, run_dir.name)
        plot_terms_evolution(indices, term_names, f"{prefix}_terms_evolution.png")


def plot_terms_evolution(
    indices: List[int],
    labels: List[str],
    output_png_path: str,
):
    assert len(indices) == len(labels)
    with np.load("results.npz") as data:
        term_values = data["term_values"]
        generations = data["generations"]
        fitness_values = data["fitness_values"]
        best = np.zeros((generations.shape[0], len(indices)))
        last_fitness = 0
        for i in range(generations.shape[0]):
            j = np.argmax(fitness_values[i, :])
            if fitness_values[i, j] > last_fitness:
                best[i, :] = term_values[i, j, indices]
                last_fitness = fitness_values[i, j]
            else:
                best[i, :] = best[i - 1, :]
    plt.figure(figsize=(10, 6))
    plt.plot(best)
    plt.savefig(output_png_path)
    plt.clf()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input_folder_path", type=str)
    parser.add_argument("output_folder_path", type=str)
    parser.add_argument("term_indices", nargs="+", type=int)
    args = parser.parse_args()
    terms_evolution_from_folder(
        args.term_indices, args.input_folder, args.output_folder
    )
