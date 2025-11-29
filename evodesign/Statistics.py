from Utils.ChainFactory import ChainFactory, AMINO_ACIDS
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import numpy as np
import numpy.typing as npt
from typing import Optional, Tuple, List


ALPHABET_SIZE = 20
NUM_SERIES = 4


def best_sequence_indices(
    generations: npt.NDArray[np.int64],
    fitness_values: npt.NDArray[np.float64],
) -> Tuple[int, int]:
    (num_generations, population_size, sequence_length) = generations.shape
    assert (
        num_generations == fitness_values.shape[0]
        and population_size == fitness_values.shape[1]
    )
    highest_per_generation = np.max(fitness_values, axis=1)
    i = np.argmax(highest_per_generation, axis=0)
    j = np.argmax(fitness_values[i, :], axis=0)
    return i, j


def best_fitness_per_generation(
    fitness_values: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    assert len(fitness_values.shape) == 2
    highest_per_generation = np.max(fitness_values, axis=1)
    return np.maximum.accumulate(highest_per_generation)


def amino_acid_loss_per_generation(
    generations: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    sequence_length = generations.shape[2]
    unique_counts = np.array(
        [
            np.array(
                [
                    ALPHABET_SIZE - len(np.unique(pop[:, col]))
                    for col in range(sequence_length)
                ]
            )
            for pop in generations
        ]
    )
    avg_losses = np.mean(unique_counts, axis=1)
    return avg_losses


def population_identity_per_generation(
    generations: npt.NDArray[np.int64],
    sample_size: Optional[int] = None,
) -> npt.NDArray[np.float64]:
    num_generations, population_size = generations.shape[0], generations.shape[1]
    avg_identity = np.zeros(num_generations, dtype=np.float64)
    if sample_size is None:
        sample = lambda pop: pop
        sample_size = population_size
    for i in range(num_generations):
        pop = sample(generations[i])
        identities = np.sum(pop[:, None, :] == pop[None, :, :], axis=2)
        sum_identity = np.triu(identities, 1).sum()
        num_pairs = (sample_size * (sample_size - 1)) // 2
        avg_identity[i] = sum_identity / num_pairs
    return avg_identity


def num_new_sequences_per_generation(
    generations: npt.NDArray[np.int64],
) -> npt.NDArray[np.int64]:
    num_generations = generations.shape[0]
    counts = np.zeros(num_generations, dtype=np.int64)
    unique_solutions = set()
    for i in range(num_generations):
        # convert sequences to hashable tuples
        sequences = map(tuple, generations[i])
        num_new_solutions = 0
        for seq in sequences:
            if seq not in unique_solutions:
                num_new_solutions += 1
                unique_solutions.add(seq)
        counts[i] = num_new_solutions
    return counts


def convergence_plot(
    generations: npt.NDArray[np.int64],
    fitness_values: npt.NDArray[np.float64],
    color_palette_name: str = "colorblind",
) -> Tuple[Axes, pd.DataFrame]:
    (num_generations, population_size, sequence_length) = generations.shape
    assert (
        num_generations == fitness_values.shape[0]
        and population_size == fitness_values.shape[1]
    )
    data = {
        "Generation": np.arange(1, generations.shape[0] + 1),
        "Best fitness": best_fitness_per_generation(fitness_values),
        "Amino acid loss": amino_acid_loss_per_generation(generations),
        "Population identity": population_identity_per_generation(generations),
        "New sequences": num_new_sequences_per_generation(generations),
    }
    norm_data = {
        "Generation": data["Generation"],
        "Best fitness": data["Best fitness"],
        "Amino acid loss": data["Amino acid loss"] / ALPHABET_SIZE,
        "Population identity": data["Population identity"] / sequence_length,
        "New sequences": data["New sequences"] / population_size,
    }
    columns = list(data.keys())
    sns_data = {
        "Generation": NUM_SERIES * norm_data["Generation"].tolist(),
        "Series values": np.concatenate(
            [norm_data[col] for col in columns[1:]]
        ).tolist(),
        "Series": sum(([col] * num_generations for col in columns[1:]), []),
    }
    df = pd.DataFrame.from_dict(sns_data)
    ax = sns.lineplot(
        data=df,
        x="Generation",
        y="Series values",
        hue="Series",
        palette=sns.color_palette(color_palette_name),
    )
    return ax, pd.DataFrame.from_dict(data)


def mono_term_names(settings_json_path: str) -> List[str]:
    import json

    with open(settings_json_path, "rt", encoding="utf-8") as f:
        settings = json.load(f)
    algorithm = list(settings.keys())[0]
    fitness_fn = list(settings[algorithm]["fitness_fn"].keys())[0]
    return ["Fitness"] + settings[algorithm]["fitness_fn"][fitness_fn]["terms"]


def mono_statistics_from_folder(
    folder_dir: str,
    output_dir: str,
) -> None:
    from pathlib import Path
    import os

    df = None
    root_dir = Path(folder_dir)
    for run_dir in root_dir.iterdir():
        if not run_dir.is_dir():
            continue
        npz_path = run_dir / "results.npz"
        json_path = run_dir / "settings.json"
        pdb_path = list(root_dir.glob("*.pdb"))[0]
        if not npz_path.exists() or not json_path.exists() or not pdb_path.exists():
            continue
        data = np.load(npz_path)
        generations = data["generations"]
        fitness_values = data["fitness_values"]
        ax, _ = convergence_plot(generations, fitness_values)
        filename_prefix = os.path.join(output_dir, run_dir.name)
        fig = ax.get_figure()
        fig.savefig(f"{filename_prefix}.svg", format="svg")
        fig.clf()
        plt.clf()
        term_names = mono_term_names(json_path)
        term_values = data["term_values"]
        i, j = best_sequence_indices(generations, fitness_values)
        sequence = "".join(AMINO_ACIDS[k] for k in generations[i, j, :])
        ref_chain = ChainFactory.create_from_pdb(pdb_path)
        identity = np.average(generations[i, j, :] == ref_chain.sequence_numpy)
        base = {
            "ID": run_dir.name,
            "Sequence": sequence,
            "ReferenceSequence": ref_chain.sequence,
            "SequenceLength": len(sequence),
            "Identity": identity,
        }
        terms = {k: v for k, v in zip(term_names, term_values[i, j, :])}
        terms["Fitness"] = fitness_values[i, j]
        row = pd.DataFrame([{**base, **terms}])
        if df is not None:
            df = pd.concat([df, row], ignore_index=True)
        else:
            df = row
    filename = os.path.join(output_dir, f"{root_dir.name}_solutions.csv")
    df.to_csv(filename, index=False)
