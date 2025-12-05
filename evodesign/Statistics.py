
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
import numpy as np
import numpy.typing as npt
from typing import Optional, Tuple


ALPHABET_SIZE = 20
NUM_SERIES = 4


def get_final_solution_indices(
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


def get_best_fitness_evolution(
    fitness_values: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    assert len(fitness_values.shape) == 2
    highest_per_generation = np.max(fitness_values, axis=1)
    return np.maximum.accumulate(highest_per_generation)


def get_amino_acid_loss_evolution(
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


def get_population_identity_evolution(
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


def get_new_sequences_evolution(
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


def create_convergence_plot(
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
        "Best fitness": get_best_fitness_evolution(fitness_values),
        "Amino acid loss": get_amino_acid_loss_evolution(generations),
        "Population identity": get_population_identity_evolution(generations),
        "New sequences": get_new_sequences_evolution(generations),
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
