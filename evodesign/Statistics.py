import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
import numpy as np
import numpy.typing as npt
from typing import Optional, Tuple


ALPHABET_SIZE = 20
NUM_SERIES = 3


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


def get_population_amino_acid_loss(
    population: npt.NDArray[np.int64],
) -> float:
    sequence_length = population.shape[1]
    return np.mean(
        [
            ALPHABET_SIZE - len(np.unique(population[:, col]))
            for col in range(sequence_length)
        ]
    )


def get_amino_acid_loss_evolution(
    generations: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    return np.array([get_population_amino_acid_loss(pop) for pop in generations])


def get_population_identity(
    population: npt.NDArray[np.int64],
    sample_size: Optional[int] = None,
) -> float:
    population_size = population.shape[0]
    sample = lambda: np.arange(population_size)
    if sample_size is not None and sample_size > 0:
        sample = lambda: np.random.choice(
            population_size, size=sample_size, replace=False
        )
    else:
        sample_size = population_size
    indices = sample()
    pop_sample = population[indices, :]
    identities = np.sum(pop_sample[:, None, :] == pop_sample[None, :, :], axis=2)
    sum_identity = np.triu(identities, 1).sum()
    num_pairs = (sample_size * (sample_size - 1)) // 2
    return sum_identity / num_pairs


def get_population_identity_evolution(
    generations: npt.NDArray[np.int64],
    sample_size: Optional[int] = None,
) -> npt.NDArray[np.float64]:
    return np.array([get_population_identity(pop, sample_size) for pop in generations])


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
        "Diversity loss": None,
        "New sequences": data["New sequences"] / population_size,
        "Amino acid loss": data["Amino acid loss"] / ALPHABET_SIZE,
        "Population identity": data["Population identity"] / sequence_length,
    }
    norm_data["Diversity loss"] = (
        norm_data["Amino acid loss"] + norm_data["Population identity"]
    ) / 2
    columns = list(norm_data.keys())
    sns_data = {
        "Generation": NUM_SERIES * norm_data["Generation"].tolist(),
        "Series values": np.concatenate(
            [norm_data[col] for col in columns[1:4]]
        ).tolist(),
        "Series": sum(([col] * num_generations for col in columns[1:4]), []),
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
