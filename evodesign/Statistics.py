from typing import Optional, Tuple

import numpy as np
import numpy.typing as npt
import pandas as pd
from matplotlib.axes import Axes

from .Chemistry.Sequences import NUM_AMINO_ACIDS


NUM_SERIES = 3




def get_final_solution_indices(
    generations: npt.NDArray[np.int64],
    fitness_values: npt.NDArray[np.float64],
) -> Tuple[int, int]:
    num_generations, population_size, sequence_length = generations.shape
    assert (
        num_generations == fitness_values.shape[0]
        and population_size == fitness_values.shape[1]
    )
    return np.unravel_index(np.argmax(fitness_values), fitness_values.shape)



def get_best_fitness_evolution(
    fitness_values: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    assert len(fitness_values.shape) == 2
    highest_per_generation = np.max(fitness_values, axis=1)
    return np.maximum.accumulate(highest_per_generation)



def get_population_missing_amino_acids(
    population: npt.NDArray[np.int64],
) -> float:
    sequence_length = population.shape[1]
    return np.mean(
        [
            NUM_AMINO_ACIDS - len(np.unique(population[:, col]))
            for col in range(sequence_length)
        ]
    )



def get_amino_acid_loss_evolution(
    generations: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    return np.array([get_population_missing_amino_acids(pop) for pop in generations])



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
    import seaborn as sns

    num_generations, population_size, sequence_length = generations.shape
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
        "Amino acid loss": data["Amino acid loss"] / NUM_AMINO_ACIDS,
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



def get_population_diversity_loss(
    population: npt.NDArray[np.int64],
    sample_size: Optional[int] = None,
) -> float:
    sequence_length = population.shape[1]
    l1 = get_population_missing_amino_acids(population) / NUM_AMINO_ACIDS
    l2 = get_population_identity(population, sample_size) / sequence_length
    diversity_loss = (l1 + l2) / 2
    return diversity_loss



def get_reason_for_termination(
    generations: npt.NDArray[np.int64], diversity_loss_tol: float = 1.0
) -> dict:
    num_generations, population_size, sequence_length = generations.shape
    unrun_mask = generations[:, 0, 0] == -1
    valid_indices = np.where(~unrun_mask)[0]
    if len(valid_indices) == 0:
        return {
            "status_code": -1,
            "status_message": "Error: No execution (Array entirely -1)",
            "last_generation": None,
            "final_diversity_loss": None,
        }
    last_gen_idx = valid_indices[-1]
    last_population = generations[last_gen_idx]
    final_diversity_loss = get_population_diversity_loss(last_population)
    if final_diversity_loss >= diversity_loss_tol:
        status = 1
        message = "Converged (Expected Stop): Diversity dropped below threshold."
    elif last_gen_idx == num_generations - 1:
        status = 0
        message = "Completed: Reached max generations without losing diversity."
    else:
        status = -2
        message = "Crashed (Unexpected Stop): Halted prematurely, but diversity was still healthy."
    return {
        "status_code": status,
        "status_message": message,
        "last_generation": last_gen_idx,
        "final_diversity_loss": final_diversity_loss,
    }
