from argparse import ArgumentParser
from evodesign.Settings import read_json, parse
from evodesign.Chemistry.ChainFactory import ChainFactory
import numpy as np
import os
from typing import Tuple


def get_params_from_settings(
    settings: dict,
) -> Tuple[int, int, int]:
    algorithm = list(settings.keys())[0]
    num_generations = settings[algorithm]["num_generations"]
    population_size = settings[algorithm]["population_size"]
    fitness_fn = list(settings[algorithm]["fitness_fn"].keys())[0]
    F = parse(fitness_fn)
    num_fitness_fn_terms = F.num_terms()
    return num_generations, population_size, num_fitness_fn_terms


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("fasta_path", type=str)
    parser.add_argument("output_dir", type=str)
    parser.add_argument("settings_json_path", type=str)
    args = parser.parse_args()

    sequences = []
    curr_seq = None
    for line in open(args.fasta_path, "rt", encoding="utf-8"):
        if line[0] == ">":
            if len(curr_seq) > 0:
                sequences.append(curr_seq)
            if len(sequences) >= 2:
                assert len(sequences[-2]) == len(sequences[-1])
            curr_seq = ""
            continue
        curr_seq += line.strip()

    settings = read_json(args.settings_json_path)
    num_generations, population_size, num_fitness_fn_terms = get_params_from_settings(settings)
    assert(population_size == len(sequences))
    sequence_length = len(sequences[0])
    shape_3d = (num_generations, population_size, sequence_length)
    shape_2d = (num_generations, population_size)
    generations = np.full(shape_3d, -1, np.int64)
    fitness_values = np.zeros(shape_2d, np.float64)
    term_values = np.full(
        (num_generations, population_size, num_fitness_fn_terms), -1, np.float64
    )

    for i, seq in enumerate(sequences):
        generations[0, i, :] = ChainFactory.sequence_str_to_numpy(seq)

    np.savez_compressed(
        os.path.join(args.output_dir, "results.npz"),
        generations=generations,
        fitness_values=fitness_values,
        term_values=term_values,
    )
