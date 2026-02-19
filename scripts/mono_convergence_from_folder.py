from argparse import ArgumentParser
from evodesign.Statistics import create_convergence_plot, get_final_solution_indices
from evodesign.Settings import parse
from evodesign.Chemistry.ChainFactory import ChainFactory
from evodesign.Chemistry.Sequences import AMINO_ACIDS
from evodesign.Settings import read_json
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from typing import List


def get_fitness_fn_name(settings: dict) -> str:
    algorithm = list(settings.keys())[0]
    f = parse(settings[algorithm]["fitness_fn"])
    return f.name()


def get_term_names(settings: dict) -> List[str]:
    algorithm = list(settings.keys())[0]
    fitness_fn = list(settings[algorithm]["fitness_fn"].keys())[0]
    return settings[algorithm]["fitness_fn"][fitness_fn]["terms"]


def get_predictor_name(settings: dict) -> str:
    algorithm = list(settings.keys())[0]
    predictor_name = list(settings[algorithm]["predictor"].keys())[0]
    return predictor_name.split(".")[-1]


def compute_statistics_from_folder(
    folder_dir: str,
    output_dir: str,
) -> None:
    from pathlib import Path
    import os

    df = None
    root_dir = Path(folder_dir)
    for i, run_dir in enumerate(root_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        npz_path = run_dir / "results.npz"
        json_path = run_dir / "settings.json"
        pdb_path = list(run_dir.glob("*.pdb"))[0]
        if not npz_path.exists() or not json_path.exists() or not pdb_path.exists():
            continue
        data = np.load(npz_path)
        generations = data["generations"]
        fitness_values = data["fitness_values"]
        ax, _ = create_convergence_plot(generations, fitness_values)
        filename_prefix = os.path.join(output_dir, run_dir.name)
        fig = ax.get_figure()
        fig.savefig(f"{filename_prefix}.svg", format="svg")
        fig.clf()
        plt.clf()
        settings = read_json(json_path)
        term_names = get_term_names(settings)
        term_values = data["term_values"]
        i, j = get_final_solution_indices(generations, fitness_values)
        sequence = "".join(AMINO_ACIDS[k] for k in generations[i, j, :])
        ref_chain = ChainFactory.create_from_pdb(pdb_path)
        identity = np.average(generations[i, j, :] == ref_chain.sequence_numpy)
        base = {
            "ID": f"{root_dir.stem}_{i}",
            "Designer": root_dir.stem,
            "Predictor": get_predictor_name(settings),
            "FitnessFunction": get_fitness_fn_name(settings),
            "ReferencePDBName": pdb_path.stem,
            "ReferencePDBPath": pdb_path.resolve(),
            "ReferenceSequence": ref_chain.sequence,
            "Sequence": sequence,
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
    filename = os.path.join(output_dir, f"{root_dir.stem}_terms.csv")
    df.to_csv(filename, index=False)


if __name__ == "__main__":
    parser = ArgumentParser("mono_convergence_from_folder")
    parser.add_argument("folder_dir", type=str)
    parser.add_argument("output_dir", type=str)
    args = parser.parse_args()
    compute_statistics_from_folder(args.folder_dir, args.output_dir)
