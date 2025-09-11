from evodesign.Algorithms.MonoObjective.GA.Generational import Generational
from evodesign.Utils.StorageManager import StorageManager
from evodesign.Utils.DirectoryManager import DirectoryManager
from evodesign.Utils.Chain import ChainFactory
from evodesign.Prediction.ESMFoldRemoteApi import ESMFoldRemoteApi
from evodesign.Fitness.WeightedMean import WeightedMean
from evodesign.Metrics.RMSD import RMSD
from evodesign.Metrics.GDT import GDT
from evodesign.GA.Mutation.RandomResetting import RandomResetting
from evodesign.Utils.Exceptions import *
from requests import ConnectTimeout
import numpy as np
import os


if __name__ == "__main__":
    example_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    target_pdb_path = os.path.join(example_folder, "1y32.pdb")
    max_generations = 2
    population_size = 5
    ref_chain = ChainFactory.create(target_pdb_path)
    directory = DirectoryManager(example_folder)
    ga = Generational(
        predictor=ESMFoldRemoteApi(),
        fitness_fn=WeightedMean(
            [RMSD(), GDT(cutoffs=[0.5, 1, 2, 4])],
            ["Metrics.RMSD.rmsd", "Metrics.GDT.gdt"],
            [0, 1],
        ),
        mutation=RandomResetting(sequence_mutation_prob=0.16667),
        max_generations=max_generations,
        population_size=population_size,
    )
    storage = StorageManager(
        directory,
        max_generations,
        population_size,
        len(ref_chain.sequence),
        ga.fitness_fn.num_terms(),
    )
    storage.save_git_commit_hash()
    storage.save_settings(ga.settings())
    storage.save_target_pdb(target_pdb_path)
    storage.save_rng_state(np.random.get_state(), directory.initial_rng_state_path)
    while True:
        try:
            ga.run(ref_chain, storage)
            storage.delete_non_essential_files_and_folders()
            break
        except (
            HttpInternalServerError,
            HttpGatewayTimeout,
            HttpForbidden,
            ConnectTimeout,
        ):
            storage.delete_non_essential_files_and_folders()
            continue
        except Exception as e:
            storage.delete_non_essential_files_and_folders()
            raise e
