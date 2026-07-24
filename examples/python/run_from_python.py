from evodesign.Algorithms.MonoObjective.GA.Generational import Generational
from evodesign.Callbacks.StorageManager import StorageManager
from evodesign.System.PathsContainer import PathsContainer
from evodesign.Chemistry.ChainFactory import ChainFactory
from evodesign.Prediction.ESMFoldRemoteAPI import ESMFoldRemoteAPI
from evodesign.Fitness.WeightedMean import WeightedMean
from evodesign.Metrics.RMSD import RMSD
from evodesign.GA.Mutation.RandomResetting import RandomResetting
from pymoo.optimize import minimize
from evodesign.System.Exceptions import *
from requests import ConnectTimeout
import numpy as np
import os


if __name__ == "__main__":
    example_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    target_pdb_path = os.path.join(example_folder, "1y32.pdb")
    max_generations = 2
    population_size = 5
    ref_chain = ChainFactory.create_from_pdb(target_pdb_path)
    directory = PathsContainer.create(example_folder)
    ga = Generational(
        predictor=ESMFoldRemoteAPI(),
        fitness_fn=WeightedMean(
            [RMSD()],
            ["Metrics.RMSD.rmsd"],
            [1],
        ),
        mutation=RandomResetting(sequence_mutation_prob=0.1),
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
            algorithm = ga.create_algorithm()
            problem = ga.create_problem(
                ref_chain, 
                storage.predictor_directory)
            callbacks = ga.create_callbacks(storage)
            minimize(problem=problem, 
                     algorithm=algorithm, 
                     callback=callbacks,
                     verbose=True,)
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
