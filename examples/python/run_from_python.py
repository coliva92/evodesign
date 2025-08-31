from evodesign.Algorithms.MonoObjective.GA.Generational import Generational
from evodesign.Utils.SavingManager import SavingManager
from evodesign.Utils.WorkingFolder import WorkingFolder
from evodesign.Utils.Chain import ChainFactory
from evodesign.Prediction.ESMFoldRemoteApi import ESMFoldRemoteApi
from evodesign.Fitness.WeightedMean import WeightedMean
from evodesign.Metrics.RMSD import RMSD
import numpy as np
import os


if __name__ == "__main__":
    example_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    target_pdb_path = os.path.join(example_folder, "1y32.pdb")
    max_generations = 2
    population_size = 5
    ref_chain = ChainFactory.create(target_pdb_path)
    working_folder = WorkingFolder(example_folder)
    ga = Generational(
        predictor=ESMFoldRemoteApi(),
        fitness_fn=WeightedMean([RMSD()], ["Metrics.RMSD.rmsd"], [1]),
        max_generations=max_generations,
        population_size=population_size,
    )
    saving = SavingManager(
        working_folder,
        max_generations,
        population_size,
        len(ref_chain.sequence),
        ga.fitness_fn.num_terms(),
    )
    saving.save_git_commit_hash()
    saving.save_settings(ga.settings())
    saving.save_target_pdb(target_pdb_path)
    saving.save_rng_state(np.random.get_state(), working_folder.initial_rng_state_path)
    ga.run(ref_chain, saving)
    saving.delete_prediction_pdb()
