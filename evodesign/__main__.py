from argparse import ArgumentParser

import numpy as np
from pymoo.config import Config
from pymoo.optimize import minimize
from requests.exceptions import ConnectTimeout

import evodesign.Settings as Settings
from evodesign.Callbacks.StorageManager import StorageManager
from evodesign.Chemistry.ChainFactory import ChainFactory
from evodesign.Chemistry.Sequences import load_profile, to_str
from evodesign.System.Exceptions import *
from evodesign.System.PathsContainer import PathsContainer

Config.warnings["not_compiled"] = False
parser = ArgumentParser(
    prog="evodesign",
    description="A rudimentary suite of evolutionary " "algorithms for protein design.",
)
parser.add_argument(
    "target_pdb_path", help="path to the PDB file of the target protein backbone"
)
parser.add_argument(
    "settings_path",
    help="path to the JSON file describing the configuration "
    "of the evolutionary algorithm to be executed",
)
parser.add_argument("output_dir", help="path to the output folder")
parser.add_argument(
    "-s",
    "--save_every_nth_generation",
    type=int,
    default=10,
    help="save every time this number of generations have passed",
)
parser.add_argument(
    "-m",
    "--model_id",
    type=int,
    default=0,
    help="the zero-based index of the model to read from the target PDB file",
)
parser.add_argument(
    "-c",
    "--chain_id",
    type=str,
    default=None,
    help="the chain ID to read from the target PDB file",
)
parser.add_argument(
    "-j", "--jobname", type=str, default=None, help="the name for the current job"
)
parser.add_argument(
    "-p",
    "--profile_path",
    type=str,
    default=None,
    help="path to the text file describing the probability of "
    "allowed each amino acids on each positions in the "
    "designed sequences",
)
args = parser.parse_args()
algorithm_factory = Settings.load(args.settings_path)
ref_chain = ChainFactory.create_from_pdb(
    args.target_pdb_path, args.model_id, args.chain_id
)
storage = StorageManager(
    PathsContainer.create(args.output_dir, args.jobname),
    algorithm_factory.max_generations,
    algorithm_factory.population_size,
    len(ref_chain.sequence),
    algorithm_factory.num_terms(),  # TODO: we're assuming mono-objective only
    args.save_every_nth_generation,
)
algorithm = algorithm_factory.create_algorithm()
try:
    # resuming from a previous execution
    storage.load_results_npz()
    algorithm = storage.load_pymoo_algorithm()
    algorithm.n_gen += 1  # TODO: check if this is correct
    if algorithm.termination.n_max_gen < algorithm_factory.max_generations:
        # extending from a previously completed execution
        storage.extend_result_arrays(algorithm_factory.max_generations)
        algorithm.termination.n_max_gen = algorithm.max_generations
        algorithm.termination.perc = float(
            algorithm.n_gen / algorithm_factory.max_generations
        )
    state = storage.load_rng_state(storage.directory.last_rng_state_path)
    np.random.set_state(state)
except FileNotFoundError:
    try:
        # starting fresh but with a previous RNG seed
        state = storage.load_rng_state(storage.directory.initial_rng_state_path)
        np.random.set_state(state)
    except FileNotFoundError:
        # starting with a fresh RNG seed
        storage.save_rng_state(
            np.random.get_state(), storage.directory.initial_rng_state_path
        )
    storage.save_version()
    storage.save_settings(algorithm_factory.settings())
    storage.save_target_pdb(args.target_pdb_path)
while True:
    try:
        aa_profile = None
        if args.profile_path is not None:
            aa_profile = load_profile(args.profile_path)
        problem = algorithm_factory.create_problem(
            ref_chain, storage.predictor_directory, aa_profile
        )
        callbacks = algorithm_factory.create_callbacks(storage)
        result = minimize(
            problem=problem,
            algorithm=algorithm,
            callback=callbacks,
            verbose=True,
            copy_algorithm=False,
        )
        solution = to_str(result.X)
        fitness_value = result.F[0]
        print("Sequence:", solution)
        print("Fitness:", fitness_value)
    except (HttpInternalServerError, HttpGatewayTimeout, HttpForbidden, ConnectTimeout):
        # Cleanup, then retry
        storage.delete_non_essential_files_and_folders()
        continue
    except Exception as e:
        # Cleanup, then propagate
        storage.delete_non_essential_files_and_folders()
        raise e
    else:
        # Normal completion: cleanup, then exit loop
        storage.delete_non_essential_files_and_folders()
        break
