from argparse import ArgumentParser
from pymoo.config import Config
import evodesign.Settings as Settings
from evodesign.Utils.DirectoryManager import DirectoryManager
from evodesign.Callbacks.StorageManager import StorageManager
from evodesign.Utils.ChainFactory import ChainFactory
from .Utils.Exceptions import *
from requests.exceptions import ConnectTimeout
import numpy as np
import json
import os

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
# parser.add_argument('-r', '--sequence_restrictions',
#                     type=str,
#                     default=None,
#                     help='path to the JSON file describing the allowed amino acids '
#                          'for certain positions in the designed sequences')
args = parser.parse_args()
with open(args.settings_path, "rt", encoding="utf-8") as json_file:
    settings = json.load(json_file)
algorithm = Settings.parse(settings)
ref_chain = ChainFactory.create_from_pdb(args.target_pdb_path, args.model_id, args.chain_id)
storage = StorageManager(
    DirectoryManager(os.path.abspath(args.output_dir), args.jobname),
    algorithm.max_generations,
    algorithm.population_size,
    len(ref_chain.sequence),
    algorithm.num_terms(),
    args.save_every_nth_generation,
)
try:
    # resuming from a previous execution
    algorithm._algorithm = storage.load_pymoo_algorithm()
    storage.load_results_npz()
    algorithm._algorithm.n_gen += 1
    if algorithm._algorithm.termination.n_max_gen < algorithm.max_generations:
        # extending from a previously completed execution
        storage.extend_result_arrays(algorithm.max_generations)
        algorithm._algorithm.termination.n_max_gen = algorithm.max_generations
        algorithm._algorithm.termination.perc = float(
            algorithm._algorithm.n_gen / algorithm.max_generations
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
    storage.save_git_commit_hash()
    storage.save_settings(algorithm.settings())
    storage.save_target_pdb(args.target_pdb_path)
while True:
    try:
        algorithm.run(ref_chain, storage)
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
