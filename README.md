# EvoDesign

This is a Python framework for implementing and testing different evolutionary algorithms for protein design.
EvoDesign can easily integrate AI and non-AI models for various tasks (e.g., AlphaFold and ESMFold for structure prediction, ESM-2 for protein descriptors, etc.). 
EvoDesign is built on top of the [PyMOO library](https://pymoo.org/).

## Installing with Anaconda or Miniconda

The following instructions assume you are using a Linux system.

1. Clone this repository into your system.
2. Navigate to the directory of the cloned repository; e.g.: `cd evodesign/`.
3. Create a new conda environment from the provided `environment.yml` file:
```
conda env create -f environment.yml
```
4. When running EvoDesign from the terminal, don't forget to activate the installed environment:
```
conda activate evodesign
```

## Installing with Docker

If you prefer, you can create a Docker image with EvoDesign pre-installed.
It is highly recommended you install [Rootless Docker](https://docs.docker.com/engine/security/rootless/#distribution-specific-hint) *before* installing EvoDesign's docker image.

1. Clone this repository into your system.
2. Navigate to the directory of the cloned repository; e.g.: `cd evodesign/`.
3. Create the docker image:
```
docker build -t evodesign:latest .
```

## Running the examples

Regardless of the installation method used (Anaconda or Docker), you can verify the installation was successfull by running the examples provided in the `examples/` folder. 
These examples run a simple genetic algorithm on a dummy PDB file.
We recommend you check the source code of these examples if you want to get a quick instruction on how to run EvoDesign.

To run the `python` example, navigate that folder (e.g., using `cd examples/python`) and use the following command:
```
python run_from_python.py
```

To run the `terminal` example, navigate to that folder (e.g., using `cd examples/terminal`) and use the following command:
```
./run_from_terminal.sh
```

To run the `docker` example, navigate to that folder (e.g., using `cd examples/docker`) and use the following command:
```
./run_from_docker.sh
```

## Installation of the external models

The EvoDesign framework provides wrapping code for using external models for various tasks (e.g., AlphaFold for structure prediction, ESM2 for protein descriptors, etc.).
These models need to be installed separately.
We recommend installing only the models that will be used in your design tasks.
Bellow are links to each model's official installation instructions:

- [AlphaFold 2](https://github.com/google-deepmind/alphafold).
- [AlphaFold 3](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md)
- [ESM2 and ESMFold](https://github.com/facebookresearch/esm).
- [PyRosetta](https://www.pyrosetta.org/downloads#h.iwt5ktel05jc).
- [iLearn](https://github.com/Superzchen/iLearn).

### Tips for installing these models

**Create a separate Anaconda environment for each model** to avoid dependency conflicts. 
After a model is installed in its own environment, we recommend installing EvoDesign in that same environment. 
This can be done using the `pip` provided in this repository.
For example, say you installed AlphaFold 3 in an environment called `alphafold3`:

1. Clone this repository if you haven't already.
2. Navigate to EvoDesign's repository, e.g., using `cd ~/evodesign/`.
3. Activate the model's environment:
```
conda activate alphafold3
```
4. Install EvoDesign:
```
pip install -e .
```

**To install the iLearn descriptors, simply clone the official GitHub repository.** After that, modify your `PYTHONPATH` environment variable to add the path to the repository.
For example, say you cloned the iLearn repository to your home folder, use these commands to modify the `PYTHONPATH`:

```
export PYTHONPATH=${PYTHONPATH}:~/Documents/iLearn
source ~/.bashrc
```

You can also edit the `~/.bashrc` file to make these changes permanent.

## Running EvoDesign from the Linux terminal

This is the most common method for running EvoDesign. 
The `examples/terminal` folder contain an example code on how to use EvoDesign in this manner.
Before you run EvoDesign, don't forget to activate the conda environment were you installed it; e.g., using: `conda activate evodesign`.
You can use this command to review EvoDesign's arguments:
```
python3 -m evodesign --help
```

Follow these steps to run a custom evolutionary algorithm using EvoDesign:
1. Write a JSON configuration file describing the characteristics of the algorithm to be executed (check the `examples/terminal/settings.json` for an example). 
2. Gather the PDB file of the target protein. 
3. Run the following command: 
```
python3 -m evodesign path/to/target_protein.pdb \
                     path/to/settings.json \
                     path/to/output_folder \
                     -j custom_jobname 
```

The `custom_jobname` can be any name you want to use to identify the results from this execution.
After finishing execution, an folder named with the especified jobname will be created in the especified output folder.
The created folder will contain the following structure:

```
.
├─ custom_jobname/
|  ├─ initial_rng_state.txt
|  ├─ last_rng_state.txt
|  ├─ pymoo_algorithm.bin
|  ├─ results.npz
|  ├─ settings.json
|  ├─ target_protein.pdb
|  ├─ version.txt
```

These files are stored mainly for reproduction purposes. 
- The `version.txt` file contains the specific EvoDesign version used to produced these results.
- The `initial_rng_state.txt` and `last_rng_state.txt` files contain the RNG state at the beginning and at the end of the last iteration of the algorithm's execution, respectively.
- The `pymoo_algorithm.bin` file stores the state of the PyMOO evolutionary algorithm in the last iteration.
- The `results.npz` file stores the sequences and fitness values produced by the evolutionary algorithm in each iteration. See below for more details.
- The `settings.json` file is a copy of the settings file that was provided as input to the `python3 -m evodesign` command.
- The `*.pdb` file is a copy of the target protein PDB file that was provided as input to the `python3 -m evodesign` command. In this example, this file is named `target_protein.pdb`.

By default, these files are updated every 10 iterations.

### The `results.npz` file

For classic single-objective genetic algorithms, the `results.npz` file stores three NumPy arrays called `generations`, `fitness_values` and `term_values`.

- The `generations` array store the amino acid sequences (represented as integers in the \[0, 19\] interval) produced by the algorithm on each iteration. The shape of this array is `num_generations x population_size x sequence_length`.
- The `fitness_values` array stores the fitness values of each sequence produced by the algorithm during its execution. The shape of this array is `num_generations x population_size`.
- The `term_values` array stores the individual values of the terms composed to calculate the final fitness value of each generated sequence. The shape of this array is `num_generations x population_size x num_fitness_fn_terms`. 

### Example of a settings JSON file

The JSON file provided in `examples/terminal/settings.json` describe a genetic algorithm with the following characteristics:

- **Individual representation**: a string where each character represents a distinct amino acid type. 
- **Fitness function**: minimize the RMSD value after superimposing the modeled structure with the target structure.
- **Population size**: 5 individuals.
- **Selection**: by binary tournament (omitted in the `settings.json`).
- **Number of parents and children**: 2.
- **Recombination**: Uniform Crossover (omitted in the `settings.json`).
- **Mutation**: Random Resetting.
- **Mutation probability**: 0.1.
- **Replacement**: generational (omitted in the `settings.json`).
- **Maximum number of generations**: 2.

The omitted items in the JSON file will be filled with the default behavior as described above.

## Resuming EvoDesign from a previous execution

In case EvoDesign stopped running, you can resume it from the last checkpoint by running again the `python3 -m evodesign ...` command.
Just make sure you provide exactly the same input values as last time (i.e., same target PDB, settings JSON, output directory and jobname).
The program will automatically check if there is a saved checkpoint and resume from there. 

## Extending a finished EvoDesign execution

If EvoDesign finished running, but you wish to run it for additional generations, edit the settings JSON file and increase the `max_generations` value to the desired amount.
Then, re-run the `python3 -m evodesign ...` command.
The program will automatically detect that the algorithm ran for less generations than the value written in `max_generations` and continue the execution from the last generation.
