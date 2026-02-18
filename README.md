# EvoDesign

This is a Python framework for implementing and testing different evolutionary algorithms for protein design.
In the context of this framework, the goal of protein design is to find the amino acid sequence that folds into a given target structure and satisfies other design requirements. 
To determine if a given sequence would indeed fold into the desired structure, the evolutionary algorithms implemented in EvoDesign use a deep learning structure prediction model (e.g. AlphaFold, ESMFold, etc.).
EvoDesign is built on top of the [PyMOO library](https://pymoo.org/).

## Installation

It is highly recommended to install the framework's dependencies using [Anaconda](https://www.anaconda.com/).
Once Anaconda is installed and activated, clone this repository into your system and run the following command in the terminal:

```
conda env create -f <path to the cloned repository>/environment.yml
```

This command will create an Anaconda environment named `evodesign` with all the dependencies required by EvoDesign already installed.
This environment can then be activated using the following command in the terminal:

```
conda activate evodesign
```

Additionally, to be able to call EvoDesign's functions from any directory in your local file system, you need to add the path to the evodesign source code to the `PYTHONPATH` environment variable.
To do this, use the following commands in the terminal:

```
export PYTHONPATH=${PYTHONPATH}:<path to the cloned repository>/evodesign
source ~/.bashrc
```

To make these changes permanent, you can manually edit the `~/.bashrc` file and add the `export PYTHONPATH=...` command into that file.
If you don't save the `~/.bashrc` file with these changes, you will have to edit the `PYTHONPATH` environment variable every time you turn on the computer. 

Lastly, all the installation instructions just described can be done automatically by running the `install_evodesign.sh` script provided in this repository.
To run it, simply do the follwing in the terminal:

```
./install_evodesign.sh
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
After a model is installed in its own environment, we recommend installing EvoDesign's dependencies in that same environment. 
This can be done using the `environment.yml` provided in this repository.
For example, say you installed AlphaFold 3 in an environment called `alphafold3`. 
After finishing the installation, EvoDesign can be installed into the `alphafold3` environment using these commands in the terminal:

```
conda activate alphafold3
conda env update -f <path to the evodesign respository>/environment.yml
```

**To install the iLearn descriptors, simply clone the official GitHub repository.** After that, modify the `PYTHONPATH` environment variable to add the path to the repository.
For example, say you cloned the iLearn repository to the `Documents` folder in your home directory.
Use these commands to modify the `PYTHONPATH`:

```
export PYTHONPATH=${PYTHONPATH}:~/Documents/iLearn
source ~/.bashrc
```

You can also edit the `~/.bashrc` file to make these changes permanent.

## Using EvoDesign from the terminal

EvoDesign is typically run from the terminal. 
Assuming the path to EvoDesign's source code was added to the `PYTHONPATH`, you can check the instructions using this command:

```
python -m evodesign --help
```

The `examples` folder contain example code showcasing how to use EvoDesign.
To run a custom evolutionary algorithm using EvoDesign, follow these steps:

1. Write a JSON file describing the characteristics of the algorithm to be executed (check the `examples/terminal/settings.json` for an example). 
2. Gather the PDB file of the target protein. 
3. Run the following command: 

```
python -m evodesign path/to/target.pdb \
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
|  ├─ git_commit_hash.txt
|  ├─ initial_rng_state.txt
|  ├─ last_rng_state.txt
|  ├─ pymoo_algorithm.bin
|  ├─ results.npz
|  ├─ settings.json
|  ├─ target_protein.pdb
```

These files are stored mainly for reproduction purposes. 

- The `git_commit_hash.txt` file contains the commit hash corresponding to the specific EvoDesign version used to produced these results.
- The `initial_rng_state.txt` and `last_rng_state.txt` files contain the RNG state at the beginning and at the end of the last iteration of the algorithm's execution, respectively.
- The `pymoo_algorithm.bin` file stores the state of the PyMOO evolutionary algorithm in the last iteration.
- The `results.npz` file stores the sequences and fitness values produced by the evolutionary algorithm in each iteration. See below for more details.
- The `settings.json` file is a copy of the settings file that was provided as input to the `python -m evodesign` command.
- The `*.pdb` file is a copy of the target protein PDB file that was provided as input to the `python -m evodesign` command. In this example, this file is named `target_protein.pdb`.

By default, these files are updated every 10 iterations.

### The `results.npz` file

For classic genetic algorithms, the `results.npz` file stores three NumPy arrays called `generations`, `fitness_values` and `term_values`.

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

In case EvoDesign stopped running, you can resume it from the last checkpoint by running again the `python -m evodesign ...` command.
Just make sure you provide exactly the same input values as last time (i.e., same target PDB, settings JSON, output directory and jobname).
The program will automatically check if there is a saved checkpoint and resume from there. 

## Extending a finished EvoDesign execution

If EvoDesign finished running, but you wish to run it for additional generations, edit the settings JSON file and increase the `max_generations` value to the desired amount.
Then, re-run the `python -m evodesign ...` command.
The program will automatically detect that the algorithm ran for less generations than the value written in `max_generations` and continue the execution from the last generation.
