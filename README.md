# EvoDesign

This is a rudimentary framework written in Python 3.8 for implementing and testing different evolutionary algorithms for protein design.
In the protein design problem, the goal is to find the amino acid sequence that folds into a given target structure. 
To determine if a given sequence would indeed fold into the desired structure, the evolutionary algorithms implemented with EvoDesign use a structure prediction model (e.g. AlphaFold2 or ESMFold).

## Content

1. [Installation](#installation)
2. [Using EvoDesign as a module from the console](#module-usage)
  1. [Example of a settings JSON file](#example)
  2. [Resuming execution from an interruption](#resuming-execution)
3. [Using EvoDesign as a library](#library-usage)

<a name="installation"></a>
## Installation

For running EvoDesign, it is highly recommended to install [Anaconda](https://www.anaconda.com/).
Using the `conda` command, you need to install the following packages:

- [NumPy](https://numpy.org/install/).
- [BioPython](https://biopython.org/wiki/Packages).
- [Pandas](https://pypi.org/project/pandas/).
- [Matplotlib](https://matplotlib.org/stable/users/getting_started/index.html#installation-quick-start).
- [SciPy](https://scipy.org/install/#pip-install).

Next, clone this repository using the following command:

```
git clone https://github.com/coliva92/evodesign.git
```

Finally, it is recommended that you add the `evodesign` repo folder to the `PYTHONPATH` environment variable so you can run the program from any directory in your system. 
If you're using bash, you can achieve this by running the following command:

```
export PYTHONPATH=${PYTHONPATH}:path/to/evodesign/repo
source ~/.bashrc
```

Installation instructions for the external components used by EvoDesign are provided separately:

- [AlphaFold2](https://github.com/google-deepmind/alphafold).
- [ESM2 and ESMFold](https://github.com/facebookresearch/esm).
- [iLearn](https://github.com/Superzchen/iLearn).
- [PyRosetta](https://www.pyrosetta.org/downloads#h.iwt5ktel05jc).

<a name="module-usage"></a>
## Using EvoDesign as a module from the console

You can run EvoDesign as a Python module from the command line interface. You can check the command instructions using:

```
python -m evodesign --help
```

Follow these steps to run a custom evolutionary algorithm using EvoDesign:

1. Write a JSON file describing the characteristics of the algorithm to be executed. 
2. Gather the PDB and FASTA files of the target protein. 
3. Run the following command: 

```
python -m evodesign path/to/settings.json \
                    path/to/target.pdb \
                    path/to/output/folder \
                    -f path/to/target.fasta 
```

**NOTE**: if the output folder, here on refered to as the _workspace_, does not exists, it will be created automathically. 

After finishing execution, the workspace folder will contain the following structure:

```
.
├─ workspace_name/
|  ├─ pdbs/
|  |  ├─ prot_0001_0000.pdb
|  |  ├─ ...
|  ├─ populations/
|  |  ├─ pop_0001.csv
|  |  ├─ ...
|  ├─ commit_hash.txt
|  ├─ initial_rng_state.json
|  ├─ settings.json
|  ├─ statistics.csv
|  ├─ target_protein.pdb
```

- The `pdbs` folder contains the PDB files for all predicted structures produced during algorithm execution. Each file is named `prot_XXXX_YYYY.pdb`, where `XXXX` is an integer with trailing zeroes indicating the generation number and `YYYY` is also an integer with trailing zeroes indicating the individual number in the given generation.
- The `populations` folder contains the population data for every generation. Each file is named `pop_XXXX`, where `XXXX` is an integer with trailing zeroes indicating the generation number.
- The `commit_hash.txt` file simply contains the number of the specific EvoDesign version used to produced all the data in the current workspace folder. This version number is simply the hash value of a specific commit in the EvoDesign GitHub repository.
- The `initial_rng_state.json` contains data needed to recreate the state that the RNG had at the very first generation. This is useful for reproducing the results of a given execution. See [here](https://numpy.org/doc/stable/reference/random/bit_generators/pcg64.html) for more details.
- The `settings.json` file is a copy of the settings file provided as input of the `python -m evodesign` command.
- Lastly, a copy of the target protein PDB file will be stored with the original name. In this example, this file is named `target_protein.pdb`.

Additional folders and files may be created if external components (e.g. iLearn or ESM2) are used, but the structure explained here should be common for all cases.

<a name="example"></a>
### Example of a settings JSON file

Let's say you want to run a simple genetic algorithm with the following characteristics:

- **Individual representation**: a string where each character represents a distinct amino acid type. 
- **Fitness function**: maximize the GDT after superimposing the modeled structure with the target structure.
- **Population size**: 100 individuals.
- **Selection**: by binary tournament.
- **Number of parents and children**: 2.
- **Recombination**: by 2-points crossover.
- **Mutation**: randomly choose 1 amino acid in the string and swap it for another one randomly chosen.
- **Mutation probability**: 10%.
- **Replacement**: combine the parents and children populations and preserve the 100 fittest individuals.
- **Maximum number of generations**: 100.

All these characteristics can be described in a settings JSON file like the following:

```json
{
    "Algorithms.GASteadyState": {
        "max_generations": 100,
        "population_size": 100,
        "predictor": {
            "Prediction.ESMFoldRemoteApi": {}
        },
        "selection": {
            "GA.Selection.Tournament": {
                "tournament_size": 2
            }
        },
        "recombination": {
            "GA.Recombination.TwoPointsCrossover": {}
        },
        "mutation": {
            "GA.Mutation.Swap": {
                "mutation_prob": 0.1,
                "num_swaps": 1
            }
        },
        "metrics": [
            {
                "Metrics.Gdt": {}
            }
        ],
        "fitnessFn": {
            "Fitness.LinearCombination": {
                "upper_bound": 0.95,
                "metric_columns": [ "Metrics.Gdt" ]
            }
        },
        "sort_columns": [ "Fitness.LinearCombination", "Metrics.Rmsd", "plddt" ],
        "sort_ascending": [ false, true, false ]
    }
}
```

Once the settings JSON file has been written, EvoDesign can be executed using the command shown in the previous section.

<a name="resuming-execution"></a>
### Resuming execution from an interruption

Algorithms run by EvoDesign tipycally take a long time and execution may be interrupted either manually or due to some error. 
The EvoDesign saves its progress after computing the fitness of any individual.
Thus, the execution can be resumed from its last saved point by simply re-running the same command shown in the [Using EvoDesign as a module from the console](#module-usage) section.

<a name="library-usage"></a>
## Using EvoDesign as a library

EvoDesign can also be imported as a library into a custom Python script.
In this case, EvoDesign code is actually written emulating the settings JSON file. 
For example, the settings described in the [Example of a settings JSON file](#example) could be reproduced in a Python script as follows:

```python
from evodesign.Algorithms import GASteadyState
from evodesign.Prediction.ESMFoldRemoteApi import ESMFoldRemoteApi
from evodesign.GA.Selection.Tournament import Tournament
from evodesign.GA.Recombination.TwoPointsCrossover import TwoPointsCrossover
from evodesign.GA.Mutation.Swap import Swap
from evodesign.Metrics.Gdt import Gdt
from evodesign.Fitness.LinearCombination import LinearCombination
from evodesign.Context import Context

algorithm = GASteadyState(max_generations=100,
                          population_size=100,
                          predictor=ESMFoldRemoteApi(),
                          selection=Tournament(tournament_size=2),
                          recombination=TwoPointsCrossover(),
                          mutation=Swap(mutation_prob=0.1,
                                        num_swaps=1),
                          metrics=[ Gdt() ],
                          fitness_fn=LinearCombination(upper_bound=0.95,
                                                       metric_columns=[ 'Metrics.Gdt' ]),
                          sort_columns=[ 
                            "Fitness.LinearCombination", "Metrics.Rmsd", "plddt" 
                          ],
                          sort_ascending=[ False, True, False ])
context = Context.create('target_protein.pdb', 'target_protein.fasta')
algorithm.setup(context, 'path/to/workspace')
algorithm() # begin algorithm execution
```
