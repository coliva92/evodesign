from typing import List
import evodesign.Sequence as Sequence
import evodesign.Utils as Utils
import pandas as pd





def create_random(size: int,
                  sequence_length: int,
                  generation_id: int = 0
                  ) -> pd.DataFrame:
    """
    Creates a table of randomly generated amino acid sequences. 
    Populations in EvoDesign are represented as Pandas DataFrames in which each 
    row corresponds to a single sequence, and where additional columns 
    typically represent fitness values for each sequence.

    Parameters
    ----------
    size : int
        The size of the population.
    sequence_length : int
        The length of each amino acid sequence in the population.
    generation_id : int, optional
        The unique identifier for the population being created. 
        The default value is 0.

    Returns
    -------
    pandas.DataFrame
        The generated table. Column `Sequence` contains the amino acid
        sequence of each individual in the population, while column 
        `Sequence_Id` contains a unique identifier for each sequence in the
        population. The `Survived` column contains a boolean flag which 
        can be used to indicate which sequences can move on to the next 
        generation.
    """
    sequences = [ 
        Sequence.create_random(sequence_length) 
        for _ in range(size) 
    ]
    return create(sequences, generation_id)
  


def create(sequences: List[str], 
           generation_id: int = 0
           ) -> pd.DataFrame:
    """
    Creates a new population DataFrame containing the specified amino acid 
    sequences.

    Parameters
    ----------
    sequences : List[str]
        The amino acid sequences that will be contained in the created 
        DataFrame. Each sequence must only contain letters corresponding to one
        of the 20 standard amino acids.
    generation_id : int
        The unique identifier for the population being created.

    Returns
    -------
    pandas.DataFrame
        A new DataFrame instance containing the specified population sequences.
    """
    data = {
      'generation_id': [],
      'sequence_id': [],
      'sequence': [],
      'survivor': []
    }
    for i, sequence in enumerate(sequences):
      seq_id = sequence_id(generation_id, i)
      data['generation_id'].append(generation_id)
      data['sequence_id'].append(seq_id)
      data['sequence'].append(sequence)
      data['survivor'].append(False)
    return pd.DataFrame(data)



def sequence_id(generation_id: int, 
                row_idx: int
                ) -> str:
    """
    Produces a unique identifier for some protein sequence.

    Parameters
    ----------
    generation_id : int
        The generation in the evolutionary algorithm at which the hypothetical 
        sequence was created.
    row_idx : int
        The row index for the hypothetical sequence in the corresponding 
        population DataFrame.

    Returns
    -------
    str
        A string containing a unique identifier to be used by a particular 
        protein sequence.
    """
    return f'prot_{Utils.pad_zeroes(generation_id)}_{Utils.pad_zeroes(row_idx)}'
  


def filename(population: pd.DataFrame) -> str:
    """
    Returns the full path to the CSV file where the given population ought to be 
    saved.

    Parameters
    ----------
    population : pd.DataFrame
        The population for which the file path will be produced.

    Returns
    -------
    str
        The file path for the given population.
    """
    generation_id = population.iloc[0]['generation_id']
    return f'pop_{Utils.pad_zeroes(generation_id)}.csv'
