from typing import List
import evodesign.Sequence as s
import evodesign.Utils as Utils
import pandas as pd





class Population:
  
  @classmethod
  def create_random(cls, 
                    size: int,
                    sequenceLength: int,
                    generationId: int = 0
                    ) -> pd.DataFrame:
    """
    Creates a table of randomly generated amino acid sequences. 
    Populations in EvoDesign are Pandas DataFrames in which each row corresponds 
    to a single sequence, and where additional columns typically represent
    fitness values for each sequence.

    Parameters
    ----------
    size : int
        The size of the population.
    sequenceLength : int
        The length of each amino acid sequence in the population.
    generationId : int, optional
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
      s.create_random(sequenceLength) 
      for _ in range(size) 
    ]
    return cls.create(sequences, generationId)
  


  @classmethod
  def create(cls, 
             sequences: List[str], 
             generationId: int = 0
             ) -> pd.DataFrame:
    """
    Creates a new population DataFrame containing the specified amino acid 
    sequences.

    Parameters
    ----------
    sequences : List[str]
        The amino acid sequences that will be contained in the created 
        DataFrame. Each sequence must only contain letters corresponding to one
        of the 20 essential amino acids.
    generationId : int
        The unique identifier for the population being created.

    Returns
    -------
    pandas.DataFrame
        _description_
    """
    data = {
      'generation_id': [],
      'sequence_id': [],
      'sequence': [],
      'survivor': []
    }
    for i, sequence in enumerate(sequences):
      seq_id = cls._sequence_id(generationId, i)
      data['generation_id'].append(generationId)
      data['sequence_id'].append(seq_id)
      data['sequence'].append(sequence)
      data['survivor'].append(False)
    return pd.DataFrame(data)
  


  @classmethod
  def filename(cls, population: pd.DataFrame) -> str:
    """
    Returns the name of the file where the given population would be saved.

    Parameters
    ----------
    population : pd.DataFrame
        The population for which the file name will be returned.

    Returns
    -------
    str
        The filename for the given population.
    """
    generation_id = population.iloc[0]['generation_id']
    return f'pop_{Utils.pad_zeroes(generation_id)}.csv'
  


  @classmethod
  def _sequence_id(cls, 
                   generationId: int, 
                   rowIdx: int
                   ) -> str:
    return f'prot_{Utils.pad_zeroes(generationId)}_{Utils.pad_zeroes(rowIdx)}'
