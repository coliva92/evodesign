from .Selection import Selection
from typing import Optional, List
from ...Random import Random
import pandas as pd





class Tournament(Selection):

  @classmethod
  def _class_name(cls) -> str:
    return 'GA.Selection.TournamentElitism'
  


  def _params(self) -> dict:
    params = super()._params()
    params['tournamentSize'] = self._tournament_size
    params['fitnessColumns'] = self._fitness_columns
    params['ascendingSort'] = self._ascending
    return params
  
  
  
  def __init__(self,
               numCouples: int,
               tournamentSize: int,
               fitnessColumns: Optional[List[str]] = None,
               ascendingSort: Optional[List[bool]] = None
               ) -> None:
    """
    Selection operator in which a random uniform sample of size 
    `tournamentSize`, without replacement, is taken from the population, and 
    the individual with higher fitness from this sample is then selected. 
    This process is repeated `2 * numCouples` times to get the final 
    subset of individuals.

    Notice that it is possible for the same individual to be chosen multiple
    times. However, it is guaranteed that for every consecutive pair of 
    individuals in the selected subset will be distinct.

    Parameters
    ----------
    numCouples : int
        The number of parent couples to be selected from the population.
    tournamentSize : int
        The number of individuals to be randomly chosen to participate in 
        a tournament. Only one of these individuals will be chosen.
    fitnessColumns : List[str], optional
        The columns containing the fitness values to be used to select the 
        winner of each tournament. If `None`, then all columns which name
        has the 'fitness_' preffix will be used. Default is `None`.
    ascendingSort : List[bool], optional
        Each flag indicates if the corresponding column specified in 
        `fitnessColumns` should be sorted in asending order or in descending 
        order. If `None`, then all columns will be sorted in ascending order.
        Default is `None`. 
    """
    super().__init__(numCouples)
    self._tournament_size = tournamentSize
    self._fitness_columns = fitnessColumns
    self._ascending = ascendingSort 



  def select_parents(self, population: pd.DataFrame) -> pd.DataFrame:
    """
    Selects a subset of individuals from the given population.

    Parameters
    ----------
    population : pandas.DataFrame
        The population to be sampled.

    Returns
    -------
    pandas.DataFrame
        The selected subset of individuals.
    """
    if not self._fitness_columns or len(self._fitness_columns) == 0:
      self._fitness_columns = [ 
        col for col in population.columns if 'fitness' in col 
      ]
    if not self._ascending:
      self._ascending = len(self._fitness_columns) * [ True ]
    
    selected_parents = pd.DataFrame(columns=population.columns)
    for i in range(self._selection_size):
      elite = population.iloc[i]
      winner = self._tournament_selection(self._tournament_size, population)
      # garantee that two consecutive parents are different sequences
      while i % 2 != 0 and \
          selected_parents.at[i - 1, 'sequence'] == winner['sequence']:
        winner = self._tournament_selection(self._tournament_size, population)
      selected_parents = pd.concat([ selected_parents, winner.to_frame().T ],
                                   ignore_index=True)
    return selected_parents
  


  def _tournament_selection(self, 
                            selectionSize: int, 
                            population: pd.DataFrame):
    n = len(self._fitness_columns)
    if not self._fitness_columns or n == 0:
      self._fitness_columns = [ 
        col for col in population.columns if 'fitness' in col 
      ]
    if not self._ascending:
      self._ascending = len(self._fitness_columns) * [ True ]

    rng = Random.generator()
    selection = rng.choice(population.index, size=selectionSize, replace=False)
    tournament = population.loc[selection]
    if n == 1 and self._fitness_columns[0] == 'pandas.DataFrame.index':
      tournament.sort_index(inplace=True)
    else:
      # we assume pandas.DataFrame.index will never be used in multiobjective
      # optimization
      tournament.sort_values(by=self._fitness_columns, 
                             ascending=self._ascending,
                             inplace=True, 
                             ignore_index=True)
    return tournament.iloc[0]
