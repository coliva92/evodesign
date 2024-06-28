from .Selection import Selection
from typing import Optional, List
import evodesign.Random as Random
import pandas as pd





class Tournament(Selection):
  
  def _params(self) -> dict:
    return {
      'elitism': self._elitism,
      'tournamentSize': self._tournament_size,
      'sortColumns': self._sort_columns,
      'sortAscending': self._ascending
    }
  
  
  
  def __init__(self,
               tournamentSize: int,
               sortColumns: Optional[List[str]] = None,
               sortAscending: Optional[List[bool]] = None,
               elitism: bool = False
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
    tournamentSize : int
        The number of individuals to be randomly chosen to participate in 
        a tournament. Only one of these individuals will be chosen.
    sortColumns : List[str], optional
        The columns containing the fitness values to be used to select the 
        winner of each tournament. If `None`, then all columns which name
        has the 'fitness_' preffix will be used. Default is `None`.
    sortAscending : List[bool], optional
        Each flag indicates if the corresponding column specified in 
        `sortColumns` should be sorted in asending order or in descending 
        order. If `None`, then all columns will be sorted in ascending order.
        Default is `None`. 
    elitism : bool, optional
        Indicates if this operator will be applying elitism or not. 
        Default is False.
    """
    # TODO actualizar documentación
    super().__init__()
    self._elitism = elitism
    self._tournament_size = tournamentSize
    self._sort_columns = sortColumns
    self._ascending = sortAscending 



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
    if not self._sort_columns:
      self._sort_columns = [ 
        col for col in population.columns if 'fitness_' in col 
      ]
    if not self._ascending:
      self._ascending = len(self._sort_columns) * [ False ]
    
    selected_parents = pd.DataFrame(columns=population.columns)
    for i in range(len(population)):
      elitist = population.iloc[0]
      winner = self._tournament_selection(population, elitist)
      while i % 2 != 0 and \
          selected_parents.at[i - 1, 'sequence_id'] == winner['sequence_id']:
        winner = self._tournament_selection(population)
      selected_parents = pd.concat([ selected_parents, winner.to_frame().T ],
                                   ignore_index=True)
    return selected_parents
  


  def _tournament_selection(self, 
                            population: pd.DataFrame,
                            elitist: pd.Series = pd.Series()):
    rng = Random.generator()
    tournament_size = self._tournament_size - 1 \
                      if self._elitism \
                      else self._tournament_size
    selection = rng.choice(population.index, 
                           size=tournament_size, 
                           replace=False)
    tournament = population.loc[selection]
    if self._elitism:
      tournament = pd.concat([ tournament, elitist.to_frame().T ])
    tournament.sort_values(by=self._sort_columns, 
                           ascending=self._ascending,
                           inplace=True, 
                           ignore_index=True)
    return tournament.iloc[0]
    