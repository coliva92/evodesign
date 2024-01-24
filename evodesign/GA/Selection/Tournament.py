from Selection import Selection
from ...Random import Random
import pandas as pd





class Tournament(Selection):

  @classmethod
  def _name(cls) -> str:
    return 'GA_Selection_Tournament'
  


  def _params(self) -> dict:
    params = super()._params()
    params['tournamentSize'] = self._tournament_size
    return params
  
  
  
  def __init__(self,
               selectionSize: int,
               tournamentSize: int
               ) -> None:
    """
    Selection operator in which a random uniform sample of size 
    `tournamentSize`, without replacement, is taken from the population, and 
    the individual with higher fitness from this sample is then selected. 
    This process is repeated `selectionSize` times to get the final subset 
    of individuals.

    Notice that it is possible for the same individual to be chosen multiple
    times. However, it is guaranteed that for every consecutive pair of 
    individuals in the selected subset will be distinct.

    Parameters
    ----------
    selectionSize : int
        The number of individuals to be selected from the population.
    tournamentSize : int
        The number of individuals to be randomly chosen to participate in 
        a tournament. Only one of these individuals will be chosen.
    """
    super().__init__(selectionSize)
    self._tournament_size = tournamentSize



  def select_parents(self, population: pd.DataFrame) -> pd.DataFrame:
    
    def tournament_selection(k: int, population: pd.DataFrame):
      selection = rng.choice(population.index, size=k, replace=False)
      tournament = population.iloc[selection]
      tournament.sort_values(by='Fitness', ascending=False, inplace=True)
      return tournament.iloc[0]
    
    selected_parents = pd.DataFrame(columns=population.columns)
    rng = Random.generator()
    for i in range(self._selection_size):
      winner = tournament_selection(self._tournament_size, population)
      # garantee that two consecutive parents are different sequences
      while i % 2 != 0 and \
          selected_parents.iloc[i - 1, 'Sequence'] == winner['Sequence']:
        winner = tournament_selection(self._tournament_size, population)
      selected_parents = pd.concat([ selected_parents, pd.DataFrame(winner) ],
                                   ignore_index=True)
    return selected_parents
    