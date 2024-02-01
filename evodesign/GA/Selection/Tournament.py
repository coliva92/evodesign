from Selection import Selection
from ...Random import Random
import pandas as pd





class Tournament(Selection):

  @classmethod
  def _class_name(cls) -> str:
    return 'GA.Selection.Tournament'
  


  def _params(self) -> dict:
    params = super()._params()
    params['tournamentSize'] = self._tournament_size
    params['fitnessColumn'] = self._fitness_column
    return params
  
  
  
  def __init__(self,
               numSelectedCouples: int,
               tournamentSize: int,
               fitnessColumn: str | None = None
               ) -> None:
    """
    Selection operator in which a random uniform sample of size 
    `tournamentSize`, without replacement, is taken from the population, and 
    the individual with higher fitness from this sample is then selected. 
    This process is repeated `2 * numSelectedCouples` times to get the final 
    subset of individuals.

    Notice that it is possible for the same individual to be chosen multiple
    times. However, it is guaranteed that for every consecutive pair of 
    individuals in the selected subset will be distinct.

    Parameters
    ----------
    numSelectedCouples : int
        The number of parent couples to be selected from the population.
    tournamentSize : int
        The number of individuals to be randomly chosen to participate in 
        a tournament. Only one of these individuals will be chosen.
    fitnessColumn : str, optional
        The column containing the fitness value to be used to select the winner
        of each tournament. If `None`, then the rightmost column containing
        the 'fitness_' suffix will be used. Default is `None`.
    """
    super().__init__(numSelectedCouples)
    self._tournament_size = tournamentSize
    self._fitness_column = fitnessColumn



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
    selected_parents = pd.DataFrame(columns=population.columns)
    for i in range(self._selection_size):
      winner = self._tournament_selection(self._tournament_size, population)
      # garantee that two consecutive parents are different sequences
      while i % 2 != 0 and \
          selected_parents.iat[i - 1, 'sequence'] == winner['sequence']:
        winner = self._tournament_selection(self._tournament_size, population)
      selected_parents = pd.concat([ selected_parents, pd.DataFrame(winner) ],
                                   ignore_index=True)
    return selected_parents
  


  def _tournament_selection(self, 
                            selectionSize: int, 
                            population: pd.DataFrame):
    if not self._fitness_column:
      temp = [ col for col in population.columns if 'fitness' in col ]
      self._fitness_column = temp[-1]
    rng = Random.generator()
    selection = rng.choice(population.index, size=selectionSize, replace=False)
    tournament = population.loc[selection]
    tournament.sort_values(by=self._fitness_column, ascending=False)
    return tournament.iloc[0]
    