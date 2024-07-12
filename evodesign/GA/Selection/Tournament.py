from .Selection import Selection
from ...Context import Context
import pandas as pd





class Tournament(Selection):
  
  def _params(self) -> dict:
    params = super()._params()
    params['elitism'] = self._elitism
    params['tournament_size'] = self._tournament_size
    return params
  
  
  
  def __init__(self,
               tournament_size: int,
               elitism: bool = False
               ) -> None:
    """
    Selection operator in which a random uniform sample of size 
    `tournament_size`, without replacement, is taken from the population, and 
    the individual with higher fitness from this sample is then selected. 

    Notice that it is possible for the same individual to be chosen multiple
    times. However, it is guaranteed that each consecutive pairs of 
    individuals will be distinct.

    Parameters
    ----------
    tournament_size : int
        The number of individuals to be randomly chosen to participate in 
        a tournament. Only one of these individuals will be chosen.
    elitism : bool, optional
        Indicates if this operator will be applying elitism or not. 
        Default is False.
    """
    super().__init__()
    self._elitism = elitism
    self._tournament_size = tournament_size



  def select_parents(self, 
                     population: pd.DataFrame,
                     context: Context
                     ) -> pd.DataFrame:
    """
    Selects a subset of individuals from the given population.

    Parameters
    ----------
    population : pandas.DataFrame
        The population to be sampled.
    context : Context
        The context data used by the calling evolutionary algorithm.

    Returns
    -------
    pandas.DataFrame
        The selected subset of individuals.
    """    
    selected_parents = pd.DataFrame(columns=population.columns)
    for i in range(len(population)):
      elitist = population.iloc[0]
      winner = self._tournament_selection(population, context, elitist)
      while i % 2 != 0 and \
          selected_parents.at[i - 1, 'sequence_id'] == winner['sequence_id']:
        winner = self._tournament_selection(population)
      selected_parents = pd.concat([ selected_parents, winner.to_frame().T ],
                                   ignore_index=True)
    return selected_parents
  


  def _tournament_selection(self, 
                            population: pd.DataFrame,
                            context: Context,
                            elitist: pd.Series = pd.Series()
                            ) -> pd.Series:
    tournament_size = self._tournament_size - 1 \
                      if self._elitism \
                      else self._tournament_size
    selection = context.rng.choice(population.index, 
                                  size=tournament_size, 
                                  replace=False)
    tournament = population.loc[selection]
    if self._elitism:
      tournament = pd.concat([ tournament, elitist.to_frame().T ])
    tournament.sort_values(by=context.sort_columns, 
                           ascending=context.sort_ascending,
                           inplace=True, 
                           ignore_index=True)
    return tournament.iloc[0]
    