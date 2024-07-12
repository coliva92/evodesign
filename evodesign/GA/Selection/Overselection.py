from .Selection import Selection
from ...Context import Context
import pandas as pd





class Overselection(Selection):
  
  def _params(self) -> dict:
    params = super()._params()
    params['upper_size'] = self._upper_size
    params['upper_prob'] = self._upper_prob
    params['lower_prob'] = self._lower_prob
    params['two_children'] = self._two_children
    return params
  


  def __init__(self, 
               upper_size: int,
               upper_prob: float = 0.8,
               lower_prob: float = 0.2,
               two_children: bool = True
               ) -> None:
    """
    Selection operator where the population is divided into to groups, which
    we call the "upper bin" and "lower bin". The upper bin contains the top
    `upper_size` individuals according to their fitness, and the lower bin
    contains all remaining individuals. Then, individuals are selected in pairs,
    where each pair could be formed by two randomly selected individuals from
    the upper bin, two from the lower bin, or one from each bin. Which of 
    these three options will be taken for a given pair is chosen randomly. 
    All individuals are selected from their bins with a uniform distribution and
    without replacement.

    Parameters
    ----------
    upper_size : int
        The number of individuals in the upper bin.
    upper_prob : float, optional
        The probability for selecting a pair of individuals from the upper bin.
        The default is 0.8.
    lower_prob : float, optional
        The probability for selecting a pair of individuals from the lower bin.
        The default is 0.2. The probability for selecting  mixed pair of 
        individuals is always 1.0 - upper_prob - lower_prob.
    two_children : bool, optional.
        A flag that indicates if the recombination operation would produce two
        children or only one child per parent pair. Depending on the case, the
        number of selections performed by this operator will be different in
        order to produce the correct number of individuals for the next 
        generation. The default is True.
    """
    super().__init__()
    self._two_children = two_children
    self._upper_size = upper_size
    self._upper_prob = upper_prob
    self._lower_prob = lower_prob
    self._weights = [
      upper_prob,
      lower_prob,
      1. - upper_prob - lower_prob
    ]
  


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
    selection_size = len(population) \
                     if self._two_children \
                     else 2 * len(population)
    
    selected_parents = pd.DataFrame(columns=population.columns)
    upper_bin = population.iloc[:self._upper_size]
    lower_bin = population.iloc[self._upper_size:]
    while len(selected_parents) < selection_size:
      option = context.rng.choice([ 0, 1, 2 ], p=self._weights)
      if option == 0:
        selection = context.rng.choice(upper_bin.index, size=2, replace=False)
        parents = upper_bin.loc[selection]
      if option == 1:
        selection = context.rng.choice(lower_bin.index, size=2, replace=False)
        parents = lower_bin.loc[selection]
      if option == 2:
        m = context.rng.choice(upper_bin.index)
        f = context.rng.choice(lower_bin.index)
        mother, father = population.loc[m], population.loc[f]
        parents = pd.concat([ mother, father ], axis=1, ignore_index=True).transpose()
      while parents.iloc[0]['sequence_id'] == parents.iloc[1]['sequence_id']:
        if option == 0:
          f = context.rng.choice(upper_bin.index)
          father = upper_bin.loc[f]
        else:
          f = context.rng.choice(lower_bin.index)
          father = lower_bin.loc[f]
        parents.iloc[1] = father
      selected_parents = pd.concat([ selected_parents, parents ],
                                   ignore_index=True)
    return selected_parents
