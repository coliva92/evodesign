from .Selection import Selection
import evodesign.Random as Random
import pandas as pd





class Overselection(Selection):

  @classmethod
  def _class_name(cls) -> str:
    return 'GA.Selection.Overselection'
  


  def _params(self) -> dict:
    return {
      'upperSize': self._upper_size,
      'upperProb': self._upper_prob,
      'lowerProb': self._lower_prob
    }
  


  def __init__(self, 
               upperSize: int,
               upperProb: float = 0.8,
               lowerProb: float = 0.2,
               twoChildren: bool = True,) -> None:
    """
    Selection operator where the population is divided into to groups, which
    we call the "upper bin" and "lower bin". The upper bin contains the top
    `upperSize` individuals according to their fitness, and the lower bin
    contains all remaining individuals. Then, individuals are selected in pairs,
    where each pair could be formed by two randomly selected individuals from
    the upper bin, two from the lower bin, or one from each bin. Which of 
    these three options will be taken for a given pair is chosen randomly. 
    All individuals are selected from their bins with a uniform distribution and
    without replacement.

    Parameters
    ----------
    upperSize : int
        The number of individuals in the upper bin.
    upperProb : float, optional
        The probability for selecting a pair of individuals from the upper bin.
        The default is 0.8.
    lowerProb : float, optional
        The probability for selecting a pair of individuals from the lower bin.
        The default is 0.2. The probability for selecting  mixed pair of 
        individuals is always 1.0 - upperProb - lowerProb.
    twoChildren : bool, optional.
        A flag that indicates if the recombination operation would produce two
        children or only one child per parent pair. Depending on the case, the
        number of selections performed by this operator will be different in
        order to produce the correct number of individuals for the next 
        generation. The default is True.
    """
    self._two_children = twoChildren
    self._upper_size = upperSize
    self._upper_prob = upperProb
    self._lower_prob = lowerProb
    self._weights = [
      upperProb,
      lowerProb,
      1. - upperProb - lowerProb
    ]
  


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
    selection_size = len(population) \
                     if self._two_children \
                     else 2 * len(population)
    
    selected_parents = pd.DataFrame(columns=population.columns)
    rng = Random.generator()
    upper_bin = population.iloc[:self._upper_size]
    lower_bin = population.iloc[self._upper_size:]
    while len(selected_parents) < selection_size:
      option = rng.choice([ 0, 1, 2 ], p=self._weights)
      if option == 0:
        selection = rng.choice(upper_bin.index, size=2, replace=False)
        parents = upper_bin.loc[selection]
      if option == 1:
        selection = rng.choice(lower_bin.index, size=2, replace=False)
        parents = lower_bin.loc[selection]
      if option == 2:
        m = rng.choice(upper_bin.index)
        f = rng.choice(lower_bin.index)
        mother, father = population.loc[m], population.loc[f]
        parents = pd.concat([ mother, father ], axis=1, ignore_index=True).transpose()
      while parents.iloc[0]['sequence_id'] == parents.iloc[1]['sequence_id']:
        if option == 0:
          f = rng.choice(upper_bin.index)
          father = upper_bin.loc[f]
        else:
          f = rng.choice(lower_bin.index)
          father = lower_bin.loc[f]
        parents.iloc[1] = father
      selected_parents = pd.concat([ selected_parents, parents ],
                                   ignore_index=True)
    return selected_parents
