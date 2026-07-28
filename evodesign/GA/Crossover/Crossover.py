from abc import ABC

from pymoo.core.crossover import Crossover as PyMOOCrossover

from ...RetrievableSettings import RetrievableSettings




class Crossover(RetrievableSettings, ABC):

    def __init__(
        self,
        prob: float,
        pymoo_crossover: PyMOOCrossover,
    ) -> None:
        super().__init__()
        self.prob = prob
        self._pymoo_crossover = pymoo_crossover
        return
