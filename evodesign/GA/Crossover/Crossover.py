from abc import ABC
from ...RetrievableSettings import RetrievableSettings
from pymoo.core.crossover import Crossover as PyMOOCrossover


class Crossover(RetrievableSettings, ABC):

    def __init__(self, prob: float, pymoo_crossover: PyMOOCrossover):
        super().__init__()
        self.prob = prob
        self._pymoo_crossover = pymoo_crossover
