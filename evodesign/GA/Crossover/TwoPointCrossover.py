from pymoo.operators.crossover.pntx import TwoPointCrossover as TPX

from .Crossover import Crossover


class TwoPointCrossover(Crossover):

    def __init__(self, prob: float = 1.0):
        super().__init__(prob, TPX(prob=prob))
        return
