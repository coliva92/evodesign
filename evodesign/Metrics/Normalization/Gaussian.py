from .Normalization import Normalization
import numpy as np


class Gaussian(Normalization):

    def __init__(self, 
                 mean_val: float = 1.0, 
                 stdev: float = 1.0, 
                 scaling_factor: float = 2.0
                 ) -> None:
        super().__init__()
        self.mean_val = mean_val
        self.stdev = stdev
        self.scaling_factor = scaling_factor
        return
    
    def do(self, x: float) -> float:
        return np.exp(-((x - self.mean_val)**2) / (self.scaling_factor * np.pi * self.stdev**2))
    