from ..Metric import Metric
from Bio.PDB.Atom import Atom
from pandas import Series
from evodesign.Context import Context
from typing import Optional, List
import numpy as np





class Sigmoid(Metric):

    def _params(self) -> dict:
        params = super()._params()
        params["metric"] = self._metric.settings()
        params["scaling_factor"] = self._scaling_factor
        return params
    


    def __init__(self, 
                 metric: Metric,
                 x_scaling_factor: float = -1.0,
                 column: Optional[str] = None
                 ) -> None:
        super().__init__(column)
        self._x_scaling_factor = x_scaling_factor
        self._metric = metric



    def _compute_values(self, 
                        backbone: List[Atom], 
                        data: Series, 
                        context: Context
                        ) -> Series:
        # compute the metric if not already computed
        data = self._metric(backbone, data, context)
        value = data[self._metric.column_name()]
        data[self.column_name()] = self._normalize(value)
        return data



    def _normalize(self, x: float) -> float:
        return 1.0 / (1.0 + np.exp(self._x_scaling_factor * x))
