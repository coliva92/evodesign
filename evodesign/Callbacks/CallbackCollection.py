from pymoo.core.callback import Callback
from pymoo.core.algorithm import Algorithm as PyMOOAlgorithm
from typing import List


class CallbackCollection(Callback):

    def __init__(self, callbacks: List[Callback]):
        super().__init__()
        self.callbacks = callbacks
        return
    
    def append(self, callback: Callback) -> None:
        self.callbacks.append(callback)
        return

    def notify(
        self, 
        algorithm: PyMOOAlgorithm,
    ) -> None:
        for callback in self.callbacks:
            callback.notify(algorithm)
        return
