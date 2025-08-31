from pymoo.core.problem import Problem
from ...Fitness.FitnessFunction import FitnessFunction
from ...Prediction.Predictor import Predictor
from ...Utils.Chain import Chain, ChainFactory, numpy_sequence_to_str
import numpy as np
import numpy.typing as npt
from ...Utils.Exceptions import HttpGatewayTimeout
from requests.exceptions import ConnectTimeout


class MonoObjectiveCPD(Problem):

    def __init__(
        self,
        predictor: Predictor,
        fitness_fn: FitnessFunction,
        ref_chain: Chain,
        prediction_pdb_path: str = "prediction.pdb.tmp",
    ):
        self.fitness_fn = fitness_fn
        self.predictor = predictor
        self.ref_chain = ref_chain
        self.term_values = None
        self.prediction_pdb_path = prediction_pdb_path
        super().__init__(
            n_var=len(self.ref_chain.sequence),
            n_obj=1, # mono objective
            n_eq_constr=0,
            n_ieq_constr=0,
            xl=0,  # working on a 20-letters amino acid alphabet
            xu=19, # represented by integers from 0 to 19
            vtype=np.int64,
        )

    def compute_fitness(
        self, sequence: npt.NDArray[np.int64]
    ) -> npt.NDArray[np.float64]:
        sequence = numpy_sequence_to_str(sequence)
        while True:
            try:
                self.predictor.predict_pdb_file(sequence, self.prediction_pdb_path)
            except (HttpGatewayTimeout, ConnectTimeout):
                continue
            break
        model_chain = ChainFactory.create(self.prediction_pdb_path)
        return self.fitness_fn.do(model_chain, self.ref_chain)

    def _evaluate(self, x, out, *args, **kwargs) -> None:
        self.term_values = np.apply_along_axis(self.compute_fitness, 1, x)
        out["F"] = -1.0 * self.term_values[:, 0]
