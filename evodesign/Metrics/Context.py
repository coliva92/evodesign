from .ContextInterface import ContextInterface
from .Metric import Metric
from ..Utils.Chain import Chain
from typing import Dict, Any


class Context(ContextInterface):

    def __init__(
        self,
        model_chain: Chain,
        ref_chain: Chain,
        metric_calculators: Dict[str, Metric],
        **kwargs
    ):
        self._model_chain = model_chain
        self._ref_chain = ref_chain
        self._metric_calculators = metric_calculators
        self._cached_results = {}
        self._extra_params = {**kwargs}

    def get_metric_components(self, metric_name: str) -> Dict[str, float]:
        if metric_name not in self._cached_results:
            calculator = self._metric_calculators[metric_name]
            self._cached_results[metric_name] = calculator.do_for_fitness_fn(self)
        return self._cached_results[metric_name]

    def get_component_value(self, component_name: str) -> float:
        # format for the `component_name` argument: "metric_name.component_name"
        name_parts = component_name.split(".")
        component_name = name_parts[-1]
        metric_name = ".".join(name_parts[:-1])
        components = self.get_metric_components(metric_name)
        return components[component_name]

    def get_model_chain(self) -> Chain:
        return self._model_chain

    def get_reference_chain(self) -> Chain:
        return self._ref_chain

    def get_extra_param_value(self, param_name: str):
        if param_name not in self._extra_params:
            return None
        return self._extra_params[param_name]

    def set_extra_param_value(self, param_name: str, param_value) -> None:
        self._extra_params[param_name] = param_value
