from abc import ABC, abstractmethod
from ..Utils.Chain import Chain
from typing import Dict


class ContextInterface(ABC):

    @abstractmethod
    def get_metric_components(
        self,
        metric_name: str,
    ) -> Dict[str, float]:
        raise NotImplementedError

    @abstractmethod
    def get_component_value(
        self,
        component_name: str,
    ) -> float:
        raise NotImplementedError

    @abstractmethod
    def get_model_chain(self) -> Chain:
        raise NotImplementedError

    @abstractmethod
    def get_reference_chain(self) -> Chain:
        raise NotImplementedError

    @abstractmethod
    def get_extra_param_value(
        self,
        param_name: str,
    ):
        raise NotImplementedError

    @abstractmethod
    def set_extra_param_value(
        self,
        param_name: str,
        param_value,
    ) -> None:
        raise NotImplementedError
