from abc import ABC
from ...RetrievableSettings import RetrievableSettings
from pymoo.core.selection import Selection


class Selection(RetrievableSettings, ABC):

    def __init__(
        self,
        pymoo_selection: Selection,
    ):
        super().__init__()
        self._pymoo_selection = pymoo_selection
