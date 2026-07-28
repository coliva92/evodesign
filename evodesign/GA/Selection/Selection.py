from abc import ABC

from pymoo.core.selection import Selection

from ...RetrievableSettings import RetrievableSettings


class Selection(RetrievableSettings, ABC):

    def __init__(
        self,
        pymoo_selection: Selection,
    ):
        super().__init__()
        self._pymoo_selection = pymoo_selection
