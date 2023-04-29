__author__ = "Mateus Stano Junqueira"
__copyright__ = "Copyright 20XX, RocketPy Team"
__license__ = "MIT"

from typing import Any, Callable, List, Tuple, Union

from pydantic import Field, StrictFloat, StrictInt, StrictStr

from ..Parachute import Parachute
from .DispersionModel import DispersionModel


class McParachute(DispersionModel):
    """Monte Carlo Parachute class, used to validate the input parameters of the
    parachute, based on the pydantic library. It uses the DispersionModel class
    as a base class, see its documentation for more information. The inputs
    defined here correspond to the ones defined in the Parachute class.
    """

    # Field(...) means it is a required field, exclude=True removes it from the
    # self.dict() method, which is used to convert the class to a dictionary
    # Fields with typing Any must have the standard dispersion form of tuple or
    # list. This is checked in the DispersionModel @root_validator
    # Fields with typing that is not Any have special requirements
    parachute: Parachute = Field(..., exclude=True)
    CdS: Any = 0
    trigger: List[Union[Callable, None]] = []
    samplingRate: Any = 0
    lag: Any = 0
    name: List[Union[StrictStr, None]] = []
    noise: List[
        Union[
            Tuple[
                Union[StrictInt, StrictFloat],
                Union[StrictInt, StrictFloat],
                Union[StrictInt, StrictFloat],
            ],
            None,
        ]
    ] = []
