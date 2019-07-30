from typing import List, Dict, Iterator
from copy import deepcopy
from collections import OrderedDict
from cavecalc.configuration.model_parameters import read_parameters_file, ModelParameter, ParameterTypes
from cavecalc.file_utilities import read_config_from_csv, write_config_to_csv
import itertools


class RunConfig:
    """
    Encapsulates configuration settings for a single cavecalc model run (i.e. a collection
    of ModelParameters).

    This can be passed into CaveCalc to be executed.
    """
    required_params = set(p.Name for p in read_parameters_file())

    @classmethod
    def default(cls):
        config = cls(read_parameters_file())
        config._validate()
        return config

    @classmethod
    def from_file(cls, filepath: str = None):
        config = cls(read_parameters_file(filepath))
        config._validate()
        return config

    @classmethod
    def generate_suite(cls, all_combinations: bool = True, **kwargs) -> Iterator['RunConfig']:
        """
        Generate a suite of RunConfig objects based on requested run parameters.

        Unspecified parameters are populated with default values.

        Args:
            - all_combinations. If True (default), parameters will be combined to output
            RunConfigs for all possible permutations. If False, return only in-order comnbinations
            of listed parameters.

            - kwargs. Other arguments must match CaveCalc input parameter names. Values may be
            single values, or lists of multiple values. Multiple parameters may be given (examples
            below). All lists of multiple values must be of equal length if all_combinations is 
            False.

             By default, this function will return objects for all possible permutations of inputs.
             To combine options by list position, see the all_combinations arg.

        Example Usage:
            generate_suite()                                                        - 1 result
            generate_suite(temperature = 15)                                        - 1 result
            generate_suite(temperature = [15, 20], soil_pCO2 = 20000)               - 2 results
            generate_suite(temperature = [15, 20], soil_pCO2 = [15000, 20000])      - 4 results

            The final example would return only 2 results if all_combinations = False were passed
            alongside the existing arguments.
        """

        # Determine the function to cross-combine options as required
        fcomb = itertools.product if all_combinations else lambda *a : zip(*a)

        constant_params = {k: v for k, v in kwargs.items() if type(v) is not list}
        varying_params = {k: v for k, v in kwargs.items() if type(v) is list}

        varied_params = {}
        config_count = 1
        if any(varying_params):
            max_len = max(len(v) for v in varying_params.values())
            if not all_combinations and not all(len(v) == max_len for v in varying_params.values()):
                raise ValueError("All varying parameters must be of the same length.")

            # If required, calculate all permutations of parameters
            keys, values = zip(*varying_params.items())
            combined_values = list(fcomb(*values))
            varied_params = dict(zip(keys, zip(*combined_values)))
            config_count = len(combined_values)

        base_configs = (RunConfig.default() for _ in range(config_count))
        for i, config in enumerate(base_configs):
            for k, v in constant_params.items():
                config[k] = v
            for k, v in varied_params.items():
                config[k] = v[i]
            yield config

    def __init__(self, modelParameters: List[ModelParameter]):
        self._params: Dict[str, ModelParameter] = {p.Name: p for p in modelParameters}

    def __getitem__(self, key): 
        """
        RunConfig is subscriptable. This exposes the underlying _params dict.
        """
        try:
            return self._params[key].Value
        except KeyError:
            # check whether the key should have been present
            if key in self.required_params:
                raise KeyError(f"Parameter '{key}' should exist in the RunConfig, but is absent." +
                " This should not happen under normal usage - check your config setup.")
            else:
                raise KeyError(f"Parameter '{key}' does not exist. Typo?")

    def __setitem__(self, key, value):
        self._params[key].Value = value

    @property
    def params(self) -> Dict[str, ModelParameter]:
        self._validate()
        return deepcopy(self._params)

    @property
    def values(self) -> Dict[str, str]:
        self._validate()
        return deepcopy({k: p.Value for k, p in self.params.items()})

    def write_to_csv(self, filepath: str):
        return write_config_to_csv([p.print_format() for p in self.params.values()], filepath)
        
    def _validate(self) -> None:
        # validation of ModelParameters is handled by that class.
        # Here it is sufficient to validate we have all the necessary entries.
        if set(self._params.keys()) != self.required_params:
            raise ValueError(f"Run Config does not have the correct parameters")

    def __eq__(self, other: 'RunConfig') -> bool:
        """Two RunConfigs are equal if they contain all the same ModelParameters
        """
        if set(self._params.keys()) != set(other._params.keys()):
            return False

        for p1 in self._params.values():
            if p1 not in other._params.values():
                return False

        for p2 in other._params.values():
            if p2 not in self._params.values():
                return False
        
        return True