import csv
import os.path
from dataclasses import dataclass
from typing import List, Set, Dict
from copy import deepcopy
from enum import Enum

default_parameters_file = 'model_parameters.csv'
file_directory = os.path.dirname(__file__)
default_filepath = os.path.join(file_directory, default_parameters_file)


class ParameterTypes(Enum):
    NUMERIC = 'numeric'
    STRING = 'str'
    BOOLEAN = 'bool'

    @classmethod
    def supported_types(cls) -> Set[str]:
        return set(i.value for i in cls)


@dataclass
class ModelParameter:
    # field names must match header names in the parameters file
    Name: str
    DisplayName: str
    DataType: str
    NumericMin: str
    NumericMax: str
    AllowedStrings: str
    Value: str

    def update(self, newValue):
        self._validate(newValue)
        self.Value = newValue

    def _is_value_special_input(self) -> bool:
        return self.Value in self._split_allowed_strings()

    def _split_allowed_strings(self) -> List[str]:
        return [a.strip() for a in self.AllowedStrings.split(',')]

    def _validate(self, value: str = None):

        if value is None:
            value = self.Value

        # Checks on string inputs
        if self.DataType == ParameterTypes.STRING.value:
            if self.AllowedStrings and not self._is_value_special_input():
                raise ValueError(f"Parameter '{self.Name}' default value is not in list of allowed values.")
            if self.NumericMin or self.NumericMax:
                raise ValueError(f"Parameter '{self.Name}' cannot specify min/max numeric values - it is a string.")

        # Checks on numeric types
        elif self.DataType == ParameterTypes.NUMERIC.value:
            try:
                if self.NumericMin: 
                    float(self.NumericMin)
                if self.NumericMax:
                    float(self.NumericMax)
            except ValueError:
                raise ValueError(f"Numeric parameter '{self.Name}' has invalid min and/or max values." +
                    " Min/Max values should be numeric or blank.")

            try:
                float(value)
            except ValueError:
                if not self._is_value_special_input():
                    raise ValueError(f"Parameter '{self.Name}' default value is invalid.")
                    
            if not self._is_value_special_input():
                if self.NumericMin and float(value) < float(self.NumericMin):
                    raise ValueError(f"Parameter '{self.Name}' value ({value}) is below minimum ({self.NumericMin}).")
                if self.NumericMax and float(value) > float(self.NumericMax):
                    raise ValueError(f"Parameter '{self.Name}' value ({value}) is above maximum ({self.NumericMax}).")

        # Checks on boolean types
        elif self.DataType == ParameterTypes.BOOLEAN.value:
            if self.NumericMin or self.NumericMax:
                raise ValueError(f"Parameter '{self.Name}' cannot specify min/max numeric values - it is a bool.")
            try:
                bool(value)
            except ValueError:
                raise ValueError(f"Parameter '{self.Name}' is not boolean: {value}")

        else:
            raise ValueError(f"Parameter '{self.Name}' has unrecognised type: {self.DataType}.")

def read_parameters_file(filepath: str = default_filepath) -> List[ModelParameter]:
    with open(filepath, 'r') as csvfile:
        rows = [r for r in csv.reader(csvfile, delimiter=',', quotechar="\"")]

    header_row, param_rows = rows[0], rows[1:]
    params = [ModelParameter(**dict(zip(header_row, r))) for r in param_rows]

    for p in params:
        p._validate()

    return params

class RunConfig:
    required_params = set(p.Name for p in read_parameters_file())

    @classmethod
    def default(cls):
        config = cls(read_parameters_file())
        config._validate()
        return config

    def __init__(self, modelParameters: List[ModelParameter]):
        self._params: Dict[str, ModelParameter] = {p.Name: p for p in modelParameters}

    @property
    def params(self) -> Dict[str, ModelParameter]:
        self._validate()
        return self._params

    @property
    def values(self) -> Dict[str, str]:
        self._validate()
        return {k: p.Value for k, p in self.params.items()}

    def update(self, name: str, value: str):
        self._params[name].update(value)

    def generate_suite(self, **kwargs):
        default_config = self.default()

        max_len = max(len(o) for o in kwargs.values() if type(o) is list)
        for v in kwargs.values():
            if type(v) is list and len(v) != max_len:
                raise ValueError("All inputs must be of equal length.")

        base_configs = (deepcopy(default_config) for _ in range(max_len))
        for i, config in enumerate(base_configs):
            for name, value in kwargs.items():
                v = value[i] if type(value) is list else value
                config.update(name, v)
            yield config
        
    def _validate(self):
        for p in self._params.values():
            p._validate()

        if set(self._params.keys()) != self.required_params:
            raise ValueError(f"Run Config does not have the correct parameters")
