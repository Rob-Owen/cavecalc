from typing import List, Dict
from copy import deepcopy
from cavecalc.configuration.model_parameters import read_parameters_file, ModelParameter

class RunConfig:
    """
    Encapsulates configuration settings for a single cavecalc model run (i.e. a collection
    of ModelParameters).
    """
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
