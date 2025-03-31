from collections.abc import Mapping
from numbers import Number

import sys
import json

import numpy as np
import scipy.stats as stats

sys.path.append(snakemake.scriptdir + "/../..")
from project_parameters import default_simulation_parameters

def number_or_string(x):
    try:
        result = float(x)
    except ValueError as err:
        if isinstance(x, str):
            result = x
        else:
            raise err
    return result


def instantiate_parameter(key, value):
    if isinstance(value, Mapping):
        instantiated = draw_from_distribution(value)
    else:
        instantiated = value
    try:
        json_friendly_value = number_or_string(instantiated)
    except Exception as err:
        raise ValueError(
            f"Could not instantiate {key} from:\n{type(value)}\n{value}. Got this error:\n{err}"
        )
    return json_friendly_value

params = snakemake.config["simulation-parameters"]

results = dict()
for key, value in params.items():
    instantiated = instantiate_parameter(key, value)
    results[key] = instantiated

# We always want to start every parameter set with a seed.
# If not given explicitly in the desired parameters, we make one right here.
if "seed" not in results:
    results["seed"] = np.random.randint(1, 2 ** 32)

# Project-wide parameters are added now if not overwritten by the specific config file.
for key, value in default_simulation_parameters.items():
    if key not in results:
        results[key] = number_or_string(value)

with open(snakemake.output["params_file"], "w") as f:
    json.dump(results, f, sort_keys=True, indent=2)