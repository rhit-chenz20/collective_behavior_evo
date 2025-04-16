from pathlib import Path
from collections import defaultdict

import json

import numpy as np

with open(snakemake.input["params_file"], "r") as f:
    specified_params = json.load(f)

# Missing parameter values get a placeholder -1.
params = defaultdict(lambda: -1)
params.update(specified_params)

replacements = {
    # SLiM configuration
    "INITIAL_SEED": params["seed"],
    "LAST_GENERATION": int(params["last_generation"]),
    "LOG_INTERVAL": int(params["log_interval"]),
    # Biology
    "NUM_OF_QTLS": int(params["num_of_qtl"]),
    "BUFFER_SIZE": int(params["buffer_size"]),
    "POP_SIZE": int(params["pop_size"]),
    "NEIGHBOR_SIZE": int(params["neighbor_size"]),
    "PSI11": float(params["psi11"]),
    "OPTIMUM1": float(params["optimum"]),
    "SELECTION_GENERATION": int(params["selection_generation"]),
    # Files
    "CURRENT_WORKDIR": Path(".").resolve(),
    "FREQUENCY_OUTPUT_FILE": str(
        Path(snakemake.params["outdir"])
        / (snakemake.wildcards["sim_id"] + "_slim.tsv")
    ),
    "METRICS_OUTPUT_FILE": str(
        Path(snakemake.params["outdir"])
        / (snakemake.wildcards["sim_id"] + "_metrics.txt")
    ),
}

def find_template(name):
    template = None
    for template_option in snakemake.input["slim_templates"]:
        if Path(template_option).stem == name:
            template = template_option
            break
    if template is None:
        raise FileNotFoundError(
            f"Could not find template {name}. Template options were:\n{snakemake.input['slim_templates']}"
        )
    return template


def load_template(template_file, replacements):
    with open(template_file, "r") as f:
        text = "".join(line for line in f)
    for key, replacer in replacements.items():
        text = text.replace(key, str(replacer))
    return text


template = find_template(params["model"])
script_text = load_template(template, replacements)
with open(snakemake.output["slim_script"], "w") as f:
    f.write(script_text)