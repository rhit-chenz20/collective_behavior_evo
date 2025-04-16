
from random import choice

if 'slim' not in config:
    raise Exception("Need to provide a path to SLiM executable. Please invoke snakemake with --config slim=PATH_TO_SLIM.")

if 'use_subdirectory' not in config:
    config['use_subdirectory'] = False

def get_outdir():
    outdir = Path("simulations")
    if config['use_subdirectory']:
        vowels = 'aeiou'
        consonants = 'bcdfghjklmnprstvz'
        subdir = choice(consonants) + choice(vowels) + choice(consonants) + choice(vowels)
        return outdir / subdir
    else:
        return outdir

OUTDIR_ID = get_outdir()
OUTDIR = Path('output')/OUTDIR_ID
SIM_IDS = [
    str(i + 1).zfill(len(str(config['simulations'])))
    for i in range(config['simulations'])
]


rule all:
    input:
        OUTDIR/".slim_all.done"

rule instantiate_parameters:
    output:
        params_file = OUTDIR/"{sim_id}_simulation-params.json"
    script: "slim/instantiate-simulation-parameters.py"

rule slim_script:
    input:
        params_file = OUTDIR/"{sim_id}_simulation-params.json",
        slim_templates = expand(
            'slim/template/{model}.slim',
            model=['1_trait']
        )
    output:
        slim_script = OUTDIR/"{sim_id}_script.slim"
    params:
        outdir = OUTDIR
    # conda: "envs/simulate.yaml"
    script: "slim/make_slim_script.py"

rule slim:
    input:
        OUTDIR/"{sim_id}_script.slim"
    output:
        touch(OUTDIR/".{sim_id}_slim.done")
    params:
        slim = config['slim']
    log: 'logs/' + str(OUTDIR_ID) + '/{sim_id}_slim.log'
    shell:
        "{params.slim} {input} &> {log}"

checkpoint slim_has_run:
    input:
        expand(OUTDIR/".{sim_id}_slim.done", sim_id=SIM_IDS)
    output:
        touch(OUTDIR/".slim_all.done")

