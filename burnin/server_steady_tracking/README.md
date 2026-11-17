# steady_tracking: which shift speeds can the population track?

(The description below is of the original run, `code/runs/steady_tracking.sh`; see Runs for the others.)

Moving optimum, n = 8, psi 0, 0.2, 0.5, 0.8, rules ave / ext / fit,
speeds 0.5 to 2.9 by 0.3 (9 speeds), 60 reps each:
108 Slurm array tasks (one per speed x psi x rule) x 60 reps = **6480 sims**.

The optimum moves by `shift_speed` every tick from tick 10000, and **every sim runs to tick
12000** (2000 ticks after the shift). Failing to track is only recorded: a tick fails the
fitness check if mean fitness < 0.1 x the fitness at the optimum (a lag of about 96). The
metadata has

| column | meaning |
|---|---|
| `failed` | T if any tick failed the check |
| `fail_tick`, `fail_lag` | first failing tick and the lag then (-1 and NaN if never) |
| `n_fail_ticks` | number of failing ticks (a brief dip vs. a lasting failure) |
| `min_fit_ratio`, `min_fit_tick` | lowest mean fitness / fitness at the optimum, and when |
| `fit_ratio`, `end_lag` | the same ratio and the lag at tick 12000 |

A steady lag is also only recorded: in test runs, populations that later
failed first held a near-constant lag of 80-90 for a few hundred ticks. Every 100 ticks from
200 ticks after the shift, the mean lag (optimum - mean phenotype) over the last 100 ticks is
compared with the mean over the 100 ticks before. A check passes if they differ by less than
`steady_tol = steady_eps * shift_speed * steady_window` (0.05 x speed x 100: the population
advances at 95% or more of the optimum's speed). `steady_tick` is the first tick at which 2
checks in a row passed (-1 if never) and `steady_lag` the window mean then. The data log has
the lag every tick (`z_opt`), so other rules can be applied afterwards.

The parameters are `fail_fit_ratio`, `steady_window`, `steady_eps`, `steady_n` and `end_time`
at the top of `code/base_script.slim`. The metadata also has `lag_prev` and `lag_last` (the
last two window means). Per-individual and allele records are dense around the shift, then
every 20 ticks; the final tick is always recorded.

## Runs

The pipeline runs any of the runs defined in `code/runs/<run>.sh`. A run file sets the
combinations (a grid of speeds x n x psi x reg, or a hand-picked `COMBOS` list), the end tick
(`END_TIME`, passed to SLiM as `end_time`), replicates (`N_REPS_FULL`, split over tasks of
`REPS_PER_TASK`), an `ID_OFFSET` so sim IDs never repeat across runs, the time limits and the
output folder (`RUN_DIR`). Several runs can share one output folder as long as their
combinations differ: each task's `summary/<combo>/_RUN` records the run it belongs to, and both
the submit and the job refuse to touch another run's task (tasks from before runs existed
belong to the folder's own run). Each run's frozen code goes to `code/` for the run the folder
is named after and to `code_<run>/` for a run added to it. Test and pilot always go to
`test_<run>` / `pilot_<run>`.

| run | what | combinations | reps | to tick | output |
|---|---|---|---|---|---|
| `steady_tracking` | the finished original run (record of its grid) | 9 speeds 0.5-2.9 x 4 psi x 3 reg = 108 | 60 | 12000 | `steady_tracking` |
| `long_borderline` | longer runs of the combinations not settled at 12000, plus 2 references | 11 (list) | 20, 10 per task | 15000 | `steady_tracking_long` |
| `fast_speeds` | speeds above 2.9, added to the original dataset | 7 speeds 3.2-5.0 x 4 psi x 3 reg = 84 | 20 | 12000 | `steady_tracking` (code in `code_fast_speeds/`) |

`long_borderline` has its own folder because its combinations are the same as the original's
with longer sims: in one dataset the analyses would mix sims of two lengths. Only add a run to
an existing folder if it uses the same model and `END_TIME` and new combinations.

To add a run, copy one of these files, change it, and give it a `RUN_DIR` and an
`ID_OFFSET` that does not overlap the others (the original uses IDs 2-6481; long_borderline
1000002 and up; fast_speeds 2000002 and up). Speeds and psi must be written as in the pop file
names (1.1 not 1, 0.0 not 0); only n = 8 and psi 0.0 / 0.2 / 0.5 / 0.8 have burnins, reps 1-60.
With the original grid and `ID_OFFSET=0` the scripts give exactly the original run's sim IDs.

## Contents

    code/   scripts; submit from here
        submit_slurm.sh     checks the setup, freezes the code, makes ave/ext/fit.slim, submits the array job
        slurm_job.sh        one task: its sims -> parquet -> archive -> lab storage
        pipeline_config.sh  paths, task layout, nodes to avoid; reads runs/<run>.sh
        runs/               one file per run (grid or combination list, end tick, reps, output folder)
        base_script.slim    the model; make_model.py turns it into ave/ext/fit.slim
        convert_parquet.py  tsv -> parquet
    pop/    the 720 burnin files (cons_n_8_psi_<psi>_reg_<reg>_<rep>.pop, reps 1-60)

Uses what is already on the server: SLiM at ~/bin/slim5.0 and the simdata conda env
(python with polars, zstd) in ~/miniconda3/envs/simdata.

## Run

The package must be on a /home path (the compute nodes cannot see /local/workdir). From cbsulogin
or cbsuhein:

    cd ~/collective_behavior_evo_restructure/burnin/server_steady_tracking/code
    bash submit_slurm.sh long_borderline test    # the run's first 3 combinations (one per reg) x 2 reps
    bash submit_slurm.sh long_borderline test --status   # check, then look at logs/ in test_long_borderline
    bash submit_slurm.sh long_borderline pilot   # optional: every combination x 3 reps
    bash submit_slurm.sh long_borderline full    # every combination x N_REPS_FULL reps

The run and the mode are required and everything else comes from `code/pipeline_config.sh` and
`code/runs/<run>.sh`; variables set in your shell are ignored. pilot and full are refused until
that run's test has finished with the same code (only time limits and EXCLUDE_NODES may change),
and they check that the test's slowest sim x reps per task fits their time limit. They ask you
to type `yes` before submitting. Each submit freezes a copy of the code in
`~/slurm_code/<dir>_<time>/` (and `$OUT/code`) and the jobs run that copy, so editing `code/`
while tasks are queued is safe.

Sim times in the original run (2000 ticks after the shift): mean 640 s, max 1740 s, about the
same at every speed; nodes differ by more than 2x. long_borderline runs 5000 ticks after the
shift (~2.5x as long), so its reps are split into tasks of 10. Over 24 h, jobs go to long7 only.

Every submit asks Slurm to email zc524@cornell.edu once when the whole array has ended, and
if any task fails (`MAIL_USER` at the top of `code/submit_slurm.sh`).

Progress and reruns:

    squeue -u $USER
    bash submit_slurm.sh <run> full --status       # how many of the run's tasks are done
    bash submit_slurm.sh <run> full --unfinished   # resubmits the tasks without a _DONE that are not
                                                   # still queued or running; refused if the code changed

Output (on lab storage, /fs/cbsubscb10/storage/zc524/<RUN_DIR>):

    parquet/<type>/shift_speed=../reg=../n=8/part-psi_<psi>-0.parquet          (one task per combination)
    parquet/<type>/shift_speed=../reg=../n=8/part-psi_<psi>-reps_<a>-<b>-0.parquet  (reps split over tasks)
    raw/<combo>.tar.zst                 raw tsv + SLiM logs
    summary/<combo>/                    conversion summaries, failed_sims.txt, sim_seconds.tsv, _DONE
    logs/<combo>.<jobid>.log            each task's log (+ .convert.log); submit_<time>.log per submit
    code/, code_<run>/                  the code the tasks ran (code_<run>/: a run added to this folder)

`<combo>` is `sp_<speed>_reg_<reg>_n_8_psi_<psi>`, with `_reps_<a>-<b>` added when reps are split.
The analysis loader reads every folder the same way: `-d steady_tracking` then includes the
fast_speeds speeds once they are done (filter with `--shift-speed`), and
`-d steady_tracking_long` loads the long runs.

Slurm itself has to open its log on /home (`~/slurm_logs`) before a job can mount lab storage;
each job deletes that file once its own log is on storage, so one is only left there if a job
died before it could copy its log.
