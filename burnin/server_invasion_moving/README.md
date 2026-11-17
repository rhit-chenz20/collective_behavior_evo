# invasion_moving: can social learners invade under a moving optimum?

The steady_tracking model (n = 8, the optimum moves by `shift_speed` every tick from tick
10000), with an invasion: at tick 11000, after 1000 ticks of movement, a number of random
individuals in a psi = 0 population become **social**. The same package also runs the
invasion under **stabilizing selection** (`shift_speed` 0.0: the optimum stays at 0, as in the
burnin; invasion at tick 10000). Only social individuals copy (rule ave / fit / ext, weight psi); the others
keep phenotype = genotype. The copied value comes from the whole group, as in the other
datasets. Social status is inherited: two social parents → social, one → social with
probability 1/2, none → asocial.

A sim stops as soon as the social individuals are **lost** (0) or **fixed** (all 1000), or
2000 ticks after the invasion (tick 13000 moving, 12000 stabilizing).

## Runs

| run | invaders | grid | reps | tasks | sim IDs |
|---|---|---|---|---|---|
| `moving_inv10` | 10 (rare) | 9 speeds 1.1-3.5 x invaders' psi 0.8/0.5/0.2/0.0 x 3 reg = 108 | 200, 25 per task | 864 | 2 - 21601 |
| `moving_inv500` | 500 (half) | same 108 | 60, 20 per task | 324 | 1000002 - 1006481 |
| `stab_inv10` | 10 (rare) | speed 0.0 x 4 psi x 3 reg = 12 | 200, 25 per task | 96 | 2000002 - 2002401 |
| `stab_inv500` | 500 (half) | same 12 | 60, 20 per task | 36 | 3000002 - 3000721 |

The moving runs write to `/fs/cbsubscb10/storage/zc524/invasion_moving`, the stabilizing ones to
`invasion_stabilizing`. In the stabilizing runs the metadata's moving-optimum fitness and
steady-lag records are meaningless. **psi 0.0 is a neutral control**:
tagged invaders that do not copy, so their fixation probability is the drift-only baseline
(1% for 10 invaders, 50% for 500). The invasion tick, number of invaders, reps and speeds are
settings in `code/runs/<run>.sh`.

Residents start from the psi = 0 constant-environment burnins (`pop/` links to
`../server_steady_tracking/pop`: n = 8, 3 regs, reps 1-60). Rep r uses burnin
((r - 1) % 60) + 1 with its own random seed, so reps 61-200 reuse burnins.

## Output

Same layout as steady_tracking, with the number of invaders as a partition:

    parquet/<type>/shift_speed=../invasion_ind=../reg=../n=8/part-psi_<psi>-reps_<a>-<b>-0.parquet

Types: `data`, `metadata`, `phenotype`, `genotype`, `group_tag`, `allele` and **`social`**
(per individual: 1 = social, 0 = asocial, at the same ticks as phenotype). Per-individual records
are dense around the invasion, then every 20 ticks, plus the stop tick.

- `data` (every tick): the steady_tracking columns plus `num_social_ind` and social / asocial
  means of z, g and fitness (`social_z_mean`, `asocial_fitness_mean`, ...; NaN when a class is empty).
- `metadata`: `invasion_tick`, `invasion_ind`, **`outcome`** (lost / fixed / segregating),
  `absorb_tick` (-1 if neither), `final_social`, `stop_tick`, and the steady_tracking fitness and
  steady-lag records up to the stop.

Load with `analysis/load_data/server/load_data.R -d invasion_moving --invasion-ind 10 ...`
(or `-d invasion_stabilizing`).

## Submit

From cbsulogin (the package must be under /home):

    cd ~/collective_behavior_evo_restructure/burnin/server_invasion_moving/code
    bash submit_slurm.sh moving_inv500 test            # first 3 combinations (psi 0.8, one per reg) x 2 reps
    bash submit_slurm.sh moving_inv500 test --status
    bash submit_slurm.sh moving_inv500 full            # after the test finished with the same code
    bash submit_slurm.sh moving_inv10 test
    bash submit_slurm.sh moving_inv10 full
    bash submit_slurm.sh stab_inv500 test        # then full; likewise stab_inv10

Everything else (frozen code snapshots, `--status`, `--unfinished`, email, logs) works as in
`../server_steady_tracking/README.md`. The pre-submit time check projects from the test's
slowest sim; with 10 invaders the test sims may stop early, so that projection is optimistic
(the time limit allows for every sim running to the end).
