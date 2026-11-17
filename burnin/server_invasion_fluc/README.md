# invasion_fluc: can social learners invade under a fluctuating optimum?

The fluc_env_adjusted_shift_size model (n = 8; the optimum alternates between +fluc_opt/2 and
-fluc_opt/2, flipping every fluc_interval ticks from tick 50), with an invasion: at the first
flip at or after tick 10000, a number of random individuals in a psi = 0 population become
**social**. Only social individuals copy (rule ave / fit / ext, weight psi); the others keep
phenotype = genotype. Social status is inherited: two social parents → social, one → social with
probability 1/2, none → asocial. A sim stops as soon as the social individuals are **lost** or
**fixed**, or at tick 12000.

## Two stages

Running every invasion from scratch takes ~31 min per sim. Residents are psi = 0, so the copying rule does not affect them, and one
resident population per (fluc_opt, fluc_interval, rep) serves every rule, psi and invader count:

1. **`fluc_burnin`**: psi = 0 from scratch to tick 9949 for each of the 24 regimes x 20 reps
   (480 sims, ~250 CPU-hours) → `/fs/cbsubscb10/storage/zc524/fluc_env_pop/pop/fluc_sz_<sz>_gap_<gap>_<rep>.pop`
   plus `<...>.pop.L`, the mutable-region length L, which the pop file does not hold. Saved at
   9949 so the invasion runs log from 9950, like the dataset.
2. **`fluc_inv10`**, **`fluc_inv500`**: read a burnin, restore L and the optimum, invade.

| run | kind | grid | reps | tasks | sim IDs |
|---|---|---|---|---|---|
| `fluc_burnin` | burnin | fluc_opt 9.4-103.4 (1x-11x 9.4) x fluc_interval 20/50/80/110 = 24 | 20, 10 per task | 48 | 2 - 481 |
| `fluc_inv10` | 10 invaders | 24 x invaders' psi 0.8/0.5/0.2/0.0 x 3 reg = 288 | 20, one task each | 288 | 1000002 - 1005761 |
| `fluc_inv500` | 500 invaders | same 288 | 20, one task each | 288 | 2000002 - 2005761 |

The invasion runs write to `/fs/cbsubscb10/storage/zc524/invasion_fluc`. psi 0.0 is a neutral
control (tagged invaders that do not copy). Rep r uses burnin r (20 of each). 20 reps is a first
pass: for fixation rates with 10 invaders, raise `N_REPS_FULL` in the invasion run and in
`fluc_burnin` (or reuse burnins: rep r uses burnin ((r - 1) % N_POP_REPS) + 1).

## Output (invasion runs)

    parquet/<type>/fluc_opt=../fluc_interval=../invasion_ind=../reg=../n=8/part-psi_<psi>-reps_<a>-<b>-0.parquet

Types: `data`, `metadata`, `phenotype`, `genotype`, `group_tag`, `allele` and **`social`**
(per individual, 1 = social). Per-individual and allele records: dense around the invasion, then
every 10 ticks (at least 4 per cycle), plus the stop tick. `data` (every tick from 9950) adds
`num_social_ind` and social / asocial means of z, g and fitness; `metadata` adds
`invasion_tick`, `invasion_ind`, **`outcome`** (lost / fixed / segregating), `absorb_tick`,
`final_social`, `stop_tick` and `burnin_L`.

Load with `analysis/load_data/server/load_data.R -d invasion_fluc --invasion-ind 500 ...`.

## Submit, in this order

From cbsulogin (the package must be under /home):

    cd ~/collective_behavior_evo_restructure/burnin/server_invasion_fluc/code
    bash submit_slurm.sh fluc_burnin test      # 3 regimes x 2 reps, ~1 h
    bash submit_slurm.sh fluc_burnin full      # after the test; 48 tasks
    # when fluc_burnin full is done (bash submit_slurm.sh fluc_burnin full --status):
    bash submit_slurm.sh fluc_inv500 test
    bash submit_slurm.sh fluc_inv500 full
    bash submit_slurm.sh fluc_inv10 test
    bash submit_slurm.sh fluc_inv10 full

An invasion submit, its test included, is refused while any resident pop file it needs is
missing. Everything else (frozen code, `--status`, `--unfinished`, email, logs) works as in
`../server_steady_tracking/README.md`.
