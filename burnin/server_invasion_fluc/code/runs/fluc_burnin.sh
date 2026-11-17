# Resident populations for the fluctuating-optimum invasions: psi = 0 (nobody copies, so the rule
# does not matter: reg ave is used), from scratch to tick 9949 under each fluctuation regime of
# fluc_env_adjusted_shift_size (the optimum flips from tick 50 on, as there). One pop file per
# (fluc_opt, fluc_interval, rep) -> fluc_env_pop/pop/fluc_sz_<sz>_gap_<gap>_<rep>.pop (+ .L).
# Run this in full before any invasion run (their tests included).
RUN_KIND=burnin
RUN_DIR=fluc_env_pop
END_TIME=12000              # not used by a burnin (it stops at 9949); kept for the checks
ID_OFFSET=0                 # sim IDs 2 - 481
SZS=(9.4 28.2 47.0 65.8 84.6 103.4)   # 1x 3x 5x 7x 9x 11x 9.4, as in fluc_env_adjusted_shift_size
GAPS=(20 50 80 110)
NS=(8)
PSIS=(0.0)
REGS=(ave)
N_REPS_FULL=20              # one resident population per invasion rep (fluc_inv* run 20 reps)
# ~31 min per sim from scratch (fluc_env_adjusted_shift_size, to tick 11000; 9949 is ~10% less);
# 10 sims took 4.5-10.4 h per task depending on the node
REPS_PER_TASK=10
TIME_LIMIT_H_TEST=4
TIME_LIMIT_H_PILOT=8
TIME_LIMIT_H_FULL=24
