# Speeds above the steady_tracking grid (which went up to 2.9), for every psi and reg: where do
# psi 0 / 0.2 and ave finally run away, or do they keep settling at ever larger lags?
# Same model and span as steady_tracking (2000 ticks after the shift), so its results go into
# the same folder and extend that dataset's speed axis to 0.5-5.0 (new shift_speed= folders; no
# task of the original is touched). Its frozen code goes to steady_tracking/code_fast_speeds/,
# its test and pilot to test_fast_speeds / pilot_fast_speeds. Keep END_TIME at 12000 while it
# shares that folder: sims of different lengths in one dataset would be mixed by the analyses.
RUN_DIR=steady_tracking
END_TIME=12000              # 2000 ticks after the shift, as in steady_tracking
ID_OFFSET=2000000           # sim IDs 2000002 and up
SPEEDS=(3.2 3.5 3.8 4.1 4.4 4.7 5.0)
NS=(8)
PSIS=(0.0 0.2 0.5 0.8)
REGS=(ave ext fit)
N_REPS_FULL=20
# ~640 s per sim on average, up to ~1740 s: 20 reps is ~3.5 h on average, ~10 h at worst
REPS_PER_TASK=20
TIME_LIMIT_H_TEST=2
TIME_LIMIT_H_PILOT=6
TIME_LIMIT_H_FULL=24
