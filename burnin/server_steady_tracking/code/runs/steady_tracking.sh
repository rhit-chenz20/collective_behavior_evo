# The original steady_tracking run (finished 2026-11, data in /fs/cbsubscb10/storage/zc524/steady_tracking).
# Kept as the record of its grid; submitting it again would rerun all 108 tasks into the same
# folder. Its frozen code is in that folder's code/.
RUN_DIR=steady_tracking
END_TIME=12000              # 2000 ticks after the shift
ID_OFFSET=0                 # sim IDs 2 - 6481, as in the finished run
SPEEDS=(0.5 0.8 1.1 1.4 1.7 2.0 2.3 2.6 2.9)
NS=(8)
PSIS=(0.0 0.2 0.5 0.8)
REGS=(ave ext fit)
N_REPS_FULL=60
REPS_PER_TASK=60
TIME_LIMIT_H_TEST=2
TIME_LIMIT_H_PILOT=6
TIME_LIMIT_H_FULL=24
