# Longer runs for the combinations that were not settled at tick 12000 in steady_tracking:
# the lag looked flat but H was still falling (or the lag was still creeping up), so they may
# still run away. Each runs 5000 ticks after the shift (to tick 15000).
#   ext psi 0.8 speed 1.4 / 1.7      H -4% / -9% over ticks 11000-12000
#   ext psi 0.5 speed 2.3 / 2.6      H -2% / -5%
#   fit psi 0.8 speed 2.0 / 2.3      lag +6 / +15 per 1000 ticks, H -3% / -5%
#   fit psi 0.5 speed 2.9            H -2.5%
#   ave psi 0.8 speed 2.6 / 2.9      lag about 230 / 300 and flat, H -2% / -3%
# plus two settled references (psi 0.0 at speed 2.9, ave and fit) to compare over the same span.
# The first three are one of each reg, so the test (first 3 combinations) runs all three scripts.
RUN_DIR=steady_tracking_long
END_TIME=15000              # 5000 ticks after the shift
ID_OFFSET=1000000           # sim IDs 1000002 and up
COMBOS=(
    "1.7 8 0.8 ext"
    "2.3 8 0.8 fit"
    "2.9 8 0.8 ave"
    "1.4 8 0.8 ext"
    "2.3 8 0.5 ext"
    "2.6 8 0.5 ext"
    "2.0 8 0.8 fit"
    "2.9 8 0.5 fit"
    "2.6 8 0.8 ave"
    "2.9 8 0.0 ave"
    "2.9 8 0.0 fit"
)
N_REPS_FULL=20
# 2.5x the ticks of steady_tracking: ~1600 s per sim on average, up to ~4400 s on slow nodes,
# so 10 reps per task is ~4.5 h on average and ~12 h at worst
REPS_PER_TASK=10
TIME_LIMIT_H_TEST=4
TIME_LIMIT_H_PILOT=8
TIME_LIMIT_H_FULL=24
