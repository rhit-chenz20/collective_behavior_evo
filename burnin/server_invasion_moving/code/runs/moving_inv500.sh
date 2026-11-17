# Head-to-head: 500 social invaders (half of N = 1000) into a psi = 0 population whose optimum
# starts moving at tick 10000; the invaders arrive at tick 11000, after 1000 ticks of movement.
# Invaders' psi 0.8 / 0.5 / 0.2, plus 0.0 = a neutral control (tagged, but they do not copy;
# fixation probability 50%). Shares invasion_moving with moving_inv10.
RUN_DIR=invasion_moving
END_TIME=13000              # at most 2000 ticks after the invasion
INVASION_IND=500
INVASION_TICK=11000         # 1000 ticks after the optimum starts moving, so the lag has settled
ID_OFFSET=1000000           # sim IDs 1000002 - 1006481
SPEEDS=(1.1 1.4 1.7 2.0 2.3 2.6 2.9 3.2 3.5)
NS=(8)
PSIS=(0.8 0.5 0.2 0.0)       # 0.8 first, so the test (first 3 combinations) runs copying under every reg
REGS=(ave ext fit)
N_REPS_FULL=60
# with half the population social, few sims are absorbed early: 3000 ticks per sim, ~960 s on
# average and up to ~2600 s (steady_tracking times x 1.5), so 20 reps is ~5.5 h, ~14.5 h at worst
REPS_PER_TASK=20
TIME_LIMIT_H_TEST=2
TIME_LIMIT_H_PILOT=6
TIME_LIMIT_H_FULL=24
