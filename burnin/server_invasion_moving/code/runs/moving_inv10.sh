# Rare invasion: 10 social invaders (1% of N = 1000) into a psi = 0 population whose optimum
# starts moving at tick 10000; the invaders arrive at tick 11000, after 1000 ticks of movement.
# Invaders' psi 0.2 / 0.5 / 0.8, plus 0.0 = a neutral control (tagged, but they do not copy): its fixation probability is the
# drift-only baseline, 10 / 1000 = 1%. Most rare invasions are lost to drift within a few
# dozen ticks (and stop there), so many reps are cheap; 200 reps give a neutral expectation of
# ~2 fixations per combination. Shares invasion_moving with moving_inv500 (the combination
# names and partitions carry the number of invaders, so the two never collide).
RUN_DIR=invasion_moving
END_TIME=13000              # at most 2000 ticks after the invasion
INVASION_IND=10
INVASION_TICK=11000         # 1000 ticks after the optimum starts moving, so the lag has settled
ID_OFFSET=0                 # sim IDs 2 - 21601
SPEEDS=(1.1 1.4 1.7 2.0 2.3 2.6 2.9 3.2 3.5)
NS=(8)
PSIS=(0.8 0.5 0.2 0.0)       # 0.8 first, so the test (first 3 combinations) runs copying under every reg
REGS=(ave ext fit)
N_REPS_FULL=200
# every sim runs 1000 ticks before the invasion (~320 s, up to ~870 s on slow nodes, from the
# steady_tracking times); one that keeps its invaders runs 3000 ticks in all (~960 s, up to
# ~2600 s). 25 reps is at most ~18 h, much less when most invasions are lost
REPS_PER_TASK=25
TIME_LIMIT_H_TEST=2
TIME_LIMIT_H_PILOT=4
TIME_LIMIT_H_FULL=24
