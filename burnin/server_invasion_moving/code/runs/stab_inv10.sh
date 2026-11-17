# Rare invasion under stabilizing selection: 10 social invaders (1% of N = 1000) into a psi = 0
# population at equilibrium. Same model and burnins as the moving runs with shift_speed 0.0: the
# optimum stays at 0, as in the burnin (the moving-optimum fitness / steady-lag records in the
# metadata are then meaningless). Invaders' psi 0.8 / 0.5 / 0.2, plus 0.0 = a neutral control
# (fixation probability 1%). Shares invasion_stabilizing with stab_inv500.
RUN_DIR=invasion_stabilizing
END_TIME=12000              # at most 2000 ticks after the invasion
INVASION_IND=10
INVASION_TICK=10000         # nothing changes over time, so no wait before the invasion
ID_OFFSET=2000000           # sim IDs 2000002 - 2002401
SPEEDS=(0.0)
NS=(8)
PSIS=(0.8 0.5 0.2 0.0)       # 0.8 first, so the test (first 3 combinations) runs copying under every reg
REGS=(ave ext fit)
N_REPS_FULL=200
# at most 2000 ticks per sim (~640 s, up to ~1740 s); 25 reps is at most ~12 h, much less
# when most invasions are lost
REPS_PER_TASK=25
TIME_LIMIT_H_TEST=2
TIME_LIMIT_H_PILOT=4
TIME_LIMIT_H_FULL=24
