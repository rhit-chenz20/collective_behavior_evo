# Head-to-head under stabilizing selection: 500 social invaders (half of N = 1000) into a psi = 0
# population at equilibrium (shift_speed 0.0: the optimum stays at 0; see stab_inv10.sh).
# Invaders' psi 0.8 / 0.5 / 0.2, plus 0.0 = a neutral control (fixation probability 50%).
# Shares invasion_stabilizing with stab_inv10.
RUN_DIR=invasion_stabilizing
END_TIME=12000              # at most 2000 ticks after the invasion
INVASION_IND=500
INVASION_TICK=10000
ID_OFFSET=3000000           # sim IDs 3000002 - 3000721
SPEEDS=(0.0)
NS=(8)
PSIS=(0.8 0.5 0.2 0.0)       # 0.8 first, so the test (first 3 combinations) runs copying under every reg
REGS=(ave ext fit)
N_REPS_FULL=60
# ~640 s per sim, up to ~1740 s: 20 reps is ~3.5 h, ~10 h at worst
REPS_PER_TASK=20
TIME_LIMIT_H_TEST=2
TIME_LIMIT_H_PILOT=6
TIME_LIMIT_H_FULL=24
