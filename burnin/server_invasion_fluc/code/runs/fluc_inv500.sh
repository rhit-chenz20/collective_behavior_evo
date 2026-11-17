# Head-to-head under a fluctuating optimum: 500 social invaders (half of N = 1000) put into a
# psi = 0 resident population (fluc_burnin) at the first optimum flip at or after tick 10000.
# Invaders' psi 0.8 / 0.5 / 0.2, plus 0.0 = a neutral control (fixation probability 50%).
# Shares invasion_fluc with fluc_inv10.
RUN_KIND=invasion
RUN_DIR=invasion_fluc
BURNIN_RUN=fluc_burnin
N_POP_REPS=20               # fluc_burnin's N_REPS_FULL: rep r starts from pop ((r - 1) % 20) + 1
END_TIME=12000
INVASION_IND=500
ID_OFFSET=2000000           # sim IDs 2000002 - 2005761
SZS=(9.4 28.2 47.0 65.8 84.6 103.4)
GAPS=(20 50 80 110)
NS=(8)
PSIS=(0.8 0.5 0.2 0.0)      # 0.8 first, so the test (first 3 combinations) runs copying under every reg
REGS=(ave ext fit)
N_REPS_FULL=20              # for now: these run the longest
# few sims are absorbed early with half the population social: ~6-10 min per sim, so 20 reps
# is ~2-3.5 h, up to ~7 h on slow nodes
REPS_PER_TASK=20
TIME_LIMIT_H_TEST=2
TIME_LIMIT_H_PILOT=6
TIME_LIMIT_H_FULL=24
