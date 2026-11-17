# Rare invasion under a fluctuating optimum: 10 social invaders (1% of N = 1000) put into a
# psi = 0 resident population (fluc_burnin) at the first optimum flip at or after tick 10000.
# Invaders' psi 0.8 / 0.5 / 0.2, plus 0.0 = a neutral control (tagged, but they do not copy):
# its fixation probability is the drift-only baseline, 1%. Most rare invasions are lost within
# a few dozen ticks and stop there. With 20 reps a neutral fixation is rarely seen; raise
# N_REPS_FULL (and fluc_burnin's) to estimate fixation rates. Shares invasion_fluc with fluc_inv500.
RUN_KIND=invasion
RUN_DIR=invasion_fluc
BURNIN_RUN=fluc_burnin
N_POP_REPS=20               # fluc_burnin's N_REPS_FULL: rep r starts from pop ((r - 1) % 20) + 1
END_TIME=12000              # at most ~2000 ticks after the invasion
INVASION_IND=10
ID_OFFSET=1000000           # sim IDs 1000002 - 1005761
SZS=(9.4 28.2 47.0 65.8 84.6 103.4)
GAPS=(20 50 80 110)
NS=(8)
PSIS=(0.8 0.5 0.2 0.0)      # 0.8 first, so the test (first 3 combinations) runs copying under every reg
REGS=(ave ext fit)
N_REPS_FULL=20              # for now: these run the longest
# a sim that keeps its invaders runs ~2050 ticks from the burnin (~6-10 min; per-individual
# records every 10 ticks); most lose them early. 20 reps is at most ~3.5 h
REPS_PER_TASK=20
TIME_LIMIT_H_TEST=2
TIME_LIMIT_H_PILOT=4
TIME_LIMIT_H_FULL=24
