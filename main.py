import numpy as np

from cal.Optimizer import *
from language.Eventlog import *

def xp_opt(filename):

    ev = Eventlog(filename)
    pn = ev.discover_pn_inductive()
    opt = Optimizer(ev, pn)

    # L-BFGS-B, TNC, Powell, Nelder-Mead
    # For TNC, use maxfun instead of maxiter.
    m, nw0, solver, maxiter = "KLD", 100, "L-BFGS-B", np.inf
    memoized, derivatives = False, False
    opt.estimate(m, nw0, solver, maxiter, memoized, derivatives)

if __name__ == '__main__':
    filename = "rl_data/BPIC17_offerlog.xes"
    xp_opt(filename)
