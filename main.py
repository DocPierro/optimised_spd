
from cal.estimators import *
from cal.ML_estimator import *
from language.Eventlog import *
from model.Petrinet import *

from pm4py.visualization.transition_system import visualizer as rg_visualizer


def xp_mle(mle, filename, miner):

    nw0, step, maxiter, improvement_threshold = 10, 1, 50, 0.001
    minw, maxw = 0.001, 1
    method = "Powell"
    weights_mle = mle.estimate_emd(method, nw0, step, maxiter, minw, maxw, improvement_threshold, 1)
    #mle.pn.export_pnml(filename.replace("/", "_") + "_emd_" + miner + "_" + str(nw0) + "_" + str(step) + "_" + str(maxiter) + "_" + str(improvement_threshold), weights_mle)

def xp_estimators(mle, filename, miner):

    ev_language = mle.ev.get_language_pm4py()

    # frequency_estimators
    weights_frequency = frequency_estimator(mle.ev, mle.pn)
    pn_language_frequency = mle.unfold_rg(list(weights_frequency.values()),False)[1]
    print("frequency estimator:")
    print("\tloglh: " + str(mle.log_likelihood(ev_language, pn_language_frequency)))
    print("\temd: " + str(mle.emd_norm_lev(ev_language, pn_language_frequency)))
    #mle.pn.export_pnml(filename.replace("/", "_") + "_" + miner + "_frequency", weights_frequency)

    # lhpair_estimators
    weights_lhpair = lhpair_estimators(mle.ev, mle.pn)
    pn_language_lhpair = mle.unfold_rg(list(weights_lhpair.values()), False)[1]
    print("lhpair estimator:")
    print("\tloglh: " + str(mle.log_likelihood(ev_language, pn_language_lhpair)))
    print("\temd: " + str(mle.emd_norm_lev(ev_language, pn_language_lhpair)))
    #mle.pn.export_pnml(filename.replace("/", "_") + "_" + miner + "_lhpair", weights_lhpair)

    # rhpair_estimators
    weights_rhpair = rhpair_estimators(mle.ev, mle.pn)
    pn_language_rhpair = mle.unfold_rg(list(weights_rhpair.values()), False)[1]
    print("rhpair estimator:")
    print("\tloglh: " + str(mle.log_likelihood(ev_language, pn_language_rhpair)))
    print("\temd: " + str(mle.emd_norm_lev(ev_language, pn_language_rhpair)))
    #mle.pn.export_pnml(filename.replace("/", "_") + "_" + miner + "_rhpair", weights_rhpair)

    # pairscale_estimators
    weights_pairscale = pairscale_estimators(mle.ev, mle.pn)
    pn_language_pairscale = mle.unfold_rg(list(weights_pairscale.values()), False)[1]
    print("pairscale estimator:")
    print("\tloglh: " + str(mle.log_likelihood(ev_language, pn_language_pairscale)))
    print("\temd: " + str(mle.emd_norm_lev(ev_language, pn_language_pairscale)))
    #mle.pn.export_pnml(filename.replace("/", "_") + "_" + miner + "_pairscale", weights_pairscale)

    # fork_estimators
    weights_fork = fork_estimators(mle.ev, mle.pn)
    pn_language_fork = mle.unfold_rg(list(weights_fork.values()), False)[1]
    print("fork estimator:")
    print("\tloglh: " + str(mle.log_likelihood(ev_language, pn_language_fork)))
    print("\temd: " + str(mle.emd_norm_lev(ev_language, pn_language_fork)))
    #mle.pn.export_pnml(filename.replace("/", "_") + "_" + miner + "_fork", weights_fork)

if __name__ == '__main__':

    filename = "rl_data/BPIC13_closed.xes"
    miner = "inductive"

    print("##### " + filename + " with " + miner + " #####")

    ev = Eventlog(filename)
    pn = ev.discover_pn_inductive()
    mle = MaxLikelihoodEstimators(ev, pn)

    xp_mle(mle, filename, miner)
    xp_estimators(mle, filename, miner)
    print("")
