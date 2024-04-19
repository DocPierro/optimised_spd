import numpy
import time
from copy import deepcopy
from bigtree import *
from operator import *
from queue import Queue
from scipy.optimize import minimize
from Levenshtein import distance as levenshtein_distance
from pm4py.algo.evaluation.earth_mover_distance.variants import pyemd


class MaxLikelihoodEstimators:

    def __init__(self, ev, pn):
        self.ev, self.pn = ev, pn
        self.rg = pn.get_rg()
        self.traces = None

    def log_likelihood(self, l1, l2):
        lh = 0
        for tt in l2:
            lh += l1[tt] * numpy.log(l2[tt])
        return -lh

    def emd_norm_lev(self, l1, l2):

        nc1, nc2 = 0, 0
        words = set(l1) & set(l2)
        for w in words:
            if w in l1:
                nc1 = nc1 + l1[w]
            if w in l2:
                nc2 = nc2 + l2[w]

        hist1 = numpy.zeros(len(words))
        hist2 = numpy.zeros(len(words))
        distance = numpy.zeros((len(words), len(words)))
        for x, xw in enumerate(words):
            for y, yw in enumerate(words):
                if xw in l1:
                    hist1[x] = l1[xw] / nc1
                if yw in l2:
                    hist2[y] = l2[yw] / nc2
                if len(xw) == 0 or len(yw) == 0:
                    distance[x, y] = 1.
                else:
                    maxl = float(max(len(xw), len(yw)))
                    distance[x, y] = float(levenshtein_distance(xw, yw)) / maxl

        return pyemd.emd(hist1, hist2, distance)

    def compute_rg(self, w):
        weights = {transition.get_name(): w[t] for t, transition in enumerate(self.pn.get_transitions())}
        for transition in self.rg.transitions:
            denom = 0
            for enabled_transition in transition.from_state.outgoing:
                denom += weights[enabled_transition.name]
            # use property data to memorize prob
            transition.data = weights[transition.name] / denom

    def unfold_rg(self, w, p, only_loglh):

        self.compute_rg(w)
        ev_language = self.ev.get_language_pm4py()
        max_length = max([len(trace.get_seq()) for trace in self.ev.get_language()])
        prefixes = self.ev.build_prefix_set()

        # find initial state in the reachability graph
        init_state = None
        for state in self.rg.states:
            if state.name == 'start1':
                init_state = state

        q = Queue(0)
        q.put((init_state, [], 0))

        triples = {}
        triples[(init_state, tuple([]), 0)] = 1

        traces = {}  # all traces that arrives to sink
        while not q.empty():
            state, trace, level = q.get()
            prob = triples[(state, tuple(trace), level)]
            if len(state.outgoing) > 0:
                for transition in state.outgoing:
                    new_state = transition.to_state
                    new_trace = trace.copy()
                    label = self.pn.get_tr2lab()[transition.name]
                    if label is not None:
                        new_trace.append(label)
                    if len(new_trace) <= max_length:
                        if tuple(new_trace) in prefixes:
                            new_prob = prob * transition.data
                            key = (new_state, tuple(new_trace), level + 1)
                            if key not in triples:
                                triples[key] = new_prob
                                q.put((new_state, new_trace, level + 1))
                            else:
                                triples[key] = triples[key] + new_prob
            else:
                tt = tuple(trace)
                if tt in ev_language:
                    if tt not in traces:
                        traces[tt] = prob
                    else:
                        traces[tt] = traces[tt] + prob

        param = 0
        if p == "loglh":
            param = self.log_likelihood(ev_language, traces)
        elif p == "emd":
            param = self.emd_norm_lev(ev_language, traces)

        if only_loglh:
            return param
        else:
            return param, traces

    def estimate_loglh(self, method, nw0, step, maxiter, minw, maxw, improvement_threshold, v=0):

        print("Minimize lohlh | method: " + method +
              " | nw0: " + str(nw0) + " | step: " + str(step) +
              " | maxiter: " + str(maxiter) +
              " | minw: " + str(minw) + " | maxw: " + str(maxw) +
              " | improvement_threshold: " + str(improvement_threshold))

        start_time = time.time()
        ev_language = self.ev.get_language_pm4py()

        best_w_start, best_loglh_start, best_pn_language_start = [], numpy.inf, {}
        for _ in range(nw0):
            w_start = (maxw - minw) * numpy.random.rand(len(self.pn.get_transitions())) + minw
            loglh_start, pn_language_start = self.unfold_rg(w_start, "loglh", False)
            if loglh_start < best_loglh_start:
                best_w_start, best_loglh_start, best_pn_language_start = w_start, loglh_start, pn_language_start
        w_start, loglh_start, pn_language_start = best_w_start, best_loglh_start, best_pn_language_start
        emd_start = self.emd_norm_lev(ev_language, pn_language_start)
        print("\tBest starting loglh value: " + str(loglh_start))
        print("\tBest starting emd value: " + str(emd_start))
        print("\t#############################################")

        bnds = ()
        for i in range(len(self.pn.get_transitions())):
            bnds += ((minw, maxw),)

        w_before, loglh_before = w_start, loglh_start
        w_after, loglh_after = [], numpy.inf
        for i in range(maxiter//step):
            w_after = minimize(self.unfold_rg, w_before, ("loglh", True),
                            bounds=bnds, method=method,
                            options={'maxiter': step, 'disp': False}).x
            loglh_after, pn_language_after = self.unfold_rg(w_after, "loglh", False)
            improvement = abs(loglh_after - loglh_before)/loglh_after
            if v >= 1:
                emd_after = self.emd_norm_lev(ev_language, pn_language_after)
                print("\tloglh value at step " + str((i+1)*step) + " : " + str(loglh_after) + " (" + str(improvement) + ")")
                print("\temd value at step " + str((i+1)*step) + " : " + str(emd_after))

            if improvement <= improvement_threshold:
                break
            w_before = w_after
            loglh_before = loglh_after

        print("\t#############################################")
        emd_after = self.emd_norm_lev(ev_language, pn_language_after)
        print("\tloglh value with optimized values: " + str(loglh_after))
        print("\temd value with optimized values: " + str(emd_after))
        print("\tExecution time: " + str(time.time() - start_time) + " seconds")
        return {transition.name: w_after[i] for i, transition in enumerate(self.pn.get_transitions())}

    def estimate_emd(self, method, nw0, step, maxiter, minw, maxw, improvement_threshold, v=0):

        print("Minimize emd | method: " + method +
              " | nw0: " + str(nw0) + " | step: " + str(step) +
              " | maxiter: " + str(maxiter) +
              " | minw: " + str(minw) + " | maxw: " + str(maxw) +
              " | improvement_threshold: " + str(improvement_threshold))

        start_time = time.time()
        ev_language = self.ev.get_language_pm4py()

        best_w_start, best_emd_start, best_pn_language_start = [], numpy.inf, {}
        for _ in range(nw0):
            w_start = (maxw - minw) * numpy.random.rand(len(self.pn.get_transitions())) + minw
            emd_start, pn_language_start = self.unfold_rg(w_start, "emd", False)
            if emd_start < best_emd_start:
                best_w_start, best_emd_start, best_pn_language_start = w_start, emd_start, pn_language_start
        w_start, emd_start, pn_language_start = best_w_start, best_emd_start, best_pn_language_start
        loglh_start = self.log_likelihood(ev_language, pn_language_start)
        print("\tBest starting emd value: " + str(emd_start))
        print("\tBest starting loglh value: " + str(loglh_start))
        print("\t#############################################")

        bnds = ()
        for i in range(len(self.pn.get_transitions())):
            bnds += ((minw, maxw),)

        w_before, emd_before = w_start, emd_start
        w_after, emd_after = [], numpy.inf
        for i in range(maxiter//step):
            w_after = minimize(self.unfold_rg, w_before, ("emd", True),
                            bounds=bnds, method=method,
                            options={'maxiter': step, 'disp': False}).x
            emd_after, pn_language_after = self.unfold_rg(w_after, "emd", False)
            improvement = abs(emd_after - emd_before)/emd_after
            if v >= 1:
                loglh_after = self.log_likelihood(ev_language, pn_language_after)
                print("\temd value at step " + str((i+1)*step) + " : " + str(emd_after) + " (" + str(improvement) + ")")
                print("\tloglh value at step " + str((i+1)*step) + " : " + str(loglh_after))

            if improvement <= improvement_threshold:
                break
            w_before = w_after
            emd_before = emd_after

        print("\t#############################################")
        loglh_after = self.log_likelihood(ev_language, pn_language_after)
        print("\temd value with optimized values: " + str(emd_after))
        print("\tloglh value with optimized values: " + str(loglh_after))
        print("\tExecution time: " + str(time.time() - start_time) + " seconds")
        return {transition.name: w_after[i] for i, transition in enumerate(self.pn.get_transitions())}
