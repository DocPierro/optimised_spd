import numpy as np
import pm4py
import sys
import time
import random
import functools
import multiprocessing
from matplotlib import pyplot as plt
from queue import Queue
from scipy.optimize import minimize
from Levenshtein import distance as levenshtein_distance
from pm4py.algo.evaluation.earth_mover_distance.variants import pyemd

class Optimizer:

    def __init__(self, ev, pn):
        self.ev, self.pn = ev, pn
        self.rg = pn.get_rg()
        self.ev_language = self.ev.get_language_pm4py()
        self.pn_language_memoization = None

    def emd_norm_lev(self, l1, l2, unit=False):

        nc1, nc2 = 0, 0
        words = set(l1) & set(l2)
        for w in words:
            if w in l1:
                nc1 = nc1 + l1[w]
            if w in l2:
                nc2 = nc2 + l2[w]

        hist1 = np.zeros(len(words))
        hist2 = np.zeros(len(words))
        distance = np.zeros((len(words), len(words)))
        for x, xw in enumerate(words):
            for y, yw in enumerate(words):
                if xw in l1:
                    hist1[x] = l1[xw] / nc1
                if yw in l2:
                    hist2[y] = l2[yw] / nc2
                if len(xw) == 0 or len(yw) == 0:
                    distance[x, y] = 1.
                else:
                    if unit:
                        distance[x, y] = 0 if xw == yw else 1  # uEMSC
                    else:
                        maxl = float(max(len(xw), len(yw)))
                        distance[x, y] = float(levenshtein_distance(xw, yw)) / maxl


        return pyemd.emd(hist1, hist2, distance)

    def log_likelihood(self, ev_language, pn_language):
        lh = 0
        for trace in pn_language:
            lh += ev_language[trace] * np.log(pn_language[trace])
        return -lh

    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################

    def compute_rg(self, w):
        weights = {transition.get_name(): w[t] for t, transition in enumerate(self.pn.get_transitions())}
        for transition in self.rg.transitions:
            denom = 0
            for enabled_transition in transition.from_state.outgoing:
                denom += weights[enabled_transition.name]
            # use property data to memorize prob
            transition.data = weights[transition.name] / denom

    def unfold_rg(self, w):

        self.compute_rg(w)
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
                if tt in self.ev_language:
                    if tt not in traces:
                        traces[tt] = prob
                    else:
                        traces[tt] = traces[tt] + prob

        return traces

    def compute_unfolding(self, w, d):

        pn_language = self.unfold_rg(w)

        if d == "KLD":
            return self.log_likelihood(self.ev_language, pn_language)
        elif d == "EMD":
            return self.emd_norm_lev(self.ev_language, pn_language, unit=False)
        elif d == "both":
            return self.log_likelihood(self.ev_language, pn_language), self.emd_norm_lev(self.ev_language, pn_language, unit=False)

    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################

    def unode_function_memoized(self, w, traversed, incoming_data, key):

        if len(incoming_data) == 0:
            return 1
        if key in traversed:
            return traversed[key]

        probability = 0
        for firing_transition, enabling_transitions, prob_function, _ in incoming_data:
            probability += ((w[firing_transition] / sum([w[enabling_transition] for enabling_transition in enabling_transitions]))
                            * prob_function(w, traversed))
        traversed[key] = probability
        return probability

    def create_unode_function(self, incoming_data, key):
        return functools.partial(self.unode_function_memoized, incoming_data=incoming_data, key=key)

    def unode_derivative_memoized(self, w, traversed_prob, traversed_der, l, incoming_data, key):

        if len(incoming_data) == 0:
            return 0
        if not(traversed_der[key][l] is None):
            return traversed_der[key][l]

        der = 0
        for firing_transition, enabling_transitions, prob_function, der_function in incoming_data:
            sum_enabled = sum([w[j] for j in enabling_transitions])
            der += der_function[l](w, traversed_prob, traversed_der) * (w[firing_transition] / sum_enabled)
            if not(l in enabling_transitions):
                der += 0
            elif l != firing_transition:
                der += prob_function(w, traversed_prob) * -(w[firing_transition] / sum_enabled**2)
            else:
                der += prob_function(w, traversed_prob) * (sum([w[j] for j in enabling_transitions if j!=l]) / sum_enabled**2)
        traversed_der[key][l] = der
        return der

    def create_unode_derivatives(self, incoming_data, key):
        unode_derivatives = [functools.partial(self.unode_derivative_memoized, l=l, incoming_data=incoming_data, key=key) for l in range(len(self.pn.get_transitions()))]
        return unode_derivatives

    def unfold_rg_memoization(self):

        max_length = max([len(trace.get_seq()) for trace in self.ev.get_language()])
        prefixes = self.ev.build_prefix_set()

        transitions_names = list(map(lambda transition: transition.name, self.pn.get_transitions()))

        # find initial state in the reachability graph
        init_state = None
        for state in self.rg.states:
            if state.name == 'start1':
                init_state = state

        q = Queue(0)
        q.put((init_state, [], 0))

        unodes_data, unodes_functions, unodes_derivatives = {}, {}, {}
        unodes_data[(init_state, tuple([]), 0)] = []

        traces = {}  # all traces that arrives to sink
        while not q.empty():

            state, trace, level = q.get()
            key = (state, tuple(trace), level)
            unodes_functions[key] = self.create_unode_function(unodes_data[key], key)
            unodes_derivatives[key] = self.create_unode_derivatives(unodes_data[key], key)
            if len(state.outgoing) > 0:
                for transition in state.outgoing:
                    new_state = transition.to_state
                    new_trace = trace.copy()
                    label = self.pn.get_tr2lab()[transition.name]
                    if label is not None:
                        new_trace.append(label)
                    if len(new_trace) <= max_length:
                        if tuple(new_trace) in prefixes:
                            new_key = (new_state, tuple(new_trace), level + 1)
                            if new_key not in unodes_data:
                                unodes_data[new_key] = [(transitions_names.index(transition.name),
                                                         [transitions_names.index(tr.name) for tr in state.outgoing],
                                                         unodes_functions[key], unodes_derivatives[key])]
                                q.put((new_state, new_trace, level + 1))
                            else:
                                unodes_data[new_key].append((transitions_names.index(transition.name),
                                                             [transitions_names.index(tr.name) for tr in state.outgoing],
                                                             unodes_functions[key], unodes_derivatives[key]))
            else:
                tt = tuple(trace)
                if tt in self.ev_language:
                    if tt not in traces:
                        traces[tt] = [[unodes_functions[key], unodes_derivatives[key]]]
                    else:
                        traces[tt].append([unodes_functions[key], unodes_derivatives[key]])

        return traces, unodes_functions, unodes_derivatives

    ###################################################################################################################

    def compute_unfolding_memoized(self, w, d):

        pn_language, traversed_prob = {}, {}

        for key in self.pn_language_memoization:
            prob = 0
            for f in self.pn_language_memoization[key]:
                prob += f[0](w, traversed_prob)
            pn_language[key] = prob

        if d == "KLD":
            return self.log_likelihood(self.ev_language, pn_language)
        elif d == "EMD":
            return self.emd_norm_lev(self.ev_language, pn_language, unit=False)
        elif d == "both":
            return self.log_likelihood(self.ev_language, pn_language), self.emd_norm_lev(self.ev_language, pn_language, unit=False)

    def compute_unfolding_memoized_der(self, w):

        pn_language, traversed_prob = {}, {}

        for key in self.pn_language_memoization:
            prob = 0
            for f in self.pn_language_memoization[key]:
                prob += f[0](w, traversed_prob)
            pn_language[key] = prob

        derivatives = np.zeros(len(self.pn.get_transitions()))
        traversed_der = {key: [None] * len(self.pn.get_transitions()) for key in traversed_prob}
        for l in range(len(derivatives)):
            for key in self.pn_language_memoization:
                derivatives[l] += self.ev_language[key] * (1 / pn_language[key]) * sum(
                    [self.pn_language_memoization[key][i][1][l](w, traversed_prob, traversed_der) for i in
                     range(len(self.pn_language_memoization[key]))])

        return self.log_likelihood(self.ev_language, pn_language), -derivatives

    def jac(self, w):

        pn_language, traversed_prob = {}, {}

        for key in self.pn_language_memoization:
            prob = 0
            for f in self.pn_language_memoization[key]:
                prob += f[0](w, traversed_prob)
            pn_language[key] = prob

        derivatives = np.zeros(len(self.pn.get_transitions()))
        traversed_der = {key: [None] * len(self.pn.get_transitions()) for key in traversed_prob}
        for l in range(len(derivatives)):
            for key in self.pn_language_memoization:
                derivatives[l] += self.ev_language[key] * (1 / pn_language[key]) * sum([self.pn_language_memoization[key][i][1][l](w, traversed_prob, traversed_der) for i in range(len(self.pn_language_memoization[key]))])

        return -derivatives

    ###################################################################################################################

    def compute_starting_vector(self, i, memoized):
        np.random.seed(i)
        w = (1 - 0.001) * np.random.rand(len(self.pn.get_transitions())) + 0.001
        KLD, EMD = self.compute_unfolding_memoized(w, "both") if memoized else self.compute_unfolding(w, "both")
        return (w, KLD, EMD)

    def optimize(self, param):

        i, m, solver, maxiter, memoized, derivatives, bnds, w_before, KLD_before, EMD_before = param

        start_time = time.time()
        if memoized:
            if derivatives:
                min = minimize(self.compute_unfolding_memoized_der, w_before, bounds=bnds, method=solver, jac=True, options={'disp': False, 'maxiter': maxiter})
                mbis_after = self.compute_unfolding(min.x, "EMD" if m=="KLD" else "KLD")
            else:
                min = minimize(self.compute_unfolding_memoized, w_before, m, bounds=bnds, method=solver, options={'disp': False, 'maxiter': maxiter})
                mbis_after = self.compute_unfolding(min.x, "EMD" if m=="KLD" else "KLD")
        else:
            min = minimize(self.compute_unfolding, w_before, m, bounds=bnds, method=solver, options={'disp': False, 'maxiter': maxiter})
            mbis_after = self.compute_unfolding(min.x, "EMD" if m=="KLD" else "KLD")
        time_elasped = time.time() - start_time

        print(f"KLD / EMD value with starting values: {KLD_before:<6} / {EMD_before:<6}")
        print(f"Optimize " + str(m) + " with " + ("memoized" if memoized else "non-memoized") + " unfolding (" + ("wth" if derivatives else "w/o") + " derivatives):")
        print(f"\t\t " + ("KLD / EMD" if m=="KLD" else "EMD / KLD") + f" value with optimized values ({m:<3} | {maxiter:<3}): {min.fun:<6} / {mbis_after:<6} in {time_elasped:<10.2f} seconds ({min.nit:<4} iterations)")

        return True

    def estimate(self, m, nw0, solver, maxiter, memoized=False, derivatives=False):
        np.set_printoptions(threshold=np.inf, linewidth=np.inf)
        sys.setrecursionlimit(8000)

        start_time = time.time()

        print("OPT-based SPD | Minimize " + str(m) + " | nw0: " + str(nw0) + " | solver: " + str(solver) + " | maxiter: " + str(maxiter) + " | memoized: " + str(memoized) + " | derivatives: " + str(derivatives) + "\n")

        if derivatives and m=="EMD":
            raise Exception("Impossible to use derivative mode with EMD.")
        if derivatives and not(memoized):
            raise Exception("Impossible to use derivative mode with non-memoized unfolding.")

        if memoized:
            print("####################################################################################################")
            print("Memoized unfolding setup started at " + str(time.time() - start_time) + " seconds")
            self.pn_language_memoization, self.unodes_functions, self.unodes_derivatives = self.unfold_rg_memoization()
            print("Memoized unfolding setup ended at " + str(time.time() - start_time) + " seconds")
            print("####################################################################################################\n")

        print("####################################################################################################")
        print("Looking for starting points (" + str(nw0) + " vectors) started at " + str(time.time() - start_time) + " seconds")
        bnds = ((0.001, 1),) * len(self.pn.get_transitions())
        ws = [(index, m, solver, maxiter, memoized, derivatives, bnds, *tup) for index, tup in enumerate(sorted([self.compute_starting_vector(i, memoized) for i in range(nw0)], key=lambda x:x[1]))]
        print("Looking for starting points (" + str(nw0) + " vectors) ended at " + str(time.time() - start_time) + " seconds")
        print("####################################################################################################\n")

        print("####################################################################################################")
        print("Optimising vectors started at " + str(time.time() - start_time) + " seconds")
        result = self.optimize(ws[0])
        print("Optimising vectors ended at " + str(time.time() - start_time) + " seconds")
        print("####################################################################################################\n")

        print("Execution time: " + str(time.time() - start_time) + " seconds\n")

        return result
