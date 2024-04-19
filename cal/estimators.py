import pm4py


def frequency_estimator(ev, pn):
    weights = {}
    for transition in pn.get_transitions():
        weights[transition.get_name()] = 0
        for trace in ev.get_language()[:-1]:
            weights[transition.get_name()] += trace.get_seq().count(transition.get_label()) * trace.get_freq()
        if weights[transition.get_name()] == 0:
            weights[transition.get_name()] = 1
    return weights


def lhpair_estimators(ev, pn):
    weights = {}
    for transition in pn.get_transitions():
        weights[transition.get_name()] = 0
        left_pair_transitions = [outarc.source for outarc in pn.get_outarcs() if outarc.target in
                                 [inarc.source for inarc in pn.get_inarcs() if inarc.target == transition]]
        for trace in ev.get_language()[:-1]:
            if trace.get_seq()[0] == transition.get_label():
                weights[transition.get_name()] += trace.get_freq()
            if trace.get_seq()[-1] == transition.get_label():
                weights[transition.get_name()] += trace.get_freq()
            for left_pair_transition in left_pair_transitions:
                for i in range(len(trace.get_seq()) - 1):
                    if trace.get_seq()[i] == left_pair_transition.get_label() and trace.get_seq()[i+1] == transition.get_label():
                        weights[transition.get_name()] += trace.get_freq()
        if weights[transition.get_name()] == 0:
            weights[transition.get_name()] = 1
    return weights


def rhpair_estimators(ev, pn):
    weights = {}
    for transition in pn.get_transitions():
        weights[transition.get_name()] = 0
        right_pair_transitions = [inarcs.target for inarcs in pn.get_inarcs() if inarcs.source in
                                  [outarc.target for outarc in pn.get_outarcs() if outarc.source == transition]]
        for trace in ev.get_language()[:-1]:
            if trace.get_seq()[0] == transition.get_label():
                weights[transition.get_name()] += trace.get_freq()
            if trace.get_seq()[-1] == transition.get_label():
                weights[transition.get_name()] += trace.get_freq()
            for right_pair_transition in right_pair_transitions:
                for i in range(len(trace.get_seq()) - 1):
                    if trace.get_seq()[i] == transition.get_label() and trace.get_seq()[i+1] == right_pair_transition.get_label():
                        weights[transition.get_name()] += trace.get_freq()
        if weights[transition.get_name()] == 0:
            weights[transition.get_name()] = 1
    return weights


def pairscale_estimators(ev, pn):
    aveg_trans_freq = sum([len(trace.get_seq()) * trace.get_freq() for trace in ev.get_language()[:-1]]) / len(pn.get_transitions())
    weights = {}
    for transition in pn.get_transitions():
        weights[transition.get_name()] = 0
        right_pair_transitions = [inarcs.target for inarcs in pn.get_inarcs() if inarcs.source in
                                  [outarc.target for outarc in pn.get_outarcs() if outarc.source == transition]]
        for trace in ev.get_language()[:-1]:
            if trace.get_seq()[0] == transition.get_label():
                weights[transition.get_name()] += trace.get_freq()
            if trace.get_seq()[-1] == transition.get_label():
                weights[transition.get_name()] += trace.get_freq()
            for right_pair_transition in right_pair_transitions:
                for i in range(len(trace.get_seq())-1):
                    if trace.get_seq()[i] == transition.get_label() and trace.get_seq()[i+1] == right_pair_transition.get_label():
                        weights[transition.get_name()] += trace.get_freq()
        weights[transition.get_name()] = weights[transition.get_name()] / aveg_trans_freq
        if weights[transition.get_name()] == 0:
            weights[transition.get_name()] = 1
    return weights


def fork_estimators(ev, pn):

    places_weights = {}
    for place in pn.get_places():
        if place.get_name() == "source":
            places_weights[place.get_name()] = len(ev.get_language())-1
        else:
            places_weights[place.get_name()] = 0
            incomming_transitions = [outarc.source for outarc in pn.get_outarcs() if outarc.target == place]
            outcomming_transitions = [inarc.target for inarc in pn.get_inarcs() if inarc.source == place]
            for trace in ev.get_language()[:-1]:
                for incomming_transition in incomming_transitions:
                    for outcomming_transition in outcomming_transitions:
                        for i in range(len(trace.get_seq())-1):
                            if trace.get_seq()[i] == incomming_transition.get_name() and trace.get_seq()[i+1] == outcomming_transition.get_name():
                                places_weights[place.get_name()] += trace.get_freq()
            if places_weights[place.get_name()] == 0:
                places_weights[place.get_name()] = 1

    weights = {}
    frequency_weights = frequency_estimator(ev, pn)
    for transition in pn.get_transitions():
        weights[transition.get_name()] = 0
        incomming_places = [inarc.source for inarc in pn.get_inarcs() if inarc.target == transition]
        for incomming_place in incomming_places:
            weights[transition.get_name()] += (places_weights[incomming_place.get_name()] *
                                               (frequency_weights[transition.get_name()] /
                                                sum([frequency_weights[inarc.target.get_name()] for inarc in pn.get_inarcs() if inarc.source == incomming_place])))

    return weights

def align_estimator(log, net):
    weights = {transition.get_name(): 0 for transition in net.get_transitions()}
    for i, align in enumerate(pm4py.algo.conformance.alignments.petri_net.algorithm.apply_log(log.get_log(), net.get_net(), net.get_im(), net.get_fm())):
        for move in align["alignment"]:
            for transition in net.get_transitions():
                if transition.get_name() == move[1]:
                    weights[transition.get_name()] += log.get_traces()[i].get_freq()
                    break
    return weights
