import copy
import pm4py
import sympy
import xml.etree.ElementTree as et
from pm4py.visualization.transition_system import visualizer as rg_visualizer

import language.Eventlog


class Place:

    def __init__(self, name, tokens=0):
        self.name = name
        self.tokens = tokens

    def get_name(self):
        return self.name

    def get_tokens(self):
        return self.tokens

    def __repr__(self):
        return "(" + str(self.name) + "," + str(self.tokens) + ")"


class Transition:

    def __init__(self, name, label, weight=1):
        self.name = name
        self.label = label
        self.weight = weight
        self.silent = self.label is None

    def get_name(self):
        return self.name

    def get_label(self):
        return self.label

    def get_weight(self):
        return self.weight

    def set_weight(self, weight):
        self.weight = weight

    def is_silent(self):
        return self.silent

    def __repr__(self):
        return "(" + str(self.name) + "," + str(self.label) + "," + str(self.weight) + ")"


class Arc:

    def __init__(self, source, target, multiplicity):
        self.source = source
        self.target = target
        self.multiplicity = multiplicity

    def get_source(self):
        return self.source

    def get_target(self):
        return self.target

    def get_multiplicity(self):
        return self.multiplicity

    def __repr__(self):
        return ("(" + str(self.source) + "," +
                str(self.target) + "," +
                str(self.multiplicity) + ")")


# P -> T
class InArc(Arc):

    def __init__(self, source, target, multiplicity=1):
        super().__init__(source, target, multiplicity)


# T -> P
class OutArc(Arc):

    def __init__(self, source, target, multiplicity=1):
        super().__init__(source, target, multiplicity)


class Petrinet:

    def __init__(self, net, im, fm, places, transitions, inarcs, outarcs, dir_places, dir_transitions, tr2lab, start, end,
                 max_activities):

        self.net, self.im, self.fm = net, im, fm
        self.places, self.transitions = places, transitions
        self.inarcs, self.outarcs = inarcs, outarcs
        self.dir_places, self.dir_transitions, self.tr2lab = dir_places, dir_transitions, tr2lab
        self.start, self.end = start, end
        self.max_activities = max_activities

        self.rg = pm4py.objects.petri_net.utils.reachability_graph.construct_reachability_graph(self.net, self.im, True)

    def get_net(self):
        return self.net

    def get_im(self):
        return self.im

    def get_fm(self):
        return self.fm

    def get_places(self):
        return self.places

    def get_transitions(self):
        return self.transitions

    def get_inarcs(self):
        return self.inarcs

    def get_outarcs(self):
        return self.outarcs

    def get_dir_places(self):
        return self.dir_places

    def get_dir_transitions(self):
        return self.dir_transitions

    def get_tr2lab(self):
        return self.tr2lab

    def get_rg(self):
        return self.rg

    def set_weights(self, weights):
        for transition in self.transitions:
            transition.set_weight(weights[transition.get_name()])

    def __str__(self):
        if self.net is not None:
            pm4py.view_petri_net(self.net, self.im, self.fm)
        result = "Places: ["
        for pl in self.places:
            result += pl.__str__() + ","
        result = result[:-1] + "]\n" + "Transitions: ["
        for tr in self.transitions:
            result += tr.__str__() + ","
        result = result[:-1] + "]\n" + "Inarcs: ["
        for ia in self.inarcs:
            result += ia.__str__() + ","
        result = result[:-1] + "]\n" + "Outarcs: ["
        for oa in self.outarcs:
            result += oa.__str__() + ","
        return result[:-1] + "]"

    def show_rg(self):
        gviz = rg_visualizer.apply(self.rg, parameters={rg_visualizer.Variants.VIEW_BASED.value.Parameters.FORMAT: "png"})
        rg_visualizer.view(gviz)

    ####################################################################################################

    def get_enabled_transitions(self, marking):
        enabled_transitions = []
        for transition in self.get_transitions():
            flag = True
            for place in [ia.source for ia in self.inarcs if ia.target == transition]:
                if place not in marking:
                    flag = False
            if flag:
                enabled_transitions.append(transition)
        return enabled_transitions

    def update(self, cpt_activities, mass, marking, enabled_transitions, enabled_transition):
        new_cpt_activities = copy.deepcopy(cpt_activities)
        if not enabled_transition.is_silent():
            new_cpt_activities[enabled_transition.get_label()] += 1
        new_mass = (mass if mass != 0 else 1) * (enabled_transition.get_weight() / sum(
            [enabled_transition.get_weight() for enabled_transition in enabled_transitions]))
        new_marking = [oa.target for oa in self.outarcs if oa.source == enabled_transition]
        for place in marking:
            if place not in [ia.source for ia in self.inarcs if ia.target == enabled_transition]:
                new_marking.append(place)
        return new_cpt_activities, new_mass, new_marking

    def build_mfc(self, max_mass):

        probabilities, markings, firing_sequences, conflicts = [], [], [], []

        def build_mfc_rec(cpt_activities, mass, marking, firing_sequence, conflict):
            nonlocal probabilities, markings, firing_sequences, conflicts

            if marking not in markings:
                markings.append(marking)
            enabled_transitions = self.get_enabled_transitions(marking)

            if len(enabled_transitions) == 0:
                if firing_sequence in firing_sequences:
                    idx = firing_sequences.index(firing_sequence)
                    probabilities[idx] += mass
                    conflicts[idx] += ["+"] + conflict
                else:
                    probabilities.append(mass)
                    firing_sequences.append(firing_sequence)
                    conflicts.append(conflict)

            elif sum(probabilities) < max_mass:
                enabled_transitions.sort(key=lambda transition: transition.get_weight(), reverse=True)
                for enabled_transition in enabled_transitions:
                    new_cpt_activities, new_mass, new_marking = self.update(cpt_activities, mass, marking,
                                                                            enabled_transitions, enabled_transition)
                    if enabled_transition.is_silent():
                        build_mfc_rec(cpt_activities, new_mass, new_marking,
                                      firing_sequence,
                                      conflict + [(enabled_transition, enabled_transitions)])
                    else:
                        if new_cpt_activities[enabled_transition.get_label()] <= self.max_activities[
                            enabled_transition.get_label()]:
                            build_mfc_rec(new_cpt_activities, new_mass, new_marking,
                                          firing_sequence + [enabled_transition],
                                          conflict + [(enabled_transition, enabled_transitions)])

        build_mfc_rec({activity: 0 for activity in self.max_activities}, 1, [self.start], [], [])
        return probabilities, markings, firing_sequences, conflicts

    ####################################################################################################

    def get_language(self, max_mass):
        probabilities, markings, firing_sequences, conflicts = self.build_mfc(max_mass)
        pn_language = []
        for i, firing_sequence in enumerate(firing_sequences):
            pn_language.append(
                language.Eventlog.Trace(i, "".join([transition.get_label() for transition in firing_sequence]), None,
                                        probabilities[i], None))
        pn_language.append(language.Eventlog.Trace(len(pn_language), "", 0, 0, -1))
        return pn_language

    def build_eqs(self, max_mass):
        probabilities, markings, firing_sequences, conflicts = self.build_mfc(max_mass)
        eqs = {}
        for i, conflict in enumerate(conflicts):
            trace, eq = "", ""
            for c in conflict:
                if c == "+":
                    eq = "+"
                else:
                    trace += (str(c[0].get_label()) if c[0].get_label() != None else "")
                    eq += "(" + str(c[0].get_name()) + "/("
                    for term in c[1]:
                        eq += str(term.get_name()) + "+"
                    eq = eq[:-1] + "))*"
            eq = eq[:-1]
            eqs[trace] = sympy.sympify(eq)
        return eqs

    def compute_eqs(self, max_mass, weights):
        pn_language = []
        eqs = self.build_eqs(max_mass)
        for i, trace in enumerate(eqs):
            pn_language.append(language.Eventlog.Trace(i, trace, None, (eqs[trace]).subs(weights), None))
        return pn_language

    ####################################################################################################

    def export_pnml(self, filename, weights):

        pm4py.write_pnml(self.net, self.im, self.fm, "gspn/" + filename)

        tree, t = et.parse("gspn/" + filename + ".pnml"), 0
        root = tree.getroot()
        for node in root[0][1]:
            if node.tag == "transition":
                if len(node) > 1:
                    node.remove(node.find("toolspecific"))
                toolspecific = et.SubElement(node, "toolspecific")
                toolspecific.set("tool", "StochasticPetriNet")
                toolspecific.set("version", "0.2")
                pr_distributionType = et.SubElement(toolspecific, "property")
                pr_distributionType.set("key", "distributionType")
                pr_distributionType.text = "IMMEDIATE"
                pr_weight = et.SubElement(toolspecific, "property")
                pr_weight.set("key", "weight")
                pr_weight.text = str(weights["t"+str(t)])
                pr_invisible = et.SubElement(toolspecific, "property")
                pr_invisible.set("key", "invisible")
                pr_invisible.text = "true" if node[0][0].text == ("t"+str(t)) else "false"
                pr_priority = et.SubElement(toolspecific, "property")
                pr_priority.set("key", "priority")
                pr_priority.text = "1"
                t+=1
        tree.write("gspn/" + filename + ".pnml")

    ####################################################################################################

    def gen_gspn(self):

        gspn = "\n"

        # Const
        for transition in self.transitions:
            gspn += "const double w_" + transition.name + " = 1;\n"
        gspn += "\n"

        # NbPlaces
        gspn += "NbPlaces = " + str(len(self.places)) + ";\n"
        # NbTransitions
        gspn += "NbTransitions = " + str(len(self.transitions)) + ";\n\n"

        # PlacesList
        gspn += "PlacesList = {"
        for place in self.places:
            gspn += place.name + ","
        gspn = gspn[:-1] + "};\n"
        # TransitionsList
        gspn += "TransitionsList = {"
        for transition in self.transitions:
            gspn += transition.name + ","
        gspn = gspn[:-1] + "};\n\n"

        # Marking
        gspn += "Marking = {"
        for place in self.places:
            if place.name == "start":
                gspn += "(start,1);"
            else:
                gspn += "(" + place.name + ",0);"
        gspn += "};\n\n"

        # Transitions
        gspn += "Transitions = {"
        for transition in self.transitions:
            gspn += "(" + transition.name + ",EXPONENTIAL(w_" + transition.name + "),1,1,SINGLE);"
        gspn += "};\n\n"

        # InArcs
        gspn += "InArcs = {"
        for arc in self.inarcs:
            gspn += "(" + str(arc.source.name) + "," + str(arc.target.name) + "," + str(arc.multiplicity) + ");"
        gspn += "};\n"

        # OutArcs
        gspn += "OutArcs = {"
        for arc in self.outarcs:
            gspn += "(" + str(arc.source.name) + "," + str(arc.target.name) + "," + str(arc.multiplicity) + ");"
        gspn += "};\n"

        return gspn
