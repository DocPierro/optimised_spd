import pm4py
import xml.etree.ElementTree as ET
from pm4py.statistics.variants.log import get as variants_module

from model.Petrinet import Petrinet, Place, Transition, InArc, OutArc


class Event:

    def __init__(self, case, activity, timestamp):
        self.case = case
        self.activity = activity
        self.timestamp = timestamp

    def get_case(self):
        return self.case

    def get_activity(self):
        return self.activity

    def get_timestamp(self):
        return self.timestamp

    def __repr__(self):
        return ("Case: " + str(self.case) +
                "; Activity: " + str(self.activity) +
                "; Timestamp: " + str(self.timestamp) + ";")


class Trace:

    def __init__(self, id, seq, word, freq, prob, mapping):
        self.id = id
        self.seq = seq
        self.word = word
        self.freq = freq
        self.prob = prob
        self.mapping = mapping

    def get_id(self):
        return self.id

    def get_word(self):
        return self.word

    def get_seq(self):
        return self.seq

    def get_freq(self):
        return self.freq

    def get_prob(self):
        return self.prob

    def map(self):
        return self.mapping

    def __str__(self):
        return "(" + str(self.id) + ",<" + str(self.word) + ">," + str(self.freq) + "," + str(self.prob) + "," + str(self.mapping) + ")"

    def __repr__(self):
        return self.__str__()


class Eventlog:

    def __init__(self, filename):

        self.log = pm4py.read.read_xes(filename)
        self.activities = {activity: mapping+1 for mapping, activity in
                           enumerate(pm4py.get_event_attribute_values(self.log, "concept:name"))}

        self.language, i = [], 0
        self.language_pm4py = variants_module.get_language(self.log)
        if "case:frequency" in self.log:
            #frequencies = pm4py.get_event_attribute_values(self.log, "case:frequency")
            frequencies = [int(trace.find("int").attrib["value"]) for trace in ET.parse(filename).getroot().findall("trace")]
            for i, seq in enumerate(self.language_pm4py):
                trace_mapping = 0
                for c, activity in enumerate(seq):
                    trace_mapping += self.activities[activity] * pow(len(self.activities), c)
                self.language.append(Trace(i, seq, "".join(seq), int(frequencies[i]), int(frequencies[i])/sum(frequencies), trace_mapping))
        else:
            cases = pm4py.get_event_attribute_values(self.log, "case:concept:name")
            for i, seq in enumerate(self.language_pm4py):
                trace_mapping = 0
                for c, activity in enumerate(seq):
                    trace_mapping += self.activities[activity] * pow(len(self.activities), c)
                self.language.append(Trace(i, seq, "".join(seq), self.language_pm4py[seq]*len(cases), self.language_pm4py[seq], trace_mapping))
        self.language.append(Trace(len(self.language), (), "", 0, 0, -1))

    def get_log(self):
        return self.log

    def get_activities(self):
        return self.activities

    def get_language_pm4py(self):
        return self.language_pm4py

    def get_language(self):
        return self.language

    def __str__(self):
        result = "["
        for trace in self.language:
            result += trace.__str__() + ",\n"
        return result[:-2] + "]"

    ####################################################################################################

    def discover_pn_alpha(self):

        net, im, fm = pm4py.discovery.discover_petri_net_alpha(self.log,
                                                               case_id_key="case:concept:name",
                                                               activity_key="concept:name",
                                                               timestamp_key="time:timestamp")

        places, dir_places = [], {}
        pl_start, pl_end = None, None
        for p, place in enumerate(net.places):
            if place.name == "start":
                place.name = "start"
                pl_start = Place("start")
                dir_places[place] = pl_start
                places.append(pl_start)
            elif place.name == "end":
                place.name = "sink"
                pl_end = Place("end")
                dir_places[place] = pl_end
                places.append(pl_end)
            else:
                place.name = "p" + str(p)
                pl_i = Place("p" + str(p))
                dir_places[place] = pl_i
                places.append(pl_i)

        transitions, dir_transitions, tr2lab = [], {}, {}
        for t, transition in enumerate(net.transitions):
            transition.name = "t" + str(t)
            tr_i = Transition("t" + str(t), transition.label)
            dir_transitions[transition] = tr_i
            transitions.append(tr_i)
            tr2lab["t" + str(t)] = transition.label

        inarcs, outarcs = [], []
        for arc in net.arcs:
            if isinstance(arc.source, pm4py.objects.petri_net.obj.PetriNet.Place):
                source = dir_places[arc.source]
                target = dir_transitions[arc.target]
                inarcs.append(InArc(source, target))
            elif isinstance(arc.source, pm4py.objects.petri_net.obj.PetriNet.Transition):
                source = dir_transitions[arc.source]
                target = dir_places[arc.target]
                outarcs.append(OutArc(source, target))

        return Petrinet(net, im, fm, places, transitions, inarcs, outarcs, dir_places, dir_transitions, tr2lab,
                        pl_start, pl_end, self.get_max_activities())

    def discover_pn_heuristics(self):

        net, im, fm = pm4py.discovery.discover_petri_net_heuristics(self.log,
                                                                    case_id_key="case:concept:name",
                                                                    activity_key="concept:name",
                                                                    timestamp_key="time:timestamp")

        places, dir_places = [], {}
        pl_start, pl_end = None, None
        for p, place in enumerate(net.places):
            if place.name == "source0":
                place.name = "start"
                pl_start = Place("start")
                dir_places[place] = pl_start
                places.append(pl_start)
            elif place.name == "sink0":
                place.name = "sink"
                pl_end = Place("end")
                dir_places[place] = pl_end
                places.append(pl_end)
            else:
                place.name = "p" + str(p)
                pl_i = Place("p" + str(p))
                dir_places[place] = pl_i
                places.append(pl_i)

        transitions, dir_transitions, tr2lab = [], {}, {}
        for t, transition in enumerate(net.transitions):
            transition.name = "t" + str(t)
            tr_i = Transition("t" + str(t), transition.label)
            dir_transitions[transition] = tr_i
            transitions.append(tr_i)
            tr2lab["t" + str(t)] = transition.label

        inarcs, outarcs = [], []
        for arc in net.arcs:
            if isinstance(arc.source, pm4py.objects.petri_net.obj.PetriNet.Place):
                source = dir_places[arc.source]
                target = dir_transitions[arc.target]
                inarcs.append(InArc(source, target))
            elif isinstance(arc.source, pm4py.objects.petri_net.obj.PetriNet.Transition):
                source = dir_transitions[arc.source]
                target = dir_places[arc.target]
                outarcs.append(OutArc(source, target))

        return Petrinet(net, im, fm, places, transitions, inarcs, outarcs, dir_places, dir_transitions, tr2lab,
                        pl_start, pl_end, self.get_max_activities())

    def discover_pn_inductive(self):

        net, im, fm = pm4py.discovery.discover_petri_net_inductive(self.log,
                                                                   case_id_key="case:concept:name",
                                                                   activity_key="concept:name",
                                                                   timestamp_key="time:timestamp")

        places, dir_places = [], {}
        pl_start, pl_end = None, None
        for p, place in enumerate(net.places):
            if place.name == "source":
                place.name = "start"
                pl_start = Place("start")
                dir_places[place] = pl_start
                places.append(pl_start)
            elif place.name == "sink":
                place.name = "sink"
                pl_end = Place("end")
                dir_places[place] = pl_end
                places.append(pl_end)
            else:
                place.name = "p" + str(p)
                pl_i = Place("p" + str(p))
                dir_places[place] = pl_i
                places.append(pl_i)

        transitions, dir_transitions, tr2lab = [], {}, {}
        for t, transition in enumerate(net.transitions):
            transition.name = "t" + str(t)
            tr_i = Transition("t" + str(t), transition.label)
            dir_transitions[transition] = tr_i
            transitions.append(tr_i)
            tr2lab["t" + str(t)] = transition.label

        inarcs, outarcs = [], []
        for arc in net.arcs:
            if isinstance(arc.source, pm4py.objects.petri_net.obj.PetriNet.Place):
                source = dir_places[arc.source]
                target = dir_transitions[arc.target]
                inarcs.append(InArc(source, target))
            elif isinstance(arc.source, pm4py.objects.petri_net.obj.PetriNet.Transition):
                source = dir_transitions[arc.source]
                target = dir_places[arc.target]
                outarcs.append(OutArc(source, target))

        return Petrinet(net, im, fm, places, transitions, inarcs, outarcs, dir_places, dir_transitions, tr2lab,
                        pl_start, pl_end, self.get_max_activities())

    def discover_pn_ilp(self):

        net, im, fm = pm4py.discovery.discover_petri_net_ilp(self.log,
                                                             case_id_key="case:concept:name",
                                                             activity_key="concept:name",
                                                             timestamp_key="time:timestamp")

        places, dir_places = [], {}
        pl_start, pl_end = None, None
        for p, place in enumerate(net.places):
            if place.name == "source":
                place.name = "start"
                pl_start = Place("start")
                dir_places[place] = pl_start
                places.append(pl_start)
            elif place.name == "sink":
                place.name = "sink"
                pl_end = Place("end")
                dir_places[place] = pl_end
                places.append(pl_end)
            else:
                place.name = "p" + str(p)
                pl_i = Place("p" + str(p))
                dir_places[place] = pl_i
                places.append(pl_i)

        transitions, dir_transitions, tr2lab = [], {}, {}
        for t, transition in enumerate(net.transitions):
            transition.name = "t" + str(t)
            tr_i = Transition("t" + str(t), transition.label)
            dir_transitions[transition] = tr_i
            transitions.append(tr_i)
            tr2lab["t" + str(t)] = transition.label

        inarcs, outarcs = [], []
        for arc in net.arcs:
            if isinstance(arc.source, pm4py.objects.petri_net.obj.PetriNet.Place):
                source = dir_places[arc.source]
                target = dir_transitions[arc.target]
                inarcs.append(InArc(source, target))
            elif isinstance(arc.source, pm4py.objects.petri_net.obj.PetriNet.Transition):
                source = dir_transitions[arc.source]
                target = dir_places[arc.target]
                outarcs.append(OutArc(source, target))

        return Petrinet(net, im, fm, places, transitions, inarcs, outarcs, dir_places, dir_transitions, tr2lab,
                        pl_start, pl_end, self.get_max_activities())

    ####################################################################################################

    def build_prefix_set(self):
        prefixes = set()
        prefixes.add(tuple([]))
        for trace in self.language:
            prefix = []
            for activity in trace.get_seq():
                prefix.append(activity)
                if tuple(prefix) not in prefixes:
                    prefixes.add(tuple(prefix))
        return prefixes

    ####################################################################################################

    def get_max_activities(self):
        max_activities = dict(zip(self.activities, [0 for _ in range(len(self.activities))]))
        for trace in self.language:
            for activity in self.activities:
                count = trace.get_seq().count(activity)
                if count > max_activities[activity]:
                    max_activities[activity] = count
        return max_activities

    def gen_lha(self, pn):

        activities = self.get_activities()
        max_activities = self.get_max_activities()

        lha = "\n"

        # NbVariables
        lha += "NbVariables = " + str(3 + len([transition.get_name() for transition in pn.get_transitions()
                                               if transition.label is not None])) + ";\n"
        # NbLocations
        lha += "NbLocations = 3;\n\n"

        # Const
        lha += "const int n = " + str(len(activities)) + ";\n\n"

        # VariablesList
        lha += "VariablesList = {word,c,id"
        for transition in pn.get_transitions():
            if transition.get_label() is not None:
                lha += ",c_" + str(transition.get_name())
        lha += "};\n"
        # LocationsList
        lha += "LocationsList = {li,lfa,lfr};\n\n"

        # Variables
        lha += "PDF(Last(id), 1, 0, " + str(len(self.language)) + ");\n\n"

        # InitialLocations
        lha += "InitialLocations = {li};\n"
        # FinalLocations
        lha += "FinalLocations = {lfa};\n\n"

        # Locations
        lha += "Locations = {(li,TRUE);(lfa,(end=1));(lfr,TRUE);};\n\n"

        # Edges
        lha += "Edges = {\n"
        for transition in pn.get_transitions():
            if not transition.is_silent():
                lha += "((li,li),{" + str(transition.get_name()) + "},#,{word=word+" + str(
                    activities[transition.label]) + "*(n^c), c=c+1, c_" + str(transition.get_name()) + "=c_" + str(
                    transition.get_name()) + "+1});\n"
        if len(pn.get_transitions()) != len(activities):
            lha += "((li,li),{"
            for transition in pn.get_transitions():
                if transition.is_silent():
                    lha += str(transition.get_name()) + ","
            lha = lha[:-1] + "},#,#);\n"

        # Accepting Edges
        for trace in self.language[:-1]:
            lha += "((li,lfa),#,word=" + str(trace.map()) + ",{id=" + str(trace.get_id()) + "});\n"
        lha += "((li,lfa),#,#,{id=" + str(len(self.language) - 1) + "});\n"

        # Rejecting Edges
        for transition in pn.get_transitions():
            if not transition.is_silent():
                lha += ("((li,lfr),#,c_" + str(transition.get_name()) + ">=" +
                        str(max_activities[transition.get_label()] + 1) + ",#);\n")

        lha += "};\n"
        return lha
