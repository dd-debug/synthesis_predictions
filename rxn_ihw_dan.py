from pymatgen.entries import Entry
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError
import numpy as np
import matplotlib.pyplot as plt
from monty.json import MSONable
from functools import lru_cache
from copy import deepcopy
from matplotlib.path import Path
from collections import Counter
from adjustText import adjust_text
from itertools import combinations
from collections import defaultdict
from pymatgen.entries.computed_entries import ComputedEntry

def same_compound(comp1_str, comp2_str):
    """ Test if two composition strings are equivalent.

        Args:
            comp1_str (str): String of the first composition
            comp2_str (str): String of the second composition

        Returns:
            A boolean value, True if the strings correspond to the same composition
    """
    comp1 = Composition(comp1_str)
    comp2 = Composition(comp2_str)
    TF = []
    for el in comp1.elements:
        if el in comp2.elements:
            if comp1.get_atomic_fraction(el) != comp2.get_atomic_fraction(el):
                return False
        else:
            return False
    return True


def get_best_reactants(target_entry, entries, pd, formEperatom):
    combos = [2, 3]

    def get_best_reactant_and_energy(combo_num):
        rxns = []
        for reactant_combo in combinations(entries, combo_num):
            reactant_comps = [entry.composition for entry in reactant_combo]
            energies = {entry.composition: pd.get_form_energy_per_atom(entry) for entry in reactant_combo}
            energies[target_entry.composition] = formEperatom
            try:
                rxn = Reaction(reactant_comps, [target_entry.composition])
                if len(rxn.products) == 1 and rxn.products[0] == target_entry.composition:
                    E = rxn.calculate_energy(energies)
                    if E < -0.05:
                        # [print(entry.name + ': ' + str(round(pd.get_form_energy_per_atom(entry), 3))) for entry in reactant_combo]
                        # print(round(E, 3))
                        # print()
                        rxns.append([list(reactant_combo),
                                    [rxn.get_coeff(reactant.composition) for reactant in reactant_combo], 
                                    E])
            except ReactionError:
                pass
        if len(rxns) > 0:
            return sorted(rxns, key=lambda x:x[2])[0]
        else:
            return [[], [], 0.0]

    best_two_compound_rxn = get_best_reactant_and_energy(2)
    best_three_compound_rxn = get_best_reactant_and_energy(3)
    best = {}
    best_rxn = []
    if len(best_two_compound_rxn[0]) == 0:
        if len(best_three_compound_rxn[0]) > 0:
            best_rxn = best_three_compound_rxn
    elif len(best_two_compound_rxn[0]) > 0:
        best_rxn = best_two_compound_rxn
        if best_two_compound_rxn[2] > -0.2 and best_two_compound_rxn[2] < -0.4:
            best_rxn = best_three_compound_rxn
    if len(best_rxn) == 3:
        entries = best_rxn[0]
        fracs = best_rxn[1]
        best = {entries[i]: -1*fracs[i] for i in range(len(entries))}
    return best

    for combo_num in combos:
        for reactant_combo in combinations(entries, combo_num):
            reactant_comps = [entry.composition for entry in reactant_combo]
            energies = {entry.composition: pd.get_form_energy_per_atom(entry) for entry in reactant_combo}
            energies[target_entry.composition] = formEperatom
            try:
                rxn = Reaction(reactant_comps, [target_entry.composition])
                if len(rxn.products) == 1 and rxn.products[0] == target_entry.composition:
                    E = rxn.calculate_energy(energies)
                    if E < -0.05:
                        # [print(entry.name + ': ' + str(round(pd.get_form_energy_per_atom(entry), 3))) for entry in reactant_combo]
                        # print(round(E, 3))
                        # print()
                        rxns.append(CombinationRxn(list(reactant_combo), target_entry, E))
            except ReactionError:
                pass
    return rxns


class IHWEntry(Entry):
    """
    An object encompassing all relevant data for phase diagrams.

    .. attribute:: composition

        The composition associated with the IHWEntry.

    .. attribute:: formation_energy

        The formation energy associated with the entry.

    .. attribute:: inverse_hull_energy

        The inverse hull energy associated with the entry.

    .. attribute:: reactants

        Dictionary of {entry: fraction} of the entry's best reactants.

    .. attribute:: name

        String name of the entry. By default, this is the reduced formula for the composition, 
        but can be set to some other string for display purposes.

    .. attribute:: phase_type

        The entry's phase type, 'SS' for solid solutions or 'IM' for intermetallics
    """

    def __init__(
        self,
        composition: Composition,
        formation_energy: float,
        inverse_hull_energy: float,
        reactants: dict,
        name: str = None,
        phase_type: str = "IM",
    ):
        """
        Args:
            composition (Composition): Composition of the entry
            formation_energy (float): Formation energy of the entry
            inverse_hull_energy (float): Inverse hull energy of the entry
            reactants (dict): Dictionary of {entry: fraction} of the entry's best reactants
            name (str): String name of the entry
            phase_type (str): The entry's phase type, 'SS' or 'IM'
        """
        super().__init__(composition, formation_energy)
        self.name = name if name else self.composition.reduced_formula
        self.phase_type = phase_type
        self.inverse_hull_energy = inverse_hull_energy
        self.reactants = reactants
        self.formation_energy = formation_energy

    @property
    def energy(self) -> float:
        """
        :return: the energy of the entry.
        """
        # return self.formation_energy
        return self._energy + 0.0

    def __repr__(self):
        return "IHWEntry : {} with energy = {:.4f}".format(self.composition, self.formation_energy)

    def as_dict(self):
        """
        :return: MSONable dict.
        """
        return_dict = super().as_dict()
        return_dict.update({"name": self.name, "phase_type": self.phase_type, "inverse_hull_energy": self.inverse_hull_energy, "reactants": self.reactants})
        return return_dict

    def __eq__(self, other):
        # NOTE Scaled duplicates are not equal unless normalized separately
        if isinstance(other, self.__class__):
            return self.as_dict() == other.as_dict()
        return False

    def __hash__(self):
        # NOTE This hashing operation means that equivalent entries
        # hash to different values. This has implications on set equality.
        return id(self)

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: IHWEntry
        """
        return cls(
            Composition(d["composition"]),
            d["formation_energy"],
            d["inverse_hull_energy"],
            d["reactants"],
            d["name"] if "name" in d else None,
            d["attribute"] if "attribute" in d else None,
        )


class InverseHullWeb(MSONable):
    """
    Inverse Hull Web class. Inverse hull webs transform a chemical system's
    thermodynamic convex hull into a representation when phases are identified
    by their formation and inverse hull energies along with their best reactants.

    .. attribute: elements

        Elements in the inverse hull web.

    ..attribute: entries

        All entries provided for Inverse Hull Web construction. Note that this
        does not mean that all these entries are actually used in the inverse hull web. 
        For example, this includes the unstable intermetallic entries that are filtered 
        out before Inverse Hull Web construction.

    .. attribute: pd

        Phase Diagram used when calculating energies in the Inverse Hull Web

    .. attribute: target_comp

        Target composition

    .. attribute: target_comp_str

        String of target composition

    .. attribute: target_entry

        IHWEntry object of the target composition and phase

    .. attribute: stable_entries

        List of IHWEntry objects that are stable in the Phase Diagram corresponding to the Inverse Hull Web

    .. attribute: unstable_entries

        List of IHWEntry objects that are unstable in the Phase Diagram corresponding to the Inverse Hull Web

    .. attribute: web_entries

        List of IHWEntry objects that are actually used in the construction of the Inverse Hull Web
    """

    # Tolerance for determining if formation energy is positive.
    formation_energy_tol = 1e-11
    numerical_tol = 1e-8

    def __init__(
        self,
        entries: list,
        target_comp_str: str = None,
        target_phase: str = 'SS'
    ):
        """
        Constructor for inverse hull web.

        Args:
            entries ([PDEntry]): A list of PDEntry-like objects having an
                energy, energy_per_atom and composition.
            target_comp_str (str): String of the target composition. Defaults to the equimolar composition of all elements found in entries.
            target_phase (str): String of the target phase type, 'SS' or 'IM'. Defaults to 'SS'.
        """
        elements = set()
        for entry in entries:
            if "Phase Type" not in entry.data:
                entry.data["Phase Type"] = "IM"
            elements.update(entry.composition.elements)
        elements = sorted(list(elements))
        if target_comp_str is None:
            target_comp_str = ''.join([el.symbol for el in elements])
        target_comp = Composition(target_comp_str)
        pd = PhaseDiagram(entries)

        self.elements = elements
        self.entries = entries
        self.pd = pd
        self.target_comp = target_comp
        self.target_comp_str = target_comp_str
        self.target_entry = None
        for entry in pd.all_entries:
            if same_compound(self.target_comp_str, entry.name) and entry.data["Phase Type"] == target_phase:
                self.target_entry = entry
        invE_entries = list(pd.stable_entries)[::-1]
        if self.target_entry not in invE_entries:
            invE_entries.append(self.target_entry)
        invE_entries = entries[::-1]
        entry_dict = defaultdict(list)
        for entry in invE_entries:
            for i in range(len(entry.composition.elements), len(self.target_entry.composition.elements)):
                entry_dict[i].append(entry)
        self.stable_entries = [self.make_IHWEntry(entry, entry_dict[len(entry.composition.elements) - 1]) for entry in pd.stable_entries]
        self.unstable_entries = [self.make_IHWEntry(entry, entry_dict[len(entry.composition.elements) - 1]) for entry in pd.unstable_entries]
        web_entries = []
        for entry in entries:
            ie = self.make_IHWEntry(entry, entry_dict[len(entry.composition.elements) - 1])
            web_entries.append(ie)
        self.web_entries = web_entries

    def make_IHWEntry(self, entry, entry_list):
        """ Convert a PDEntry object to an IHWEntry.

            Args:
                entry (pymatgen.analysis.phase_diagram.PDEntry): Entry to convert
                entry_list (list of pymatgen.analysis.phase_diagram.PDEntry): 
                    List of entries to perform convex hull calculations against

            Returns:
                An IHWEntry object corresponding to entry
        """
        if isinstance(entry, IHWEntry):
            return entry
        phase_type = "IM"
        if "Phase Type" in entry.data:
            phase_type = entry.data["Phase Type"]
        if len(entry.composition.elements) > 1:
            formation_energy = self.pd.get_form_energy_per_atom(entry)
            try:
                entry_list.remove(entry)
            except:
                # print(entry.name + " not in entry_list")
                pass
            newpd = PhaseDiagram(entry_list)
            _, inverse_hull_energy = newpd.get_decomp_and_e_above_hull(entry,
                                                                       allow_negative=True)
            reactants = get_best_reactants(entry, entry_list, newpd, formation_energy)
            entry_list.append(entry)
        else:
            reactants = {}
            formation_energy = inverse_hull_energy = 0.0
        return IHWEntry(entry.composition,
                        formation_energy,
                        inverse_hull_energy,
                        reactants,
                        phase_type=phase_type)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        symbols = [el.symbol for el in self.elements]
        output = [
            "{} phase diagram".format("-".join(symbols)),
            "{} stable phases: ".format(len(self.stable_entries)),
            ", ".join([entry.name for entry in self.stable_entries]),
        ]
        return "\n".join(output)

    def as_dict(self):
        """
        :return: MSONAble dict
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "entries": [e.as_dict() for e in self.entries],
            "target_comp": self.target_comp,
        }

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: PhaseDiagram
        """
        entries = [MontyDecoder().process_decoded(dd) for dd in d["entries"]]
        return cls(entries)


class IHWPlotterEntry(IHWEntry):
    """
    An object encompassing all relevant data for phase diagrams.

    .. attribute:: alpha

        The alpha transparency of the entry when plotted

    .. attribute:: z

        The matplotlib zorder when plotting the entry

    .. attribute:: fontsize

         The size of the font when labeling the entry

    .. attribute:: marker

        Matplotlib marker symbol for the entry

    .. attribute:: markercolor

        RGB values to plot the entry with

    .. attribute:: markersize

         The markersize argument to pass to matplotlib when plotting the entry

    .. attribute:: textcolor

         The color of the entry's label

    .. attribute:: compdistance

         The euclidean distance between the entry's composition and the target

    .. attribute:: xy_offset

         The (x, y) offset values used to label the entry

    .. attribute:: annotate_str

        The string used to label the entry
    """

    def __init__(
        self,
        composition: Composition,
        formation_energy: float,
        inverse_hull_energy: float,
        reactants: dict,
        name: str,
        phase_type: str,
        alpha: float,
        z: int,
        fontsize: int,
        marker: object,
        markercolor: np.array,
        markersize: int,
        textcolor: str,
        compdistance: float,
        xy_offset: tuple,
        annotate_str: str
    ):
        """
        Args:
            composition (Composition): Composition of the entry
            formation_energy (float): Formation energy of the entry
            inverse_hull_energy (float): Inverse hull energy of the entry
            reactants (dict): Dictionary of {entry: fraction} of the entry's best reactants
            name (str): String name of the entry
            phase_type (str): The entry's phase type, 'SS' or 'IM'
            alpha (float): The alpha transparency of the entry when plotted
            z (int): The matplotlib zorder when plotting the entry
            fontsize (int): The size of the font when labeling the entry
            marker (str or tuple): Matplotlib marker symbol for the entry
            markercolor (np.array): RGB values to plot the entry with
            markersize (int): The markersize argument to pass to matplotlib when plotting the entry
            textcolor (str): The color of the entry's label
            compdistance (float): The euclidean distance between the entry's composition and the target
            xy_offset (tuple): The (x, y) offset values used to label the entry
            annotate_str (str): The string used to label the entry
        """
        super().__init__(composition, formation_energy, inverse_hull_energy,
                         reactants, name, phase_type)
        self.alpha = alpha
        self.z = z
        self.fontsize = fontsize
        self.marker = marker
        self.markercolor = markercolor
        self.markersize = markersize
        self.textcolor = textcolor
        self.compdistance = compdistance
        self.xy_offset = xy_offset
        self.annotate_str = annotate_str

    def __repr__(self):
        return "IHWPlotterEntry : {} with energy = {:.4f}".format(self.composition, self.formation_energy)

    def as_dict(self):
        """
        :return: MSONable dict.
        """
        return_dict = super().as_dict()
        return_dict.update({"alpha": self.alpha, "z": self.z, "fontsize": self.fontsize, "marker": self.marker, "markercolor": self.markercolor, "markersize": self.markersize, "textcolor": self.textcolor, "compdistance": self.compdistance, "xy_offset": self.xy_offset, "annotate_str": self.annotate_str})
        return return_dict

    def __eq__(self, other):
        # NOTE Scaled duplicates are not equal unless normalized separately
        if isinstance(other, self.__class__):
            return self.as_dict() == other.as_dict()
        return False

    def __hash__(self):
        # NOTE This hashing operation means that equivalent entries
        # hash to different values. This has implications on set equality.
        return id(self)



class IHWPlotter:
    """
    A class to plot Inverse Hull Webs

    .. attribute:: target_comp_str

        String name of the target composition

    .. attribute:: target_comp

        Composition object of the target phase

    .. attribute:: target_entry

        IHWEntry corresponding to target phase

    .. attribute:: inverse_hull_web

        InverseHullWeb to plot

    .. attribute:: elements

        List of names of elements in the chemical system

    .. attribute:: entries

        List of IHWPlotter entries to plot

    .. attribute:: x_bounds

        x bounds of the plot, used to position the color polygon

    .. attribute:: y_bounds

        y bounds of the plot, used to position the color polygon
    """

    def __init__(
        self,
        inverse_hull_web: InverseHullWeb,
        target_comp_str: str = None,
        target_phase: str = 'SS',
    ):
        """
        Constructor for inverse hull web plotter.

        Args:
            inverse_hull_web (InverseHullWeb): PInverseHullWeb object to plot.
            target_comp_str (str): String name of the target composition. Defaults to inverse_hull_web's target composition
            target_phase (str): String of the target phase type, 'SS' or 'IM'. Defaults to 'SS'.
        """
        if target_comp_str is None:
            target_comp_str = inverse_hull_web.target_comp_str
        self.target_comp_str = target_comp_str
        self.target_comp = Composition(self.target_comp_str)
        init_entries = inverse_hull_web.stable_entries
        self.inverse_hull_web = inverse_hull_web
        self.elements = [el.symbol for el in self.inverse_hull_web.elements]
        self.target_entry = None
        entries = []
        for entry in inverse_hull_web.web_entries:
            if same_compound(self.target_comp_str, entry.name) and entry.phase_type == target_phase:
                self.target_entry = entry
        if self.target_entry not in init_entries:
            init_entries.append(self.target_entry)
        for entry in init_entries:
            if entry is not None:
                alpha, z, fontsize = self.get_alpha_z_font_size(entry, init_entries)
                marker = self.get_plot_symbol(len(entry.composition.elements))
                entry_els = [el.symbol for el in entry.composition.elements]
                markercolor = self.get_marker_color(entry.composition)
                textcolor = self.get_text_color(entry)
                markersize = self.get_marker_size(entry)
                compdistance = self.get_compdistance(entry.composition)
                annotate_str = self.get_annotate_str(entry)
                xy_offset = self.get_xy_offset(entry.inverse_hull_energy,
                                               entry.formation_energy)
                entries.append(IHWPlotterEntry(entry.composition, entry.formation_energy,
                                               entry.inverse_hull_energy, entry.reactants,
                                               entry.name, entry.phase_type, alpha, z,
                                               fontsize, marker, markercolor, markersize,
                                               textcolor, compdistance, xy_offset, annotate_str))
        self.entries = entries
        entries_copy = entries*1
        for i in range(len(self.entries)):
            entry = self.entries[i]
            d = {}
            for reactant_entry, frac in entry.reactants.items():
                for sub_entry in entries_copy:
                    # if self.same_compound(decomp_entry.name, sub_entry.name) and decomp_entry.data["Phase Type"] == sub_entry.phase_type:
                    if same_compound(reactant_entry.name, sub_entry.name):
                        reactant_entry = sub_entry
                        d[sub_entry] = frac
            entry.reactants = d
            self.entries[i] = entry
            if same_compound(self.target_comp_str, entry.name) and entry.phase_type == 'SS':
                self.target_entry = entry
            self.labels = []
            self.x = []
            self.y = []
            self.objs = []

    def isequiatomic(self, comp):
        """ Returns True if a composition is equiatomic

            Args:
                comp (pymatgen.core.composition.Composition): Composition to check

            Returns:
                A boolean value
        """
        if isinstance(comp, str):
            comp = Composition(comp)
        if isinstance(comp, Composition):
            return len(set([comp.get_atomic_fraction(el) for el in comp])) == 1

    def get_plot_symbol(self, nels):
        """ Returns the matplotlib marker symbol corresponding 
            to the number of elements

            Args:
                nels (int): Number of elements in a composition

            Returns:
                A string or tuple corresponding to the matplotlib marker symbol
        """
        if nels == 1:
            return 'o'
        elif nels == 2:
            return 's'
        elif nels == 3:
            return '^'
        elif nels == 4:
            return 'd'
        elif nels == 5:
            return "p"
        elif nels == 6:
            return (5, 1)
        elif nels > 6:
            return (nels, 0)

    def _cot(self, a, b, c):
        return (np.dot(c - b, a - b) / np.linalg.norm(np.cross(c - b, a - b)))

    def compute_weights(self, point, vertices):
        """ Returns barycentric weights for each vertice corresponding to a point
            in a barycentric polygon defined by vertices.

            Meyer, M., Barr, A., Lee, H., & Desbrun, M. (2002)
            
            Args:
                point (list or tuple): Point in 2D space in polygon
                vertices (list): List of vertices defining a polygon

            Returns:
                A numpy.array representing the weights
        """
        point = np.asarray(point)
        total_weight = 0
        weights = []
        n = len(vertices)
        # This loops through each vertex and finds the corresponding weight
        for i in range(len(vertices)):
            previ = (i + n - 1) % n # index of previous vertex
            nexti = (i + 1) % n # index of next vertex
            if i == len(vertices) - 1:
                nexti = 0 # if the last element, index of next = 0

            # The numerator is the sum of the cotangents formed between the point of interest,
            # the vertex in the loop, and the previous / next vertices
            numerator = self._cot(point, vertices[i], vertices[previ]) + self._cot(point, vertices[i], vertices[nexti])
            # The denominator is the distance between the point of interest and the vertex in the loop
            denominator = np.linalg.norm(point - vertices[i])**2
            weights.append(numerator / denominator)
        return np.asarray(weights) / sum(weights)

    def get_rgb_from_weights(self, weights):
        """ Returns the rgb value corresponding to a set of element weights 

            Args:
                weights (list): List of element weights

            Returns:
                A numpy.array representing the rgb color
        """
        rgb = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0],
                        [1, 1, 0], [0, 0, 0], [0, 1, 1],
                        [1, 0, 1], [1, 1, 1], [0.7, 0.2, 0.2]])
        rgb = rgb[:len(weights)]
        # Multiply each weight by the corresponding rgb values and sum the array
        weights = np.asarray(weights).reshape(len(weights), 1)
        prod = np.multiply(weights, rgb)
        return np.sum(prod, axis=0)

    def get_rgb_from_point(self, point, vertices, nels):
        """ Returns the rgb value corresponding to a point in a color polygon 

            Args:
                point (list or tuple): Point in 2D space
                vertices (list of lists or list of tuples): List of 2D points defining polygon vertices
                nels (int): Number of elements to include in the color polygon

            Returns:
                A numpy.array representing the rgb color
        """
        rgb = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0],
                        [1, 1, 0], [0, 0, 0], [0, 1, 1],
                        [1, 0, 1], [1, 1, 1], [0.7, 0.2, 0.2]])
        rgb = rgb[:nels]

        # Get the weights corresponding to the point, multiply by the corresponding rgb values
        # and sum the array
        weights = self.compute_weights(point, vertices).reshape(len(vertices), 1)
        if nels == 2:
            new_wts = [0.0, 0.0]
            new_wts[0] = weights[0] + weights[1]
            new_wts[1] = weights[2] + weights[3]
            weights = np.array(new_wts)
        prod = np.multiply(weights, rgb)
        return np.sum(prod, axis=0)

    def get_marker_color(self, comp):
        """ Returns the color corresponding to a composition

            Args:
                comp (str or pymatgen.core.composition.Composition): composition to find color for

            Returns:
                A numpy.array representing the rgb color
        """
        el_list = self.elements
        weights = [0] * len(el_list)
        if isinstance(comp, str):
            comp = Composition(comp)
        if isinstance(comp, Composition):
            comp_els = [el.symbol for el in comp.elements]
            for i in range(len(el_list)):
                if el_list[i] in comp_els:
                    weights[i] = comp.get_atomic_fraction(el_list[i])
        return self.get_rgb_from_weights(weights)

    def get_text_color(self, entry):
        """ Returns blue for SS phases and red for IM phases.

            Args:
                entry (IHWPlotterEntry): Entry to get the text color for

            Returns:
                A string representing text color
        """
        if len(entry.composition.elements) == 1:
            return 'k'
        if entry.phase_type == "IM":
            return 'r'
        elif entry.phase_type == "SS":
            return 'b'
        else:
            return 'g'

    def get_marker_size(self, entry):
        """ Returns the marker size for a given entry

            Args:
                entry (IHWPlotterEntry): Entry to get the marker size for

            Returns:
                An int representing marker size
        """
        if entry == self.target_entry:
            return 300
        else:
            return 80

    def get_compdistance(self, comp):
        """ Returns the euclidean distance in composition space from comp and self.composition

            Args:
                comp (pymatgen.core.composition.Composition): Compositon to calculate distance for

            Returns:
                A float representing the distance between two compositions
        """

        comp_weights = [0] * len(self.elements)
        target_weights = [0] * len(self.elements)
        if isinstance(comp, str):
            comp = Composition(comp)
        if isinstance(comp, Composition):
            comp_els = [el.symbol for el in comp.elements]
            for i in range(len(self.elements)):
                if self.elements[i] in comp_els:
                    comp_weights[i] = comp.get_atomic_fraction(self.elements[i])
                target_weights[i] = self.target_comp.get_atomic_fraction(self.elements[i])
        # sqrt(delta_x^2 + delta_y^2 + delta_z^2...)
        return np.sqrt(sum([(comp_weights[i] - target_weights[i])**2 for i in range(len(self.elements))]))

    def get_arrow_xy(self, entry, reactant_entry):
        """ Returns two lists representing the x and y coordinates of the 
            head and tail of an arrow. Also returns the color of the arrow

            Args:
                entry (IHWPlotterEntry): Entry to find arrows for
                reactant_entry (IHWPlotterEntry): One of entry's best reactants

            Returns:
                A (list, list, numpy.array) tuple representing (x, y, color)
        """
        # points from reactant to entry, color of reactant
        x = np.asarray([reactant_entry.inverse_hull_energy, entry.inverse_hull_energy])
        y = np.asarray([reactant_entry.formation_energy, entry.formation_energy])
        arrow_color = reactant_entry.markercolor
        # if decomp_entry.compdistance < entry.compdistance:
        #     arrow_color = entry.markercolor
        #     x = [x[1], x[0]]
        #     y = [y[1], y[0]]
        return (x, y, arrow_color)

    def plot_reactants(self, ax, entry):
        """ Adds arrows connecting an entry to its best reactants

            Args:
                ax (matplotlib.pyplot.axes): Matplotlib axis object to plot on
                entry (IHWPlotterEntry): Entry to plot reactants for
        """
        reactants = {}
        comps = [e.composition.reduced_formula for e in entry.reactants.keys()]
        comp_count = Counter([e.composition.reduced_formula for e in entry.reactants.keys()])
        for comp, count in comp_count.items():
            tot_frac = 0.0
            lowest_E = 0.0
            stable_entry = None
            for e, frac in entry.reactants.items():
                if e.composition.reduced_formula == comp:
                    if count > 1:
                        tot_frac += frac
                        if e.formation_energy < lowest_E:
                            stable_entry = e
                            lowest_E = e.formation_energy
                    else:
                        stable_entry = e
                        tot_frac = frac
            reactants[stable_entry] = tot_frac
        tot = sum(frac for frac in reactants.values())

        for reactant_entry, phase_frac in reactants.items():
            frac = phase_frac / tot
            x, y, c = self.get_arrow_xy(entry, reactant_entry)
            alpha, z = 0.5, 0
            # print(reactant_entry.name + ': ' + str(frac))
            if entry.alpha + reactant_entry.alpha == 2:
                alpha, z = 1, 100
            if reactant_entry != self.target_entry:
                # This plots arrows connecting compounds
                pltx = np.asarray([x[0]], dtype=float)
                plty = np.asarray([y[0]], dtype=float)
                pltu = np.asarray([x[1] - x[0]], dtype=float)
                pltv = np.asarray([y[1] - y[0]], dtype=float)
                pltc = np.asarray([c])
                hw = 3
                hl = 5
                if frac < 0.2:
                    if frac < 0.01:
                        hw = 60
                        # hl = 100
                    else:
                        hw *= 0.2 / frac
                        # hl *= 0.2 / lw
                if pltu != 0 or pltv != 0:
                    self.objs.append(ax.quiver(pltx, plty, pltu, pltv,
                              units='xy', angles='xy',
                              color=c, scale=1, width=0.01 * frac,
                              scale_units='xy', edgecolor='k',
                              linewidth=0.15, headwidth=hw,
                              headlength=hl, alpha=alpha, zorder=z + 10))


    def get_reactant_tree(self, entry, reactant_dict, reactant_list, entries):
        """ Recursively finds and returns a list 
            of the reactant "tree" extending from entry
        """
        for sub_e, _ in reactant_dict.items():
            for e in entries:
                if same_compound(e.name, sub_e.name):
                    if e not in reactant_list:
                        reactant_list.append(e.name)
                    reactant_list = self.get_reactant_tree(e, e.reactants, reactant_list, entries)
        return reactant_list

    def get_alpha_z_font_size(self, entry, entries):
        """ Returns alpha transparency, zorder, and font size for an entry

            Args:
                entry (IHWPlotterEntry): Entry to get alpha transparency, zorder, and font size for

            Returns:
                A (float, int, int) tuple representing (alpha transparency, zorder, font size)
        """

        if self.target_entry is not None:
            in_reactants = False
            reactant_tree = self.get_reactant_tree(self.target_entry, self.target_entry.reactants, [], entries)
            print(reactant_tree)
            print()
            for reactant_entry in self.target_entry.reactants.keys():
                # if same_compound(entry.name, reactant_entry.name):
                if entry.name in reactant_tree:
                    in_reactants = True
            if entry == self.target_entry or in_reactants:
                return (1.0, 100, 8)
        return (0.5, 0, 6)

    def get_annotate_str(self, entry):
        """ Returns a string corresponding to an entry's name

            Args:
                entry (IHWPlotterEntry): Entry to get string for

            Returns:
                A string corresponding to the entry's name
        """
        if entry.phase_type == "IM":
            return entry.name
        else:
            if self.isequiatomic(entry.name):
                return ''.join([el.symbol for el in entry.composition.elements])
            else:
                return ''.join([el.symbol + str(round(comp.get_atomic_fraction(el), 2)) for el in entry.composition.elements])

    def get_xy_offset(self, invE, formE):
        """ Returns offset for text placement based on a point's location

            Args:
                invE (float): Inverse hull energy
                formE (float): Formation energy

            Returns:
                A (float, float) tuple representing the (x, y) offset
        """
        if abs(invE) > abs(0.5 * formE):
            return (-5, +5)
        else:
            return (+5, -5)

    def _regular_n_ngon(self, n, r=1):
        """ For n > 2, returns vertices for a regular polygon with vertices n scaled by r.
            For n = 2, returns a rectangle.
            
            Args:
                n (int): number of vertices
                r (float): factor to scale shape size by
        """
        if n == 2:
            x = r * 1.25
            y = r / 2
            return np.array([[-x, y], [-x, -y],
                             [x, -y], [x, y]])
        else:
            thetas = 2 * np.pi * np.asarray(range(n)) / n
            if n > 4:
                # If not a square or triangle, align with y axis
                diffs = thetas - (np.pi / 2)
                thetas -= diffs[np.argmin(np.abs(diffs))] + np.pi
            return np.column_stack((r * np.cos(thetas), r * np.sin(thetas)))

    @lru_cache(maxsize=None)
    def _plot_poly(self):
        """ Helper function for plot_color_polygon, plots the actual points in the 
            color polygon
        """
        n_pts = 100
        nverts = len(self.elements)
        aspect_yx = abs((self.y_bounds[0] - self.y_bounds[1]) / (self.x_bounds[0] - self.x_bounds[1]))
        mid_pt = (self.x_bounds[0] / 2, self.y_bounds[0] / 4)
        r = max(mid_pt) * 0.2
        # Set polygon vertices based on mid_pt
        n_gon_verts = self._regular_n_ngon(nverts, r=r) + mid_pt
        n_gon_verts[:, 0] += mid_pt[0] / 2
        n_gon_verts[:, 1] -= 3 * mid_pt[1] / 4
        n_gon_verts[:, 1] *= aspect_yx
        sample_rangex = np.linspace(-2 * r, 2 * r, n_pts) + 1.5 * mid_pt[0]
        sample_rangey = np.linspace(-2 * r, 2 * r, n_pts) + mid_pt[1] / 4
        sample_rangey *= aspect_yx
        initial_mask = np.meshgrid(sample_rangex, sample_rangey) # make a canvas with coordinates
        x, y = np.meshgrid(sample_rangex, sample_rangey) # make a canvas with coordinates
        x, y = x.flatten(), y.flatten()
        points = np.vstack((x,y)).T 

        p = Path(n_gon_verts) # make a polygon
        grid = p.contains_points(points)
        mask_TF = grid.reshape(n_pts, n_pts) # mask with points inside a polygon
        in_polyx = initial_mask[0][mask_TF]
        in_polyy = initial_mask[1][mask_TF]
        in_polyx, in_polyy = in_polyx.flatten(), in_polyy.flatten()
        in_poly_points = np.vstack((in_polyx,in_polyy)).T

        # Get color corresponding to each point
        c_list = []
        for i in range(len(in_poly_points)):
            x = in_poly_points[i, 0]
            y = in_poly_points[i, 1]
            c_list.append(self.get_rgb_from_point((x, y), n_gon_verts, nels=nverts))
        return in_poly_points, c_list


    def plot_color_polygon(self, fig, ax):
        """ Adds a color polygon composition legend to a matplotlib plot.
            
            Args:
                fig (matplotlib.pyplot.figure): Figure to add polygon to
                ax (matplotlib.pyplot.axes): Axis to add polygon to
        """
        nverts = len(self.elements)
        labels = self.elements
        aspect_yx = abs((self.y_bounds[0] - self.y_bounds[1]) / (self.x_bounds[0] - self.x_bounds[1]))
        mid_pt = (self.x_bounds[0] / 2, self.y_bounds[0] / 4)

        x_offset = 0
        current_ax = False
        # Set 'radius' of polygon based on mid_pt
        r = max(mid_pt) * 0.2
        n_gon_verts = self._regular_n_ngon(nverts, r=r)  #+ np.asarray(mid_pt)
        # Get a list of points within the polygon and their corresponding colors
        in_poly_points, c_list = self._plot_poly()

        # Plot each of these points
        for i in range(len(in_poly_points)):
            x = in_poly_points[i, 0]
            y = in_poly_points[i, 1]
            c = c_list[i]
            ax.scatter(x, y, color=c)

        # Add labels to the vertices
        for i in range(nverts):
            j = i
            if nverts == 2 and i == 1:
                j = 3
            vert = (self._regular_n_ngon(nverts, r=r) * 1.25)[j] + np.asarray(mid_pt)
            vert[0] += mid_pt[0] / 2
            vert[1] -= 3 * mid_pt[1] / 4
            vert[1] *= aspect_yx
            if nverts == 2:
                vert[1] += abs(r)
            ax.text(vert[0], vert[1], self.elements[i], fontsize=8, color='k', va='center', ha='center')
        ax.set_xlim(self.x_bounds)
        ax.set_ylim(self.y_bounds)

    def plot_entry(self, ax, entry):
        """ Adds a single entry to the inverse hull web

            Args:
                ax (matplotlib.pyplot.axes): Axis object to plot on
                entry (IHWPlotterEntry): Entry to plot on inverse hull web
        """
        self.objs.append(ax.scatter(entry.inverse_hull_energy, entry.formation_energy,
                   s=entry.markersize, color=entry.markercolor,
                   marker=entry.marker, alpha=entry.alpha, zorder=entry.z + 10))
        if len(entry.composition.elements) > 1:
            # ax.annotate(entry.annotate_str, xy=(entry.inverse_hull_energy, entry.formation_energy),
            #             xytext=entry.xy_offset, zorder=entry.z + 20, textcoords='offset points',
            #             fontsize=entry.fontsize, color=entry.textcolor, alpha=entry.alpha)
            self.labels.append(ax.text(entry.inverse_hull_energy, entry.formation_energy, entry.annotate_str, zorder=entry.z + 20, fontsize=entry.fontsize, color=entry.textcolor, alpha=1))
            self.x.append(entry.inverse_hull_energy)
            self.y.append(entry.formation_energy)
        self.plot_reactants(ax, entry)

    def plot(self, fig=None, ax=None, plot_color_polygon=True):
        """ Plots the inverse hull web
            
            Args:
                fig (matplotlib.pyplot.figure): Matplotlib figure to plot on. Will create new figure if None is passed.
                ax (matplotlib.pyplot.axes): Matplotlib axes to plot on. Defaults to fig.gca() when None is passed.
                plot_color_polygon (bool): Whether or not to add the color polygon legend.
        """
        if fig is None:
            fig, ax = plt.subplots()
        if ax is None:
            ax = plt.gca()
        for entry in self.entries:
            self.plot_entry(ax, entry)
        self.x_bounds = ax.get_xlim()
        self.y_bounds = ax.get_ylim()
        if plot_color_polygon:
            self.plot_color_polygon(fig, ax)
        ax.set_title(self.target_comp_str, fontsize=16)
        ax.set_xlabel('Inverse Hull Energy (eV/atom)', fontsize=14)
        ax.set_ylabel('Formation Energy (eV/atom)', fontsize=14)
        # adjust_text(self.labels, x=self.x, y=self.y, add_objects=self.objs)
        # adjust_text(self.labels, x=self.x, y=self.y)
        return fig

    def show(self, **kwargs):
        r"""
        Draw the phase diagram using Matplotlib and show it.

        Args:
            **kwargs: Arguments passed to plot.
        """
        self.plot(**kwargs)
        plt.show()

    def write_image(self, stream, image_format="svg", **kwargs):
        r"""
        Writes the phase diagram to an image in a stream.

        Args:
            stream:
                stream to write to. Can be a file stream or a StringIO stream.
            image_format (string): format for image. Can be any of matplotlib supported formats. Defaults to svg for best results for vector graphics.
            **kwargs: Pass through to plot function.
        """
        f = self.plot(**kwargs)
        f.set_size_inches((12, 10))
        plt.savefig(stream, format=image_format)

if __name__ == "__main__":
    import os
    import json
    from pymatgen.ext.matproj import MPRester

    MPR = MPRester("SZXJWLvi8njBGvA4sT")
    directory = os.path.join("C:/Users/jiadongc/dataset/stableEntries")

    def get_matproj_entries(els):
        els.sort()
        filename = '-'.join(els)
        cache = os.path.join(directory, '-'.join(els))
        if os.path.exists(cache):
            print(cache)
            print('loading from cache.','-'.join(els))
            with open(cache, 'r') as f:
                dict_entries = json.load(f)
            newentries = []
            for e in dict_entries:
                newentries.append(ComputedEntry.from_dict(e))
        else:
            print('Reading from database.')
            print('-'.join(els))
            entries = MPR.get_entries_in_chemsys(els)
            pd = PhaseDiagram(entries)
            newentries=[]
            for e in pd.stable_entries:
                newentries.append(e)
            dict_entries = []
            for e in newentries:
                dict_entries.append(e.as_dict())
            with open(cache,'w') as f:
                json.dump(dict_entries,f)
        for entry in newentries:
            entry.data["Phase Type"] = "IM"  # We will use this later to color solid solution and intermetallic phases
        return newentries


    def _make_entry_from_formEperatom(pd, c, formEperatom):
        # Revert formation energies to total energies
        EntryE = formEperatom*c.num_atoms+ sum([c[el]*pd.el_refs[el].energy_per_atom
                                       for el in c.elements]) 
        new_entry = ComputedEntry(c, EntryE)
        return new_entry


    def get_entries_from_target(comp):
        els = [el.symbol for el in comp.elements]
        entries = get_matproj_entries(els)
        pd = PhaseDiagram(entries)
        formEperatom = sum([pd.get_form_energy_per_atom(entry) * frac for entry, frac in pd.get_decomposition(comp).items()])
        entries.append(_make_entry_from_formEperatom(pd, comp, formEperatom))
        return entries

    comp_str = 'Cr2In2MnZnO8'
    comp = Composition(comp_str)
    els = [el.symbol for el in comp.elements]
    entries = get_entries_from_target(comp)
    pd = PhaseDiagram(entries)
    ihw = InverseHullWeb(entries, target_comp_str=comp_str, target_phase='IM')
    plotter = IHWPlotter(ihw, target_comp_str=comp_str, target_phase='IM')
    plotter.show(plot_color_polygon=False)
