import pymatgen
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError
from pymatgen.ext.matproj import MPRester, MPRestError
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, CompoundPhaseDiagram, PhaseDiagramError
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
import matplotlib.pyplot as plt
import numpy as np
import json
from itertools import combinations
import matplotlib
import os
import pandas
import ast
import warnings

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()


MPR = MPRester("SZXJWLvi8njBGvA4sT")
directory = os.path.join("hull_dir")


def get_matproj_entries(els):
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
        try:
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
        except:
            newentries = []
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


def is_good_rxn(reactants, target, entries, pd, formEperatom, maxE=0):
    ''' If the number of stable phases on the compound convex hull with 
    the N reactants as terminal components is equal to N + 1, then there 
    are no intermediate phases.'''
    reactant_comps = [entry.composition for entry in reactants]
    energies = {entry.composition: pd.get_form_energy_per_atom(entry) for entry in reactants}
    energies[target.composition] = formEperatom
    balanced = False
    try:
        rxn = Reaction(reactant_comps, [target.composition])
        balanced = True
    except ReactionError:
        pass
    if balanced:
        if len(rxn.products) == 1 and rxn.products[0] == target.composition:
            E = rxn.calculate_energy(energies)
            if E <= maxE:
                try:
                    entries.append(target)
                    cpd = CompoundPhaseDiagram(entries, reactant_comps)
                    entries.remove(target)
                    if len(cpd.stable_entries) == (len(reactants) + 1):
                        for entry in cpd.stable_entries:
                            if entry.original_entry.composition == target.composition:
                                formE = cpd.get_form_energy_per_atom(entry)
                                return True, [list(reactants),
                                              [rxn.get_coeff(comp) for comp in reactant_comps],
                                              formE]
                    # entries.append(entry.original_entry)
                except (PhaseDiagramError, OverflowError) as exc:
                    pass
    return False, []


def is_reachable(reactant, AA_reactants, entries, pd):
    formEperatom = pd.get_form_energy_per_atom(reactant)
    if len(reactant.composition.elements) == 1:
        return True
    if reactant.name in [AA_reactant.name for AA_reactant in AA_reactants]:
        return True
    else:
        for combo_num in [2, 3]:
            for reactant_combo in combinations(AA_reactants, combo_num):
                good_rxn, rxn = is_good_rxn(reactant_combo, reactant, entries, pd, formEperatom, -0.1)
                if good_rxn:
                    return True
    return False

AA_reactants_all = pandas.read_csv("relevant_reactants.csv")
def get_best_reactants(target_entry, entries, pd, formEperatom):
    AA_reactants = []
    AA_reactants_names = []
    els = [el.symbol for el in target_entry.composition.elements]
    relevant_reactants = []
    for i, row in AA_reactants_all.iterrows():
        reactant = row["Reactant"]
        if isinstance(reactant, str):
            if all([el.symbol in els for el in Composition(reactant).elements]):
                relevant_reactants.append(reactant)
    for entry in pd.stable_entries:
        for reactant in relevant_reactants:
            if entry.composition.reduced_formula == Composition(reactant).reduced_formula and entry not in AA_reactants:
                AA_reactants.append(entry)
                AA_reactants_names.append(entry.name)
    combos = [2, 3]
    rxns = []
    if target_entry in entries:
        entries.remove(target_entry)
    entries_copy = entries*1
    for combo_num in combos:
        for reactant_combo in combinations(AA_reactants, combo_num):
            # if all([reactant.name in AA_reactants_names for reactant in reactant_combo]):
            good_rxn, rxn = is_good_rxn(reactant_combo, target_entry, entries, pd, formEperatom, -0.3)
            if good_rxn:
                # if all([is_reachable(reactant, AA_reactants, entries, pd) for reactant in reactant_combo]):
                rxns.append(rxn)
    if len(rxns) > 0:
        return sorted(rxns, key=lambda x:x[2])[0]
    else:
        return [[], [], 0.0]

if __name__ == "__main__":
    import itertools
    import json
    with open("comp_and_neighbors_for_prof_poudeu2","r") as f:
        comp_neighbors = json.load(f)
    # comp_neighbors = list(filter(None, comp_neighbors))
    neighbor_comp = list(itertools.chain.from_iterable(comp_neighbors))
    # neighbor_comp = list(filter(None, neighbor_comp))
    neighbor_comp = sorted(neighbor_comp, key = lambda x: x[1], reverse = True)
    df = pandas.DataFrame()
    reactants = []
    energies = []
    targets = []
    # rxns = pandas.read_csv("good_reactions.csv")
    # reactants = [ast.literal_eval(row["Reactants"]) for i, row in rxns.iterrows()]
    # energies = list(rxns["Energy"])
    # targets = list(rxns["Target"])
    # i = 2000
    i = 0
    old_len = 0
    skipped = []
    for e in neighbor_comp[:100000]:
        comp_str = e[0]
        print()
        print(comp_str)
        print("target # " + str(i))
        print(str(len(targets)) + " reactions found")
        print()
        # comp_str = 'BaMgSnSeO6'
        comp = Composition(comp_str)
        els = [el.symbol for el in comp.elements]
        entries = get_matproj_entries(els)
        if len(entries) > 0:
            pd = PhaseDiagram(entries)
            formEperatom = sum([pd.get_form_energy_per_atom(entry) * frac for entry, frac in pd.get_decomposition(comp).items()]) - 0.001
            target_entry = _make_entry_from_formEperatom(pd, comp, formEperatom)


            br = get_best_reactants(target_entry, entries, pd, formEperatom)
            if br[2] < -0.3:
                reactants.append([ent.name for ent in br[0]])
                energies.append(br[2])
                targets.append(comp_str)
        else:
            print("Skipping " + comp_str)
            skipped.append(comp_str)
        if i % 1000 == 0 and len(targets) > old_len:
            old_len = len(targets)
            df["Target"] = targets
            df["Reactants"] = reactants
            df["Energy"] = energies
            df.to_csv("good_reactions.csv")
            df = pandas.DataFrame()
        i += 1
    print(skipped)
    df = pandas.DataFrame()
    df["Target"] = targets
    df["Reactants"] = reactants
    df["Energy"] = energies
    df.to_csv("good_reactions.csv")

