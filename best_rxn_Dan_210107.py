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
import multiprocessing
def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()


MPR = MPRester("SZXJWLvi8njBGvA4sT")
# directory = os.path.join("hull_dir")
directory = os.path.join("C:/Users/jiadongc/dataset/stableEntries")

def get_matproj_entries(els):
    s_els = els.copy()
    s_els.sort()
    filename = '-'.join(s_els)
    cache = os.path.join(directory, filename)
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


def is_balanced_rxn(reactants, target, entries, pd, formEperatom):
    reactant_comps = [entry.composition for entry in reactants]
    try:
        rxn = Reaction(reactant_comps, [target.composition])
    except ReactionError:
        return False, 0.0
    if len(rxn.products) == 1 and rxn.products[0] == target.composition:
        target_reduced_comp, target_factor = target.composition.get_reduced_composition_and_factor()
        prod_E = formEperatom * target_factor * abs(rxn.get_coeff(target.composition)) * target_reduced_comp.num_atoms
        reactant_E = 0.0
        for reactant in reactants:
            reduced_comp, factor = reactant.composition.get_reduced_composition_and_factor()
            reactant_E += factor * abs(rxn.get_coeff(reactant.composition)) * pd.get_form_energy_per_atom(reactant) * reduced_comp.num_atoms
        return True, (prod_E - reactant_E) / target_reduced_comp.num_atoms
    else:
        return False, 0.0



def is_good_rxn(reactants, target, entries):
    reactant_comps = [entry.composition for entry in reactants]
    try:
        entries.append(target)
        cpd = CompoundPhaseDiagram(entries, reactant_comps)
        entries.remove(target)
        if len(cpd.stable_entries) == (len(reactants) + 1):
            return True
    except (PhaseDiagramError, OverflowError) as exc:
        pass
    return False


AA_reactants_all = pandas.read_csv("relevant_reactants.csv")
def get_best_reactants(target_entry, entries, pd, formEperatom, AA_only=True, combos=[2, 3], maxE=-0.001):
    AA_reactants = []
    AA_reactants_names = []
    els = [el.symbol for el in target_entry.composition.elements]
    relevant_reactants = []
    if AA_only:
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
    rxns = []
    if target_entry in entries:
        entries.remove(target_entry)
    entries_copy = entries*1
    reactants = entries
    if AA_only:
        reactants = AA_reactants
    for combo_num in combos:
        for reactant_combo in combinations(reactants, combo_num):
            # if all([reactant.name in AA_reactants_names for reactant in reactant_combo]):
            balanced, E = is_balanced_rxn(reactant_combo, target_entry, entries, pd, formEperatom)
            if balanced and E < maxE:
                rxns.append((list(reactant_combo), E))
    if len(rxns) > 0:
        for reactants, E in sorted(rxns, key=lambda x:x[1]):
            if is_good_rxn(reactants, target_entry, entries):
                return [reactants, E]
    return [[], 0.0]

if __name__ == "__main__":
    def search_reactions(filename, AA_only, combos=[2, 3], neighbors_num = 300, maxE=-0.001):
        import itertools
        import json
        with open("comp_and_neighbors_for_samsung","r") as f:
            comp_neighbors = json.load(f)
        # comp_neighbors = list(filter(None, comp_neighbors))
        neighbor_comp = list(itertools.chain.from_iterable(comp_neighbors))
        # neighbor_comp = list(filter(None, neighbor_comp))
        neighbor_comp = sorted(neighbor_comp, key = lambda x: x[1], reverse = True)
#         neighbor_comp_ntargets = neighbor_comp[:ntargets]
        neighbor_comp_ntargets = [x for x in neighbor_comp if x[1]>=neighbors_num]
        del comp_neighbors
        del neighbor_comp
        df = pandas.DataFrame()
        reactants = []
        energies = []
        targets = []
        # rxns = pandas.read_csv(filename)
        # reactants = [ast.literal_eval(row["Reactants"]) for i, row in rxns.iterrows()]
        # energies = list(rxns["Energy"])
        # targets = list(rxns["Target"])
        # i = 2000
        i = 0
        old_len = 0
        skipped = []
        for e in neighbor_comp_ntargets:
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
                br = get_best_reactants(target_entry, entries, pd, formEperatom, AA_only=AA_only, combos=combos, maxE=maxE)
                if br[1] < maxE:
                    reactants.append([ent.name for ent in br[0]])
                    energies.append(br[1])
                    targets.append(comp_str)
                    
            else:
                print("Skipping " + comp_str)
                skipped.append(comp_str)
            if i % 50 == 0 and len(targets) > old_len:
                old_len = len(targets)
                df["Target"] = targets
                df["Reactants"] = reactants
                df["Energy"] = energies
                df.to_csv(filename)
                df = pandas.DataFrame()
            i += 1
        print(skipped)
        df = pandas.DataFrame()
        df["Target"] = targets
        df["Reactants"] = reactants
        df["Energy"] = energies
        df.to_csv(filename)

    search_reactions("good_reaction_for_samsung_nonAA.csv", False, combos=[2])

