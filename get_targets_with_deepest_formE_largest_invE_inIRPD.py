'''
Created on Jan 18, 2021

@author: jiadongc
'''
from pymatgen.analysis.interface_reactions import InterfacialReactivity
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter,CompoundPhaseDiagram
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.reaction_calculator import ComputedReaction
import os
import json
import itertools
import copy
from itertools import combinations
import pandas
import multiprocessing
def get_list_from_reactants_str(reactant_list_str):
    reactants = reactant_list_str[1:-1].split(", ")
    reactant_list = []
    for name in reactants:
        reactant_list.append(name[1:-1])
    return reactant_list
def make_stable_entry_from_comp(pd,comp):
    energy = pd.get_hull_energy(comp)
    new_entry=ComputedEntry(comp, energy-0.001)
    return new_entry

def get_lowest_entry_and_energy(cpd):
    depth = 0.1
    for e in cpd._stable_entries: #new_entries: not gona work
        ener = cpd.get_form_energy_per_atom(e)
#         print(e.name, ener)
        if ener < depth:
            depth = ener
            lowest_entry = e
#         print(depth,lowest_entry.name)
#         print()
    print(lowest_entry.name)
#     print("end")
    return lowest_entry,depth

def get_inverse_hull_energy(lowest_entry,cpd):
    mod_entries = []
    for e in cpd._stable_entries: #new_entries: not gona work
        if e.name != lowest_entry.name:
            mod_entries.append(e)
    mod_cpd = PhaseDiagram(mod_entries)
    invE = mod_cpd.get_decomp_and_e_above_hull(lowest_entry, allow_negative=True)[1]
    return invE

def remove_original_target_entry_if_exists_in_hull(product,entries):
    '''if the new composition exists in the convex hull, we will have a reaction
    like: ABCDE->ABCDE, and we do not want this. so we just remove the original 
    5-component entry if this happens.'''
    tar_in_hull =[]
    for entry in entries:
        try:
            reaction = ComputedReaction([entry], product)
            tar_in_hull.append(reaction)
        except:
            a = 1
    if len(tar_in_hull) != 0:
        for rea in tar_in_hull:
            for aaa in rea._reactant_entries:
                entries.remove(aaa)

def get_deepest_and_largest_invE_reaction(compstr):
    comp = Composition(compstr[0])
    els = [str(e) for e in comp.elements]
    entries = getOrigStableEntriesList(els)
#     for e in entries:
#         print(e.name)
#     print()
    pd = PhaseDiagram(entries)
    target_entry = make_stable_entry_from_comp(pd, comp)
    product = [target_entry]
    
    remove_original_target_entry_if_exists_in_hull(product,entries)
                
    combs = list(combinations(entries, 2))
    reactions_ener_inve = []
    reactions = []
    for reactants in combs:
        try:
            reactants = list(reactants)
            reaction = ComputedReaction(reactants, product)
            reactions.append(reaction)
#             print(reaction)
        except:
            a = 1
#             print(reactants[0].name,reactants[1].name)

    for reaction in reactions:
        print(reaction)
        reactants  = reaction._reactant_entries
        comp1 = reactants[0].composition
        comp2 = reactants[1].composition
        if reactants[0].name=="O2" or reactants[1].name == "O2":
            continue
        cricomps = pd.get_critical_compositions(comp1, comp2)

        new_entries = []
        for comp in cricomps:
            energy = pd.get_hull_energy(comp)
            new_entries.append(ComputedEntry(comp, energy))
        new_entries.append(target_entry)
        
        cpd = CompoundPhaseDiagram(new_entries,[comp1,comp2])
        lowest_entry, depth = get_lowest_entry_and_energy(cpd)
#         print(lowest_entry.name, depth)
#         print(target_entry.name)
        if lowest_entry.name == target_entry.name:
#             if depth <= -0.1:
            print(lowest_entry.name, target_entry.name, depth)
            invE = get_inverse_hull_energy(lowest_entry,cpd)
            reactions_ener_inve.append(
                [lowest_entry.name,[reactants[0].name,reactants[1].name], depth, invE]
                )
    if len(reactions_ener_inve) != 0:
        result = sorted(reactions_ener_inve,key = lambda x:x[3])[0]
        print(result)
        return result
    

if __name__ == "__main__":
#     with open("comp_and_neighbors_for_samsung","r") as f:
#         comp_neighbors = json.load(f)
#     print(len(comp_neighbors))
#     neighbors_num = 300
#     neighbor_comp = list(itertools.chain.from_iterable(comp_neighbors))
#     print(len(neighbor_comp))
#     neighbor_comp = sorted(neighbor_comp, key = lambda x: x[1], reverse = True)
#     neighbor_comp_ntargets = [x for x in neighbor_comp if x[1]>=neighbors_num]
#     print(len(neighbor_comp_ntargets))
#     for a in neighbor_comp_ntargets:
#         print(a)
#     del comp_neighbors
#     del neighbor_comp
    with open("mp_comp4s_Li.txt","r") as f:
        data = f.readlines()
    neighbor_comp_ntargets = []
    for e in data:
        if "Li" in e and "O" in e:
            neighbor_comp_ntargets.append(e.split("_"))

#     nome = ["Cl","F","S","H"]
#     for el in nome:
#         for nei in neighbor_comp_ntargets:
#             neis = [str(i) for i in Composition(nei[0]).elements]
#             if el in neis:
#                 neighbor_comp_ntargets.remove(nei)
    for e in neighbor_comp_ntargets:
        print(e[0])
    print(len(neighbor_comp_ntargets)) 
    
#     with multiprocessing.Pool(processes=12) as p:
#         results = p.map(get_deepest_and_largest_invE_reaction,neighbor_comp_ntargets)
#     results = list(filter(None, results))
#     print(len(results))
#     target_list = []
#     reactants_list = []
#     energy_list = []
#     inv_list = []
#     for re in results:
#         target_list.append(re[0])
#            
#         reactants_list.append(re[1])
#         energy_list.append(re[2])
#         inv_list.append(re[3])
#     df1 = pandas.DataFrame()
#     df1["Target"] = target_list
#     df1["Reactants"] = reactants_list
#     df1["Energy"] = energy_list
#     df1["invE"] = inv_list
#     print(len(energy_list))
#     df1.to_excel("mp_reaction_for_samsung_comp4_nonAA_n2.xlsx")
# 
# #     s =["LiAg3AlSbO6",0] 
# #     entries = get_deepest_and_largest_invE_reaction(s)
# #     print(entries)
#     
# 
#     
    
    
    
    