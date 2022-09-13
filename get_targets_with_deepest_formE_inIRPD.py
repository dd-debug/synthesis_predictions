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

import pandas
def get_list_from_reactants_str(reactant_list_str):
    reactants = reactant_list_str[1:-1].split(", ")
    reactant_list = []
    for name in reactants:
        reactant_list.append(name[1:-1])
    return reactant_list

'''
This code select reactions with reaction energy lower than -0.1 eV/atom
and the target 5-component material is the deepest entry in IR phase diagram.
'''
df = pandas.read_excel("good_reaction_for_samsung_nonAA_purify.xlsx",
                       engine='openpyxl',
                       index_col=0)
n = 0
target_list = []
reactants_list = []
energy_list = []
for i, row in df.iterrows():
    if float(row["Energy"]) < -0.1:
        reactant = get_list_from_reactants_str(row["Reactants"])
        reactant1 = reactant[0]
        reactant2 = reactant[1]
        target = row["Target"]


        target_comp = Composition(target)
        els = [str(e) for e in target_comp.elements]
        entries = getOrigStableEntriesList(els)
        pd = PhaseDiagram(entries)
        comp1 = Composition(reactant1)
        comp2 = Composition(reactant2)
        
        '''Compare to IR'''
        # IR = InterfacialReactivity(comp1,comp2,pd)
        # for i in IR.get_kinks():
        #     print(i)
        # print()
        
        for e in entries:
            if e.name == reactant1:
                entry1 = e
            if e.name == reactant2:
                entry2 = e
        
        reactants = [entry1,entry2]
        
        cricomps = pd.get_critical_compositions(comp1, comp2)
#         print(cricomps)
        new_entries = []

        for comp in cricomps:
            energy = pd.get_hull_energy(comp)
            new_entries.append(ComputedEntry(comp, energy))
            
        tar_E = pd.get_hull_energy(target_comp)-0.001
        target_entry = ComputedEntry(target_comp,tar_E)
        new_entries.append(target_entry)
        
        cpd = CompoundPhaseDiagram(new_entries,[comp1,comp2])
#         for e in cpd._stable_entries:
#             print(e.name,e.energy_per_atom)
#             print(cpd.get_form_energy_per_atom(e))

        depth = 0
        for e in cpd._stable_entries: #new_entries: not gona work
            ener = cpd.get_form_energy_per_atom(e)
            if ener < depth:
                depth = ener
                lowest_entry = e
            if target_entry.name == e.name:
                transformed_target_entry = e
        print("Target:",target_entry.name,target_entry.energy_per_atom)
        print("lowest_entry_name:",lowest_entry.name,lowest_entry.energy_per_atom)
        if depth == cpd.get_form_energy_per_atom(transformed_target_entry):
            print("find:",lowest_entry.name, target_entry.name)
            target_list.append(row["Target"])
            reactants_list.append(row["Reactants"])
            energy_list.append(row["Energy"])
df1 = pandas.DataFrame()
df1["Target"] = target_list
df1["Reactants"] = reactants_list
df1["Energy"] = energy_list
print(len(energy_list))
df1.to_excel("reactions_target_is_deepest_in_IRPD_for_samsung.xlsx")







#     Inter_PDPlotter(cpd, new_entries, reactions_list).show(filename = "reactions_htmls/"+target+".html")
