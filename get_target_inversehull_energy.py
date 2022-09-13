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
from interfacial_pdplotter import Inter_PDPlotter
import json
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
df = pandas.read_excel("reactions_target_is_deepest_in_IRPD_for_samsung.xlsx",
                       engine='openpyxl',
                       index_col=0)
n = 0
inverse_hull_energies = []
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
        mod_entries = []
        for e in cpd._stable_entries: #new_entries: not gona work
            if e.name != target_entry.name:
                mod_entries.append(e)
            else:
                new_target_entry = e
        mod_cpd = PhaseDiagram(mod_entries)
        ener = mod_cpd.get_decomp_and_e_above_hull(new_target_entry, allow_negative=True)[1]
        print(new_target_entry.name,ener)
        inverse_hull_energies.append(ener)
df["Inverse_hull_energies"]=inverse_hull_energies
df.to_excel("reactions_target_is_deepest_in_IRPD_for_samsung.xlsx")
 
 
 
 
 
 

