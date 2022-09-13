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
from myResearch.A201203_design_synthesis_way.interfacial_pdplotter import Inter_PDPlotter
import pandas
def get_list_from_reactants_str(reactant_list_str):
    reactants = reactant_list_str[1:-1].split(", ")
    reactant_list = []
    for name in reactants:
        reactant_list.append(name[1:-1])
    return reactant_list

# df = pandas.read_excel("reactions_list_using_16_most_common_precursors.xlsx",
#                        engine='openpyxl',
#                        index_col=0)
df = pandas.read_excel("mp_reaction_for_samsung_comp4_purify.xlsx",
                       engine='openpyxl',
                       index_col=0)

for i, row in df.iterrows():
    reactant = get_list_from_reactants_str(row["Reactants"])
    reactant1 = reactant[0]
    reactant2 = reactant[1]
    target = row["Target"]
#     if target != "Na2SnTeSO7":
#         continue

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
    print(cricomps)
    new_entries = []
    reactions_dict = {}
    for comp in cricomps:
        energy = pd.get_hull_energy(comp)
        new_entries.append(ComputedEntry(comp, energy))
        
        products = list(pd.get_decomposition(comp).keys())
        reaction = ComputedReaction(reactants, products)
        reactions_dict[ComputedEntry(comp, energy).name]=reaction
        print(reaction)
    tar_E = pd.get_hull_energy(target_comp)-0.001
    target_entry = ComputedEntry(target_comp,tar_E)
    new_entries.append(target_entry)
    reactions_dict[target_entry.name] = ComputedReaction(reactants, [target_entry])
    
    cpd = CompoundPhaseDiagram(new_entries,[comp1,comp2])
    print(len(reactions_dict))
    Inter_PDPlotter(cpd, reactions_dict, target_entry = target_entry).show(
        filename = "mpds_mp_samsung/"+target+".html"
        )


    