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
from myResearch.A200526_duality.convexhullpdplotter import new_PDPlotter
def get_list_from_reactants_str(reactant_list_str):
    reactants = reactant_list_str[1:-1].split(", ")
    reactant_list = []
    for name in reactants:
        reactant_list.append(name[1:-1])
    return reactant_list




reactants_str=["Na2O","MnO2"]
tarr = "Na2MnO3"
reactants_str=["LiPO3","ZnO"]
tarr = "LiZnPO4"
target_comp = Composition(tarr)
els = [str(el) for el in target_comp.elements]
entries = getOrigStableEntriesList(els)
pd = PhaseDiagram(entries)
# new_PDPlotter(pd).show()
comp1 = Composition(reactants_str[0])
comp2 = Composition(reactants_str[1])
'''normal compound phase diagram'''
# # tar_E = pd.get_hull_energy(target_comp)-0.001
# # target_entry = ComputedEntry(target_comp,tar_E)
# # entries.append(target_entry)
# import plotly.graph_objects as go
# cpd = CompoundPhaseDiagram(entries, [comp1,comp2])
# cpdp = PDPlotter(cpd)
# for coords, entry in cpdp.pd_plot_data[1].items():
#     if tarr == entry.name:
#         x_coord = coords[0]
#         y_coord = coords[1]
# fig = cpdp.get_plot()
# fig.add_trace(
#     go.Scatter(
#         mode='markers',
#         x=[x_coord],
#         y=[y_coord],
#         marker=dict(
#             color='rgba(135, 206, 250, 0.5)',
#             size=30,
#             line=dict(
#                 color='MediumPurple',
#                 width=4
#             )
#         ),
#         showlegend=False
#     )
# )
#  
# fig.show()

for e in entries:
    if e.name == reactants_str[0]:
        entry1 = e
    if e.name == reactants_str[1]:
        entry2 = e

reactants = [entry1,entry2]
# tar_E = pd.get_hull_energy(target_comp)-0.001
# target_entry = ComputedEntry(target_comp,tar_E)
# reaction = ComputedReaction(reactants, [target_entry])
# print(reaction)
# print(reaction.calculated_reaction_energy)

cricomps = pd.get_critical_compositions(comp1, comp2)
new_entries = []
reactions_dict = {}
for comp in cricomps:
    energy = pd.get_hull_energy(comp)
    new_entries.append(ComputedEntry(comp, energy))
    
    products = list(pd.get_decomposition(comp).keys())
    reaction = ComputedReaction(reactants, products)
    print(reaction)
    print(reaction.calculated_reaction_energy)
    reactions_dict[ComputedEntry(comp, energy).name]=reaction

    tar_E = pd.get_hull_energy(target_comp)-0.001
    target_entry = ComputedEntry(target_comp,tar_E)
    new_entries.append(target_entry)
    reactions_dict[target_entry.name] = ComputedReaction(reactants, [target_entry])
cpd = CompoundPhaseDiagram(new_entries,[comp1,comp2])
print(target_entry.name)
Inter_PDPlotter(cpd, reactions_dict,target_entry = target_entry).show()
    