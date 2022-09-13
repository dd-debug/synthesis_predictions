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

TOLERANCE = 1e-6
def get_coeff(comp, reaction):
    coeff = reaction.get_coeff(comp)
    if abs(coeff) - round(coeff) < TOLERANCE:
        coeff = round(coeff)
#     print(coeff)
    return abs(coeff)

reactants_str=["Ba2Mn3","N2"]


els = ["Ba","Mn","N"]
entries = getOrigStableEntriesList(els)
pd = PhaseDiagram(entries)
# new_PDPlotter(pd).show()
comp1 = Composition(reactants_str[0])
comp2 = Composition(reactants_str[1])

cricomps = pd.get_critical_compositions(comp1, comp2)
new_entries = []
reactions_dict = {}
mols=[]
entry1 = ComputedEntry(comp1,pd.get_hull_energy(comp1))
entry2 = ComputedEntry(comp2,pd.get_hull_energy(comp2))
for comp in cricomps:
    energy = pd.get_hull_energy(comp)
    new_entries.append(ComputedEntry(comp, energy))
    
    products = list(pd.get_decomposition(comp).keys())
    reaction = ComputedReaction([entry1,entry2],products)
    a = get_coeff(comp1,reaction)
    b = get_coeff(comp2,reaction)
    mols.append(b/(a+b))
    print(a,b,b/(a+b))
    print(reaction)
    reactions_dict[ComputedEntry(comp, energy).name]=reaction
mols = mols[::-1]
print(mols)
cpd = CompoundPhaseDiagram(new_entries,[comp1,comp2])
Inter_PDPlotter(cpd, reactions_dict).show()
x,y = Inter_PDPlotter(cpd, reactions_dict).get_x_y_values()
    