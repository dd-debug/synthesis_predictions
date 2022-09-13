'''
Created on 2020.12.3
 
@author: jiadongc
'''
from pymatgen.analysis.reaction_calculator import ComputedReaction
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram
from itertools import combinations
TOLERANCE = 1e-6
def make_stable_entry_from_comp(pd,comp):
    energy = pd.get_hull_energy(comp)
    new_entry=ComputedEntry(comp, energy-0.001)
    return new_entry
 
def sepa_reacts_products(reaction):
    text = reaction.__str__()
    res = text.split(" -> ")[0].split(" + ")
    prs = text.split(" -> ")[1].split(" + ")
    return (res, prs)
     
def get_reactE_versus_comp(comp, reaction):
    coeff = reaction.get_coeff(comp)
    if abs(coeff) - round(coeff) < TOLERANCE:
        coeff = round(coeff)
    return reaction.calculated_reaction_energy / coeff
 
def get_reactions_with_num_reacts(product_comp, num = 3, num_reacts_strict = True):
    comp = product_comp
    els = [str(el) for el in comp.elements]
    entries = getOrigStableEntriesList(els)
    # build fake new entry
    pd = PhaseDiagram(entries)
    new_entry = make_stable_entry_from_comp(pd, comp)
     
    product = [new_entry]
    print(len(entries),"+ 1, total entries number")     
     
    combs = list(combinations(entries, num))
    reactions = []
    for reactants in combs:
        try:
            reactants = list(reactants)
            reaction = ComputedReaction(reactants, product)
            if reaction.get_coeff(comp) > TOLERANCE:
                if num_reacts_strict:
                    if len(sepa_reacts_products(reaction)[0]) == num:
                        if reaction not in reactions:
                            reactions.append(reaction)
                else:
                    if reaction not in reactions:
                            reactions.append(reaction)
        except:
            a = 1
    return reactions
 
def find_largest_reactE(product_comp, reactions):
    reactions = sorted(reactions, key = lambda re: get_reactE_versus_comp(product_comp, re))
    return reactions[0]
 
 
product_comp = Composition("MnZnCr2In2O8")
reactions = get_reactions_with_num_reacts(product_comp, num = 4, num_reacts_strict=True)
for a in reactions:
    print(a,",", a.calculated_reaction_energy)
last_step = find_largest_reactE(product_comp, reactions)
print("\n" + "Save most energies for the last step:")
print(last_step,",", get_reactE_versus_comp(product_comp, last_step))
 

