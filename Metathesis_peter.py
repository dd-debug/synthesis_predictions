'''
Created on Aug 3, 2019

@author: BoboMomo
'''

'''
Questions

1. How to determine what the nodes on the pseudo ternary graph are? Nodes can have more than 1 compound...?
2. How to recreate the psuedo ternary diagram?
3. What to do with the vector? How to input them into functions? Do we even need at all?

4. Warning for InterfacialReactivity, for recomposed nodes: is the energy value wrong? Need to correct it by using the terminal node's energy?

5. Remove byproducts in step 1 reaction? How to account for this?

6. 2-step reactions: need to make multiple phase diagrams for second step?

Notes
1. InterfacialReactivity might be affected by coefficients? Don't think so, but not sure yet.

TOMORROW: ADD REACTANT (AND PRODUCTS) COUNTERS TO NODE LIST

2. Regarding recomposition: if iterate more than twice, will need to subtract out repeated compounds (if there are any)


'''


 #%20Al - 30Nb - 10Ta - 30Ti - 10Zr
 
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.phase_diagram import PDEntry,PhaseDiagram,CompoundPhaseDiagram,PDPlotter,uniquelines
from pymatgen import Composition
from scipy.spatial import ConvexHull, HalfspaceIntersection
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.interface_reactions import InterfacialReactivity
from pymatgen.entries.computed_entries import ComputedEntry
import matplotlib.patches as patches
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import itertools
from matplotlib.pyplot import cm
import os 
import json
from pymatgen.ext.matproj import MPRester
import re
import itertools
import pandas
import random
import timeit
import plotly
import plotly.express as px
import plotly.graph_objects as go

##############################################################################
def getOrigAllEntriesList(els,filename = None, property_data = None):
    # if filename == None:
    #     filename = '-'.join(els)
    # cache = os.path.join(directory, filename)
    # if os.path.exists(cache):
    #     print('loading from cache.','-'.join(els))
    #     with open(cache, 'r') as f:
    #         dict_entries = json.load(f)
    #     list_entries = []
    #     for e in dict_entries:
    #         list_entries.append(ComputedEntry.from_dict(e))
    #     return list_entries
    # else:
    #print('Reading from database.')
    #print('-'.join(els))
    entries = MPR.get_entries_in_chemsys(els,property_data=property_data)

    dict_entries = []
    for e in entries:
        dict_entries.append(e.as_dict())
    # with open(cache,'w') as f:
    #     json.dump(dict_entries,f)
    return entries

# function to find total composition of a reaction (recomposition)
def ReComp(x):
    x = re.sub('\+','',x)
    x = re.split(' ',x)
    x = [e for e in x if e != '']
    
    # convert "small" numbers to 0
    # define small numbers as those with negative exponentials: containing 'e-'
    y = []
    for i in x:
        if 'e-' in i:
            y.append('0')
        else:
            y.append(i)
    
    comp = {'Compound':[],'Coeff':[]}
    
    #comp['Compound'] = [e for e in x if any(c.isalpha() for c in e)]
    count = 0
    for e in range(0,len(y)):
        if count == 0:  # automatically append first entry
            if any(c.isalpha() for c in y[e]):
                comp['Compound'].append(y[e])
                comp['Coeff'].append(1)     # if first entry is compound, append 1 to amount
                count += 1
            else:   # if first entry is a number, append it as the first coefficient
                comp['Coeff'].append(float(y[e]))
                count += 1
        else:   # for all other entries, observe previous entry
            # if current entry opposite previous entry for letter check, append appropriately
            if any(c.isalpha() for c in y[e]) != any(c.isalpha() for c in y[e-1]):
                if any(c.isalpha() for c in y[e]):
                    comp['Compound'].append(y[e])
                    count += 1
                else:
                    comp['Coeff'].append(float(y[e]))
            # otherwise, previous and current entries have same result for letter check
            # ==> both entries compounds (not numbers)
            # ==> append 1 to coefficient list
            else:
                comp['Compound'].append(y[e])
                comp['Coeff'].append(1)
                count += 1
    
    # then, put everything together
    # break each component down into its elements, then multiply by the coefficient
    comp_tot = Composition('')
    for i in range(0,len(comp['Compound'])):
        comp_tot += comp['Coeff'][i] * Composition(comp['Compound'][i])
    # then, output the sum in string form
    comp_tot = comp_tot.formula     # final string missing 1's for coefficients (problem?)
    return comp_tot

def NewReactants(coeff,reactants):
    # multiply coefficient across list of reactants
    #print(coeff)
    #print('Reactants: ' + reactants)
    
    reactants_split = re.split(' \+ ',reactants)    # split reactants by '+'
    reactants_split_split = [re.split(' ', x) for x in reactants_split] # split each reactant into coefficient and compound
    # loop over sublists of [coefficient,compound]
    # for each sublist, multiply coefficient by coefficient provided
    reactants_split_split_new = []
    for i in reactants_split_split:
        #print(i)
        i = [j for j in i if j != '']   # delete any blanks, which may have arisen from double spaces
        
        # check first entry: if contain non-number, combine w/ second entry, and set first entry to 1
        if len(i) > 1:
            first_number = [(not j.isdigit()) for j in i[0]]
            if any(first_number):
                i[1] = i[0]+i[1]
                i[0] = 1
        if len(i) == 1: # if only 1 entry, coefficient is 1 (implicit)
            i.append(i[0])
            i[0] = 1
        # multiply coefficient (first entry in sublist) by coefficient provided
        i[0] = str(round(coeff * float(i[0]),4))
        reactants_split_split_new.append(i)
        
    # now, put reactant string back together
    new_reactants_split_new = [' '.join(sublist) for sublist in reactants_split_split_new]
    new_reactants = ' + '.join(new_reactants_split_new)
    return(new_reactants)

def NewNodes(pairs,pairs_split,input_supersets,checked_pairs,Nodes,flag):
    start = timeit.default_timer()
    
    # initialize new nodes, and supersets of new nodes
    # should be 1:1
    new_precs = []
    supersets = []
    
    for pair in pairs:
        #######################################################################
        #print(pairs.index(pair))
        IR=InterfacialReactivity(Composition(pair[0]),Composition(pair[1]),pd)
        
        # decompose compounds in each pair (unnecessary for initial pairs)
        decomp0 = pairs_split[pairs.index(pair)][0]
        decomp1 = pairs_split[pairs.index(pair)][1]
        
        # initialize sublist of new pairs, for given set of parent nodes
        pair = [re.sub(' ','',i) for i in pair]
        new_precs_sub = [[pair[0],pair[1]]] # append reactant nodes to avoid repeats
        
        #######################################################################
        # find kinks for given pair: find kink w/ lowest energy
        K = []
        K_e = []
        
        for i in IR.get_kinks():
            #print(i)
            K.append(i)         # node
            K_e.append(i[4])    # energy
            
        # find index of kink w/ largest driving force
        k_ind = K_e.index(min(K_e))
        
        k = K[k_ind]
        #print(k[1])
        
        # count number of tie line crossing (number of new nodes discovered - 1)
        # but can't be less than 0 (in case no new nodes uncovered)
        n_k = max(0,len(K) - 2 - 1)
        
        if K_e[k_ind] <= 0:
            # check if node is terminal; already accounted for, so don't append to list
            
            # split reaction string
            reaction_split = re.split(' -> ', str(k[3]))
            
            #print(reaction_split[1])
            
            ###################################################################
            ###################################################################
            # 1. clean reaction string: 0 coefficients, leading product coefficient
            
            ###################################################################
            # 1.1 clean reaction strings of any compounds w/ ~0 coefficients
            reaction_split_split = [re.split(' \+ ',e) for e in reaction_split]
            #print(reaction_split_split)
            for e1 in reaction_split_split:
                for e2 in e1:
                    if 'e-' in e2:
                        e1.remove(e2)
            # put strings back together
            reaction_split = [' + '.join(e) for e in reaction_split_split]
            
            ###################################################################
            # 1.2 create original product string: add leading 1
            product_split = re.split(' \+ ',reaction_split[1])
            product_split_split = [re.split(' ',i) for i in product_split]
            #print(product_split_split)
            for i in product_split_split:
                # if type(j) == str:
                #     pass
                #     jlist = []
                #     jlist.append(j)
                #     jlist.insert(0,'1')
                #     j = jlist   
                if len(i) == 1:
                    i.insert(0,'1')   
                    pass
                    
            #print(product_split_split)
            product_split = [' '.join(i) for i in product_split_split]
            product_ones = ' + '.join(product_split)
            reaction_split[1] = product_ones
            #print(reaction_split[1])
            
            ###################################################################
            # 1.3 save final products after changes
            products_og = reaction_split[1]
            
            ###################################################################
            ###################################################################
            # 2. Append reaction to list
            
            ###################################################################
            # check if reaction and product strings equal
            # if so, terminal node ==> ignore (already added manually)
            # also, check if reaction energy ~ 0
            # if so, terminal node (whose reactants have been combined) ==> ignore
            # if not, append to list
            if (reaction_split[0] != reaction_split[1]) and (k[4] < -0.01):
                # first, recreate proper reactant string (with compounds instead of recompositions)
                # recreate it from the decompositions
                
                # convert reactant string to proper compositions
                # also, multiply their coefficients by the leading ratio
                # make sure to multiply by the correct coefficient!
                # pair[0]: coeff
                # pair[1]: 1-coeff
                
                ###############################################################
                # 2.1 find coefficients: grab directly from reactant string
                # first, see if either reactant has coefficient ~0
                
                compound_test = []
                
                reactants = re.sub('\+ ','',reaction_split[0])
                
                if k[1] < 0.001:
                    coeff_a = 0
                    # find coefficient of other reactant (might not be 1)
                    reactant_split = re.split(' ',reactants)
                    if len(reactant_split) == 1:
                        coeff_b = 1
                        compound_test.append(1)
                        compound_test.append(reactant_split[0])
                    else:
                        coeff_b = reactant_split[0]  # otherwise, grab first entry in list
                        compound_test.append(1)
                        compound_test.append(reactant_split[1])
                elif k[1] > 0.999:
                    coeff_b = 0
                    # find coefficient of other reactant (might not be 1)
                    reactant_split = re.split(' ',reactants)
                    if len(reactant_split) == 1:
                        coeff_a = 1
                        compound_test.append(0)
                        compound_test.append(reactant_split[0])
                    else:
                        coeff_a = reactant_split[0]  # otherwise, grab first entry in list
                        compound_test.append(0)
                        compound_test.append(reactant_split[1])
                else:
                    # find coefficient of other reactant (might not be 1)
                    reactant_split = re.split(' ',reactants)
                    if len(reactant_split) == 2:
                        coeff_a = 1
                        coeff_b = 1
                        compound_test.append(0)
                        compound_test.append(reactant_split[0])
                    elif len(reactant_split) == 4:
                        coeff_a = float(reactant_split[0])
                        coeff_b = float(reactant_split[2])
                        compound_test.append(0)
                        compound_test.append(reactant_split[1])
                    elif len(reactant_split) == 3:   # otherwise, figure out which compound has coefficient of 1
                        if any([j.isupper() for j in reactant_split[0]]):
                            coeff_a = 1
                            coeff_b = float(reactant_split[1])
                            compound_test.append(0)
                            compound_test.append(reactant_split[0])
                        else:
                            coeff_a = float(reactant_split[0])
                            coeff_b = 1
                            compound_test.append(0)
                            compound_test.append(reactant_split[1])
                
                #print(compound_test)
                
                ###############################################################
                # 2.2 determine which coefficients correspond to which reactant (no set order)
                # check which elements compound_test contains
                elem_test = Composition(compound_test[1]).elements
                elem_test_target = set(elem_test) & set(Composition(target).elements)
                
                # see which of the decompositions contains more of the elements present
                elem_decomp0 = Composition(pair[0]).elements
                elem_decomp1 = Composition(pair[1]).elements
                
                l_decomp0 = set(elem_decomp0) & set(Composition(target).elements)
                l_decomp1 = set(elem_decomp1) & set(Composition(target).elements)
                
                # pair the correct coefficients with the decompositions
                if (l_decomp0 == elem_test_target) and (compound_test[0] == 0):
                    coeff0 = coeff_a
                    coeff1 = coeff_b
                elif (l_decomp1 == elem_test_target) and (compound_test[0] == 1):
                    coeff0 = coeff_a
                    coeff1 = coeff_b
                else:
                    coeff0 = coeff_b
                    coeff1 = coeff_a
                 
                # print(elem_test_target)
                # print(l_decomp0)
                # print(l_decomp1)
                # print('\n')
                
                # print(coeff0)
                # print(coeff1)
                # print('\n')
                
                reactant1 = NewReactants(coeff0,decomp0)   # necessary?
                reactant2 = NewReactants(coeff1,decomp1)   # necessary?
                new_reactants = ' + '.join([reactant1,reactant2])
                
                # Determine if target synthesized by given reaction
                # first, create lists of reactants and products
                # take reactants and product strings: filter out '.',' ', and numbers
                # s_react = ''.join([e for e in new_reactants if e != '.' and e!= ' ' and not e.isdigit()])
                # s_prod = ''.join([e for e in reaction_split[1] if e != '.' and e!= ' ' and not e.isdigit()])
                # l_react = re.split(' \+ ',s_react)
                # l_prod = re.split(' \+ ',s_prod)
                
                # first delete '+ '
                s_react = re.sub('\+ ','',new_reactants)
                s_prod = re.sub('\+ ','',reaction_split[1])
                # then split on remaining spaces
                l_react = re.split(' ',s_react)
                l_prod = re.split(' ',s_prod)
                # then, remove entries without letters (check if contain any uppercase letter):
                l_react_red = []
                for i in l_react:
                    if any([j.isupper() for j in i]):
                        l_react_red.append(i)
                l_prod_red = []
                for i in l_prod:
                    if any([j.isupper() for j in i]):
                        l_prod_red.append(i)
                        
                #print(l_react_red)
                #print(l_prod_red)
                
                ###############################################################

                
                # save copy of original decomposed reaction
                decomp_reaction_og = new_reactants,reaction_split[1]
                
                # (only if target in products list?... NO), delete any reactants that also appear as products
                # first, subtract off the smaller of the two coefficients from both reactant and product
                # second, filter out the reactants and products that have 0 as a coefficient
                #if l_bool == 1:
                new_reactants_split = re.split(' \+ ',new_reactants)
                new_reactants_split_split = [re.split(' ',i) for i in new_reactants_split]
                for i in new_reactants_split_split:
                    #print(len(i) == 1)
                    if len(i) == 1:
                        i.insert(0,1)
                        
                #print('Reactants: ')
                #print(new_reactants_split_split)
                
                #print(reaction_split[1])
                new_products_split = re.split(' \+ ',reaction_split[1])
                new_products_split_split = [re.split(' ',i) for i in new_products_split]
                for i in new_products_split_split:
                    #print(len(i) == 1)
                    if len(i) == 1:
                        i.insert(0,1)
                
                #print('Products: ')
                #print(new_products_split_split)
                
                r_del = []
                p_del = []
                for r in new_reactants_split_split:
                    for p in new_products_split_split:
                        if r[1] == p[1]:    # check if any reactants and products match
                            
                            # subtract the minimum of their coefficients off from both compounds
                            m = min(float(r[0]),float(p[0]))
                            
                            r[0] = str(round(float(r[0]) - m,4))
                            p[0] = str(round(float(p[0]) - m,4))
                            
                for r in new_reactants_split_split:
                    if float(r[0]) < 0.001:
                        r_del.append(r)
                        
                        # adjust other dictionary entries:
                        # remove compounds from reactant and product lists
                        l_react_red.pop(l_react_red.index(r[1]))
                
                for p in new_products_split_split:
                    if float(p[0]) < 0.001:
                        p_del.append(p)
                    
                        # adjust other dictionary entries:
                        # remove compounds from reactant and product lists
                        l_prod_red.pop(l_prod_red.index(p[1])) # r and p should be same
                        
                        # remove compounds from total composition?
                        # not necessary: only iterating twice (sensible), and total composition just an intermediate
                        # if iterate more than twice, will need to adjust
                                
                #print(new_reactants_split_split)
                
                if (len(r_del) > 0) and (len(p_del) > 0):
                    # delete appropriate reactants and products
                    for i in r_del:
                        new_reactants_split_split.pop(new_reactants_split_split.index(i))
                    
                    for i in p_del:
                        new_products_split_split.pop(new_products_split_split.index(i))
                    
                    #print(new_reactants_red_split_split)
                    #print(new_products_red_split_split)
                    
                    # now, recombine reactants and products strings
                    
                    # re-initialize reactants and product strings
                    # use same names: names must match when drop out of if statement!
                    new_reactants_split = []
                    new_reactants = []
                    new_products_split = []
                    new_products = []
                    new_reactants_split = [' '.join(sublist) for sublist in new_reactants_split_split]
                    new_reactants = ' + '.join(new_reactants_split)
                    new_products_split = [' '.join(sublist) for sublist in new_products_split_split]
                    new_products = ' + '.join(new_products_split)
                    
                    reaction_split[1] = new_products
                    
                ###############################################################
                # Now check if target in products list
                # check compositions of products
                # see if equal to composition of target
                l_prod_comp = [Composition(e) for e in l_prod_red]
                l_bool = any([e == Composition(target) for e in l_prod_comp])
                    
                ###############################################################
                # Now check if reaction repeated (Note: check earlier?)
                rep_flag = 0
                
                # check if current reactant list matches reactant list of previous entries
                # same for product lists
                rep_ind = 0    # count number of reactions NOT repeated
                while rep_flag == 0 and rep_ind < len(Nodes['Reaction']):
                    if set(l_react_red) == set(Nodes['Reactant List'][rep_ind]):
                        if set(l_prod_red) == set(Nodes['Product List'][rep_ind]):
                            rep_flag = 1
                    rep_ind = rep_ind + 1
                    #print(rep_ind)

                ###############################################################
                    #new_reactants = ''
                    #reaction_split[1] = ''
                    # print(new_reactants_split_split)
                    # print(new_reactants)
                    # print(new_products_split_split)
                    # print(reaction_split[1])
                    
                l_prod_red = [re.sub(' ','',i) for i in l_prod_red]
                    
                if reaction_split[1] != '' and rep_flag == 0: # don't append empty reactions; causes issue for CompoundPhaseDiagram
                    Nodes['Target'].append(target)
                    Nodes['Precursors'].append(prec)
                    Nodes['Precursor Index'].append(count)
                    Nodes['Target Synthesized?'].append(l_bool)
                    Nodes['Tie Line Crossings'].append(n_k)
                    Nodes['Reactant Count'].append(len(l_react_red))
                    Nodes['Reactant List'].append(l_react_red)
                    Nodes['Product List'].append(l_prod_red)
                    Nodes['Reactants'].append(new_reactants)
                    Nodes['Products'].append(reaction_split[1])
                    Nodes['Reaction'].append(' -> '.join([new_reactants,reaction_split[1]]))
                    Nodes['Original Reaction'].append([k[3],decomp_reaction_og])
                    Nodes['Energy [eV/atom]'].append(round(k[2],4))
                    Nodes['Energy [kJ/mol]'].append(round(k[4],4))
                    Nodes['Parent Nodes'].append(pair)  # final, then initial [1,0]
                    # zero out ratio if very small
                    ratio = k[1]
                    if k[1] < 0.001:
                        ratio = 0
                    Nodes['Ratio'].append(round(ratio,4))
                    recomp = ReComp(reaction_split[1])
                    Nodes['Recomposition'].append(recomp)
                    
                    Nodes['Original Products'].append(products_og)
                    Nodes['Original Recomposition'].append(ReComp(products_og))
                    
                    # find indices of parent nodes
                    #Nodes['Parent Indices'].append([Nodes['Recomposition'].index(ReComp(e)) for e in input_supersets[pairs.index(pair)]])
                    
                    #new_precs_sub.append(reaction_split[1]) # use products for next pairs of reactants
                new_precs_sub.append(l_prod_red)
                    
                ###############################################################
                ###############################################################
                # # 3. Make lists of new reactants: one list for all compounds of current product node
                supersets_sub = new_precs_sub.copy()
                # # create list of new node pairs to check
                # # make lists of nodes for each set of parent nodes
                # if new_precs_sub != []: # omit sublist if empty; no new nodes for given set of parent nodes
                #     # first, loop through new_precs_sub
                #     # for each entry with multiple compounds, add subsets of the entry
                #     additions = []
                #     for i in new_precs_sub:
                #         node_split = re.split(' \+ ',i)
                #         if len(node_split) > 1:
                #             additions += node_split   # first, append on single compounds
                #             for i2 in range(0,len(node_split)):
                #                 supersets_sub.append(i) # append single compounds onto supersets as well
                #             # next, append on higher-order pairings
                #             # OMIT: ONLY WANT 1 REACTANT PER NODE (INCLUDING DAUGHTER NODES)
                #             # for j in range(2,len(node_split):
                #             #     subsets = list(itertools.combinations(node_split,j))
                #             #     for k in subsets:
                #             #         additions.append(' +  '.join(list(k)))
                #             #         supersets_sub.append(i)
                #     if additions != []:
                #         new_precs_sub += additions
                        
                #     # append sublist of current product node to full list of all product nodes
                #     new_precs.append(new_precs_sub)
                #     supersets.append(supersets_sub)
                new_precs.append(new_precs_sub)
                supersets.append(supersets_sub)
            
    ###########################################################################
    # 4. Create new pairs of reactants from list of new reactants
    if flag == 1:   # only make new pairs when requested (otherwise wasting time)
        # make new pairs of nodes to check
        # formation: combinations from 2 of the 'new_precs' sublists; one from each sublist
        
        #######################################################################
        # 4.1 make pairs of sublists (select 2 nodes)
        
        # remove all entries that did not uncover a new node
        new_precs_red = [i for i in new_precs if len(i[1]) > 0] # only if products present
        
        # make sets of all reactants and products
        react = []
        for i in new_precs_red:
            react = react + i[0]
        set_react = set(react)
        ls_react = list(set_react)
        
        prod = []
        for i in new_precs_red:
            prod = prod + i[1]
        set_prod = set(prod)
        ls_prod = list(set_prod)
        
        # remove compounds in products set that are in reactant set (no need to make them)
        # in doing so, prevent combos with 2 of same compound
        for k in ls_prod:
            if k in ls_react:
                ls_prod.remove(k)
            
        # make nodes from these two lists
        new_pairs_1 = [(i,j) for i in set_react for j in set_prod]
        #new_pairs_2r1p = [(i,j,k) for i in set_react for j in set_prod]
        
        # makes nodes from only products
        new_pairs_2 = list(itertools.combinations(ls_prod,2))
        #new_pairs_3p = list(itertools.combinations(ls_prod,3))
        
        new_pairs = new_pairs_1 + new_pairs_2
        new_pairs = list(set(new_pairs))  # remove repeats (should do nothing, since repeats removed from react and prod lists)
    
        # convert entries from tuples to lists
        for k in range(0,len(new_pairs)):
            new_pairs[k] = list(new_pairs[k])
            
            
        # add triplets to the list of new precursors
        # only for products
    
        # remove combinations between reactants and products from same tieline
        # not sure if necessary...
        #for k in new_precs_red:
        #    x = [[i,j] for i in k[0] for j in k[1]]
        #    for m in x:
        #        if m in new_pairs:
        #            new_pairs.remove(m)
                
        # remove combinations that mirror each other
        # remove combinations already checked in previous iterations
        for k in new_pairs:
            k_flip = [k[1],k[0]]
            if k == k_flip:
                new_pairs.remove(k)
            elif k in checked_pairs:
                new_pairs.remove(k)
            elif k_flip in checked_pairs:
                new_pairs.remove(k)
                
        # create supersets for new pairs
        new_supersets = [Composition(x[0]+x[1]).reduced_formula for x in new_pairs]
                
        # find indices of new pairs
        #new_pairs_ind = [new_pairs.index(i) for i in new_pairs]
        
        #######################################################################
        # # 4.1 make pairs of sublists (select 2 sides of the triangle)
                
        # # NOT SURE WHAT STILL NEEDED
                
        # combos = list(itertools.combinations(new_precs,2))
        # superset_combos = list(itertools.combinations(supersets,2))
        
        # new_pairs = []
        # new_pair_ind = [] # index of combo from which new pair was created
        # # then, make new pairs, one from each sublist
        # # remove pairs from the same list: only make pairs of compounds not already on the same tie line!
        # for combo in combos:
        # #new_pairs = new_pairs + [[a,b] for a in combo[0] and a not in combo[1] for b in combo[1] and b not in combo[0]]
        #     for a in combo[0]:
        #         for b in combo[1]:
        #             if a not in combo[1] and b not in combo[0]:
        #                 new_pairs = new_pairs + [[a,b]]
        #                 new_pair_ind.append(combos.index(combo))
        #                 # then, append combinations of subsets of reactants (could throw out byproducts)
        #                 # first, find subsets of each reactant: split on equal sign
        
        # #######################################################################
        # # 4.2 remove bad pairs:
        #     # already checked (original and flipped)
        #     # from same node
        #     # repeats
                        
        # pairs_flip = [[p[1],p[0]] for p in pairs]
        
        # # # initialize reduced list of new pairs
        # new_pairs_red = []
        # new_pairs_red_ind = []  # combo index of each pair kept
        # for p in new_pairs:
        #     if (p not in pairs) and (p not in pairs_flip) and (p[0] != p[1]) and (p not in new_pairs_red):
        #         new_pairs_red.append(p)
        #         new_pairs_red_ind.append(new_pairs.index(p))
        
        #######################################################################
        # 4.3 append on supersets: look up from 'combos'
        # find corresponding supersets for list of new nodes
        # supersets = []
        # for i in range(0,len(new_pairs)):
        #     # for each entry in a pair, find superset of index matching entry
        #     # use 'combos'
            
        #     # first, get combo index from 'new_pairs_red_ind'
        #     ind = new_pairs_ind[new_pairs_ind[i]]    # confusing: 2 indices
        #     superset0 = superset_combos[ind][0][combos[ind][0].index(new_pairs[i][0])]
        #     superset1 = superset_combos[ind][1][combos[ind][1].index(new_pairs[i][1])]
    
        #     supersets.append([superset0,superset1])
    else:
        new_pairs = []
        new_supersets = []
    
    # stop timer
    stop = timeit.default_timer()
    print('Time: ',stop-start)
    
    # append list of pairs checked onto list 'checked_pairs'
    for i in pairs:
        checked_pairs.append(i)

    return(new_pairs,new_supersets,checked_pairs,Nodes)

###############################################################################
def FastFE(compounds,FE):
    i = compounds
    if i in FE['Compound']:
        fe = FE['Formation Energy [eV/atom]'][FE['Compound'].index(i)]
    else:
        criteria = {'pretty_formula':i}
    
        props = ['pretty_formula','formation_energy_per_atom']
        entries = MPR.query(criteria=criteria, properties=props)
        
        fepa = min([i['formation_energy_per_atom'] for i in entries])
        
        # next, find number of atoms
        n_atm = Composition(i).num_atoms
        
        fe = n_atm*fepa
        
        # append new result to dictionary FE
        FE['Compound'].append(i)
        FE['Formation Energy [eV/atom]'].append(fe)
            
    return(fe,FE)

##############################################################################
##############################################################################
# I. Initialize: String of Elements --> Entries --> Phase Diagram
MPR = MPRester("3xIoDFdMEi4WQ2Xu")
#MPR = MPRester("2d5wyVmhDCpPMAkq")

Targets = ['MnZnCr2In2O8'] # GLOBAL VARIABLE!!! used in line: elem_test_target = set(elem_test) & set(Composition(target).elements)

form_e = -2.1498238196867976-0.001

# initialize dictionary of formation energies
FE = {'Compound':[],'Formation Energy [eV/atom]':[]}

for target in Targets:
    #target_corr = 'Ba4Mo12S18'
    path = 'C:/Users/jiadongc/Dropbox/WHSun_Lab/Code/Peter Son-Bell/3. Reaction Generation/'
    df_prec = pandas.read_excel(path + 'Precursors - ' + target + ', Reduced 1.xlsx')
    
    # initialize dictionary of nodes
    Nodes = {'Target':[],'Target Synthesized?':[],'Tie Line Crossings':[],'Precursor Index':[],'Reactant Count':[],'Precursors':[],'Reactant List':[],'Product List':[],'Reaction':[],'Original Reaction':[],'Original Products':[],'Original Recomposition':[],'Energy [eV/atom]':[],'Energy [kJ/mol]':[],'Reactants':[],'Products':[],'Parent Nodes':[],'Ratio':[],'Recomposition':[]}
    
    # Mg-Cr-S
    #target = 'MgCr2S4'
    #elem = ['Mg','Cr','S','Na','Cl']
    #prec = ['MgCl2','CrCl3','Na2S']
    
    target_elem = list(df_prec.columns)[1:]
    l_prec1 = [Composition(i).reduced_formula for i in df_prec[target_elem[0]] if i == i]
    l_prec2 = [Composition(i).reduced_formula for i in df_prec[target_elem[1]] if i == i]
    l_prec3 = [Composition(i).reduced_formula for i in df_prec[target_elem[2]] if i == i]
    l_prec4 = [Composition(i).reduced_formula for i in df_prec[target_elem[3]] if i == i]
    l_prec5 = [Composition(i).reduced_formula for i in df_prec[target_elem[4]] if i == i]
    
    # # delete spaces in compositions (just in case)
    l_prec1 = [re.sub(' ','',i) for i in l_prec1]
    l_prec2 = [re.sub(' ','',i) for i in l_prec2]
    l_prec3 = [re.sub(' ','',i) for i in l_prec3]
    l_prec4 = [re.sub(' ','',i) for i in l_prec4]
    l_prec5 = [re.sub(' ','',i) for i in l_prec5]
    
    l_prec = [[i,j,k,m,n] for i in l_prec1[:] for j in l_prec2[:] for k in l_prec3[:] for m in l_prec4[:] for n in l_prec5[:]]
    
    #l_prec = [target_elem]
    
    # quick check: only first 100 combinations
    #l_prec = l_prec[0:20]
    
    count = 0
    print("Precursor Count: ", count)
    
    for prec in l_prec:
        # initialize list of checked pairs of nodes (avoid repeat checks of same pairs)
        checked_pairs = []
        
        # create element space from precursor list
        prec_combine = ''.join([str(Composition(i)) for i in prec])
        prec_combine = ''.join([i for i in prec_combine if i != ' '])
        elem = re.split(r'[0-9]+',prec_combine)
        # remove blanks from list
        elem = [i for i in elem if i != '']
        # remove element repeats
        elem = list(dict.fromkeys(elem))
        
        prec_recomp = [ReComp(x) for x in prec]
        
        prec_pairs = []
        prec_pairs_no_space = []
        
        combos = list(itertools.combinations(prec_recomp,2))
        # convert all combinations from type 'tuple' to 'list'
        # then, append to list 'subsets'
        for i1 in combos:
            prec_pairs.append(list(i1))
            prec_pairs_no_space.append([''.join(e.split()) for e in list(i1)])
        
        #elem = ['Al','Ta','Ti','Cr','O']
        #prec = ['TaTi','O']
        
        directory = os.path.join("E:\dataset\entries")
        
        #print("Retrieving from MP")
        entries=getOrigAllEntriesList(elem,property_data=['pretty_formula','e_above_hull','formation_energy_per_atom'])
        # for e in entries:
        #     print(e.data['e_above_hull'])
        #print("Calculating")
        
        pd=PhaseDiagram(entries)
        
        # must run phase_diagram.py before every run
        # path: C:\Users\Work\Anaconda3\pkgs\pymatgen-2020.1.28-py37he980bc4_1\Lib\site-packages\pymatgen\analysis\phase_diagram.py
        target_entry = PhaseDiagram.make_entry_from_formation_energy(pd,Composition('MnZnCr2In2O8'),form_e)
        entries.append(target_entry)
        
        pd=PhaseDiagram(entries)
        
        # if len(elem) <= 4:
        #     PDPlotter(pd).get_plot()   # plot only 1-4 components
        
        # get hull energy at composition (example!)
        #print(pd.get_hull_energy(Composition('Al0.2Nb0.3Ta0.1Ti0.3Zr0.1')))

        #######################################################################
        # IV. Add precursors to node list
        for i in prec:
            # create reaction string
            reaction_split = ' '.join(['1.0000',i])
            reaction = (' -> '.join([reaction_split,reaction_split]))
            
            Nodes['Target'].append(target)
            Nodes['Target Synthesized?'].append(False)
            Nodes['Tie Line Crossings'].append(0)
            Nodes['Precursor Index'].append(count)
            #Nodes['Parent Indices'].append('N/A')
            Nodes['Reactant Count'].append(1)
            Nodes['Precursors'].append(prec)
            Nodes['Reactant List'].append(i)
            Nodes['Product List'].append(i)
            Nodes['Reactants'].append(i)
            Nodes['Products'].append(i)
            Nodes['Reaction'].append(reaction)
            Nodes['Original Reaction'].append(reaction)
            Nodes['Original Products'].append('1 ' + i)
            Nodes['Original Recomposition'].append(ReComp(i))
            Nodes['Energy [eV/atom]'].append(0)
            Nodes['Energy [kJ/mol]'].append(0)
            Nodes['Parent Nodes'].append('N/A')
            Nodes['Ratio'].append(1)
            Nodes['Recomposition'].append(ReComp(i))
        
        [new_pairs,supersets,checked_pairs,Nodes] = NewNodes(prec_pairs_no_space,prec_pairs_no_space,prec_pairs_no_space,checked_pairs,Nodes,1)
        comp_tot = [[ReComp(node) for node in pair] for pair in new_pairs]
        
        print("New Pairs: ", len(comp_tot))
        #while len(comp_tot) != 0:
        it_count = 0
        for it_count in range(0,10):
            if len(comp_tot) != 0:
                [new_pairs,supersets,checked_pairs,Nodes] = NewNodes(comp_tot,new_pairs,supersets,checked_pairs,Nodes,1)
                comp_tot = [[ReComp(node) for node in pair] for pair in new_pairs]
            it_count = it_count + 1
            print("New Pairs: ", len(comp_tot))

        print("Total Iterations: ", it_count)
             
             #(pairs,pairs_split,input_supersets,checked_pairs,Nodes,flag)
            
        # comp_tot2 = [[ReComp(node) for node in pair] for pair in new_pairs2]
        # if len(comp_tot2) != 0:
        #     [new_pairs3,supersets3,checked_pairs,Nodes] = NewNodes(comp_tot2,new_pairs2,supersets2,checked_pairs,Nodes,1)
            
        # comp_tot3 = [[ReComp(node) for node in pair] for pair in new_pairs3]
        # if len(comp_tot3) != 0:
        #     [new_pairs4,supersets4,checked_pairs,Nodes] = NewNodes(comp_tot3,new_pairs3,supersets3,checked_pairs,Nodes,1)
            
        # comp_tot4 = [[ReComp(node) for node in pair] for pair in new_pairs4]
        # if len(comp_tot4) != 0:
        #     [new_pairs5,supersets5,checked_pairs,Nodes] = NewNodes(comp_tot4,new_pairs4,supersets4,checked_pairs,Nodes,0)
        
        # comp_tot5 = [[ReComp(node) for node in pair] for pair in new_pairs5]
        # if len(comp_tot5) != 0:
        #     [new_pairs6,supersets6,checked_pairs,Nodes] = NewNodes(comp_tot5,new_pairs5,supersets5,checked_pairs,Nodes,0)
        
        
        
        # Create phase diagram
        # use function 'CompoundPhaseDiagram'
        # manually enter recompositions of nodes in dictionary
        
    #     # first, make list of compositions of recompositions
    # #    comps = [PDEntry(Composition(Nodes['Recomposition'][e]),Nodes['Energy [eV/atom]'][e],Nodes['Products'][e]) for e in range(0,len(Nodes['Reaction']))]
        #comps = [PDEntry(e.composition,pd.get_form_energy_per_atom(e)) for e in entries]
        
        # # find formation energies of each nodes
        # p_split = [re.split(' \+ ',Nodes['Original Products'][i]) for i in range(0,len(Nodes['Original Products']))]
        # p_split_split = [[re.split(' ',j) for j in i] for i in p_split]
        
        # # initialize list for formation energies
        # Nodes_fe = []
        
        # # query the compounds in the original products string
        # # use function FastFE
        # for i in p_split_split:
        #     fe_tot = 0
        #     for j in i:
        #         [fe,FE] = FastFE(j[1],FE)
        #         #print(fe)
        #         fe_tot += float(j[0])*fe
        #     Nodes_fe.append(fe_tot)
            
        # # check if nodes working
        # # pick random number from list of nodes: change to very negative formation energy
        # # rand_set = random.randint(0,len(Nodes_fe)-1)
        # # for i in range(0,rand_set):
        # #     rand_indx = random.randint(0,len(Nodes_fe)-1)
        # #     Nodes_fe[rand_indx] = -1000
        
        # # create list of PD Entries: use 'original' node, before having removed repeat terms
        # comps = [PDEntry(Nodes['Original Recomposition'][i],Nodes_fe[i],Nodes['Original Products'][i]) for i in range(0,len(Nodes['Recomposition']))]
        # precs = [Composition(e) for e in prec_recomp]
        # #precs = [Composition(e) for e in ['Mg','Cr','S']]
        # pdc = CompoundPhaseDiagram(comps,precs)
        # PDPlotter(pdc).get_plot()   # NOT WORKING!!!
    
        count += 1
        print('Precursor Count: ', count)
        
    
    # Output DataFrame
    if len(Nodes['Reaction']) != 0:
        df = pandas.DataFrame(Nodes)
        df.to_excel(path + 'Reactions 4.1 ' + target + ', all,all,all.xlsx')
        #df.to_excel('TEST.xlsx')

##############################################################################
# Plots Reactions
        
R = pandas.read_excel(path + 'Reactions, Final.xlsx')
R['Enthalpy per Target [kJ/mol]'] = -1*R['Enthalpy per Target [kJ/mol]']

bm_disc = R['Category'] == 'Discovered'
bm_ref = R['Category'] == 'Binary Mix or Paper'

R_disc = R[bm_disc]
R_ref = R[bm_ref]

# for reference reactions, filter out reactions w/ >2 reactants
bm_two = R_ref['Reactants'] > 2
R_ref = R_ref.drop(R_ref.index[bm_two])

# split reactions into reactants and products
# convert to compositions
R_clean = [re.sub(' ','',r) for r in R['Reaction'][bm_disc]]
R_clean = [re.sub(r'[0-9]+','',r) for r in R_clean]
R_split = [re.split('->',i) for i in R_clean]
R_react = [re.split('\+',j[0]) for j in R_split]
R_prod = [re.split('\+',j[1]) for j in R_split]

R_react = [[Composition(m) for m in n] for n in R_react]
R_prod = [[Composition(m) for m in n] for n in R_prod]

# screen reactions

# 1. no metal byproducts
bm_metal = []
gases = [Composition('H'), Composition('O'), Composition('N'), Composition('F'), Composition('Cl')]
for i in R_prod:
    for j in i:
        if len(j.elements) == 1 and j not in gases:
            bm_metal.append(R_prod.index(i))

bm_metal = list(set(bm_metal))

# 2. unrealistic compounds: manual filter
fake_comps = [Composition('CrB4')]
bm_fake_react = []
bm_fake_prod = []
for i in R_react:
    for j in i:
        if j in fake_comps:
            bm_fake_react.append(R_react.index(i))
for i in R_prod:
    for j in i:
        if j in fake_comps:
            bm_fake_prod.append(R_react.index(i))

# 3. 
            

bm_tot = bm_metal + bm_fake_react + bm_fake_prod
bm_tot = list(set(bm_tot))

R_disc = R_disc.drop(R_disc.index[bm_tot])

# plot
R_N_disc = R_disc[R_disc['Ternary'] == 'Nitride']
R_S_disc = R_disc[R_disc['Ternary'] == 'Sulfide']

R_N_ref = R_ref[R_ref['Ternary'] == 'Nitride']
R_S_ref = R_ref[R_ref['Ternary'] == 'Sulfide']

R_N = pandas.concat([R_N_disc, R_N_ref])
R_S = pandas.concat([R_S_disc, R_S_ref])

figN = px.scatter(R_N, x='Target', y='Enthalpy per Target [kJ/mol]', color='Category', 
                  title='Nitride Reactions', 
                  log_y=True, hover_data=['Reaction'], 
                  range_y=[1,5000],
                  color_discrete_sequence = ['red','blue']
                  ) #,'binaryreaction']) #symbol='sig_R_str') # range_y=[-20000,2000]
plotly.offline.plot(figN, filename='figNscreen.html')

figS = px.scatter(R_S, x='Target', y='Enthalpy per Target [kJ/mol]', color='Category', 
                  title='Sulfide Reactions', 
                  log_y=True, hover_data=['Reaction'], 
                  range_y=[0.2,20000],
                  color_discrete_sequence = ['red','blue']
                  ) #,'binaryreaction']) #symbol='sig_R_str') # range_y=[-20000,2000]
plotly.offline.plot(figS, filename='figSscreen.html')

##############################################################################
# VECTORS...?
# # get vectors on nodes on convex hull
# pd.all_entries_hulldata

# def get_critical_mu(kinks):
#     critical_mu=[]
#     critical_comps=[]
#     kinks=list(map(list, kinks))
#     for ii in range(0,len(kinks)-1):
#         m=(kinks[ii+1][2]-kinks[ii][2])/(kinks[ii+1][1]-kinks[ii][1])
#         b=kinks[ii][2]-m*kinks[ii][1]
#         critical_mu.append(b)
#         critical_comps.append(str(kinks[ii+1][3]).split(" -> ")[1])
#     return(critical_mu, critical_comps)

# critical_mu, critical_comps=get_critical_mu(IR.get_kinks())

# print(critical_mu)
# print(critical_comps)

# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111, aspect='equal')

# plt.rcParams.update({'font.size': 9})

# name = "Dark2"
# cmap = get_cmap(name)  # type: matplotlib.colors.ListedColormap
# colors = itertools.cycle(cmap.colors) # type: list
 
# for ii in range(0,len(critical_mu)-1):
#     c=next(colors)
#     y=0 #Location of bar
#     # mpatches.Rectangle((x,y), width, height, angle=0.0)   << Syntax for the bar
     
#     print(critical_comps[ii])
#     patch=patches.Rectangle((critical_mu[ii+1],y), critical_mu[ii]-critical_mu[ii+1], 1,facecolor=c,linewidth=1,edgecolor='k')
#     ax1.add_patch(patch)
#     #plt.axhline(critical_mu[ii+1],'--')
#     midpoint=critical_mu[ii]-(critical_mu[ii]-critical_mu[ii+1])/2
# #     plt.text(midpoint, y+1.1, critical_comps[ii], rotation=270,color=c,horizontalalignment='center',verticalalignment='bottom',wrap=True)
     
# #This is for the left-most patch, to negative infinity. Most reduced version
# patch=patches.Rectangle((min(critical_mu)-3,y), 3, 1,facecolor='w',linewidth=1,edgecolor='k')
# ax1.add_patch(patch)
# plt.text(min(critical_mu)-0.5, y+1.1, critical_comps[len(critical_mu)-1], rotation=270,color='k',horizontalalignment='center',verticalalignment='bottom')
 
# plt.xlim(min(critical_mu)-2,max(critical_mu)+2)
# plt.xlim(-14,1)
# plt.xticks(rotation=270)
 
# plt.ylim(0,20)
# plt.gca().axes.get_yaxis().set_visible(False)
# plt.show()
# #plt.savefig("MPEA_Ellingham.svg")