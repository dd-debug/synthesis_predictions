'''
Created on Jan 6, 2021

@author: jiadongc
'''
import pandas
import collections
import numpy as np
from chord import Chord
import plotly.graph_objects as go

Chord.user = "jdchen1997@outlook.com"
Chord.key = "CP-c4615f36-9908-43f4-bdbe-97ff1b9f3fbd"
def get_list_from_reactants_str(reactant_list_str):
    reactants = reactant_list_str[1:-1].split(", ")
    reactant_list = []
    for name in reactants:
        reactant_list.append(name[1:-1])
    return reactant_list


# good_reactions_AA = pandas.read_csv("good_reactions_AA_purify.csv")
filename = "reactions_target_is_deepest_in_IRPD"
good_reactions_AA = pandas.read_excel(filename +".xlsx",
                       engine='openpyxl',
                       index_col=0)

def count_precursors_frequency_and_get_all_reactions(good_reactions_AA):
    reactants = []

    targets = []
    energy = []
    all_reactions = []
    for i, row in good_reactions_AA.iterrows():
        if float(row["Energy"]) < -0.1:
            reactant = get_list_from_reactants_str(row["Reactants"])
            reactants += reactant
            all_reactions.append(reactant)
            targets.append(row["Target"])
            energy.append(row["Energy"])
    
    print(len(all_reactions), "out of %s reactions with reactE < -0.1" % (len(good_reactions_AA["Reactants"])))
    ctr1 = collections.Counter(reactants)
    print(ctr1)
    return ctr1, all_reactions, targets, energy

def num_comp5_with_most_common_precursors(n):
    ctr1, all_reactions = count_precursors_frequency_and_get_all_reactions(good_reactions_AA)[0:2]
    rank_list = ctr1.most_common()[0:n]
    print(rank_list)
    precursors = []
    for pc in rank_list:
        precursors.append(pc[0])
    print("We choose most common %s" % n,"precursors")

    board = []
    for i in range(n):
        board.append([0 for j in range(n)])
        
    for pc in precursors:
        for reactants in all_reactions:
            if pc in reactants:
                rs_copy = reactants.copy()
                rs_copy.remove(pc)
                if rs_copy[0] in precursors:
                    board[precursors.index(pc)][precursors.index(rs_copy[0])] += 1 
    for bo in board:
        print(bo)
    numofcomp5 = np.count_nonzero(board)/2
    return numofcomp5

def get_pareto_frontier_plot():
    yy = []
    
    xx = []
    textlist=[]
    for n in range(1,33):
        numofcomp5 = num_comp5_with_most_common_precursors(n)
        yy.append(numofcomp5)
        xx.append(n)
        textlist.append("(%s,%s)" % (n,numofcomp5))
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=xx, y=yy,
                        mode='lines+markers+text',
                        text=textlist,
                        textposition="top left",
                        textfont_size = 15))
    fig.update_xaxes(title_text="num of most common precursors", ticks="inside")
    fig.update_yaxes(title_text="num of synthesizable 5-component materials", ticks="inside")
    fig.update_layout(font=dict(size=30, family='Calibri', color='black'),
                      title={'text': "How many precursors we should buy?"},)
    fig.show()


def make_reactions_list_with_most_common_precursors(n):
    ctr1, all_reactions, targets, energy = count_precursors_frequency_and_get_all_reactions(good_reactions_AA)
    rank_list = ctr1.most_common()[0:n]
    print(rank_list)
    precursors = []
    for pc in rank_list:
        precursors.append(pc[0])
    print("We choose most common %s" % n,"precursors")
    new_targets = []
    new_energy = []
    new_reactions = []
    board = []
    for i in range(n):
        board.append([0 for j in range(n)])
         
    for pc in precursors:
        for reactants in all_reactions:
            reactants.sort()
            if pc in reactants:
                rs_copy = reactants.copy()
                rs_copy.remove(pc)
                if rs_copy[0] in precursors:
                    board[precursors.index(pc)][precursors.index(rs_copy[0])] += 1
                    if reactants not in new_reactions:
                        new_reactions.append(reactants)
                        new_energy.append(energy[all_reactions.index(reactants)])
                        new_targets.append(targets[all_reactions.index(reactants)]) 
    df = pandas.DataFrame()
    df["Target"] = new_targets
    df["Reactants"] = new_reactions
    df["Energy"] = new_energy
    df.to_excel("reactions_list_using_%s_most_common_precursors.xlsx" % n)
    return precursors

if __name__ == "__main__":
    get_pareto_frontier_plot()
    precursors = make_reactions_list_with_most_common_precursors(10)
#     import random
#     random.shuffle(precursors)
#     print(precursors)
 
 
 







