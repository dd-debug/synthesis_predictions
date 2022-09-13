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
filename = "gagahaha"
good_reactions_AA = pandas.read_excel(filename +".xlsx",
                       engine='openpyxl',
                       index_col=0)

def count_precursors_frequency_and_get_all_reactions(good_reactions_AA):
    reactants = []

    all_reactions = []
    for i, row in good_reactions_AA.iterrows():
        if row["Inverse hull energy (eV/atom)"] < -0.010:
            reactant = get_list_from_reactants_str(row["Reactants"])
            reactants += reactant
            all_reactions.append(reactant)
    ctr1 = collections.Counter(reactants)
    print(ctr1)
    print(len(ctr1),"precursors")
    print(len(all_reactions),"reactions with invE < -0.015 eV/atom")
    return ctr1, all_reactions

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
    for n in range(1,164):
        numofcomp5 = num_comp5_with_most_common_precursors(n)
        yy.append(numofcomp5)
        xx.append(n)
        textlist.append("(%s,%s)" % (n,numofcomp5))
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=xx, y=yy,
                        mode='lines+markers+text',))
#                         text=textlist,
#                         textposition="top left",
#                         textfont_size = 15))
    fig.update_xaxes(title_text="num of most common precursors", ticks="inside")
    fig.update_yaxes(title_text="num of synthesizable 5-component materials", ticks="inside")
    fig.update_layout(font=dict(size=30, family='Calibri', color='black'),
                      title={'text': "How many precursors we should buy?"},)
    fig.show()
    return fig

def make_reactions_list_with_most_common_precursors(n):
    ctr1 = count_precursors_frequency_and_get_all_reactions(good_reactions_AA)[0]
    rank_list = ctr1.most_common()[0:n]
    print(rank_list)
    precursors = []
    for pc in rank_list:
        precursors.append(pc[0])
    print("We choose most common %s" % n,"precursors")
    df = pandas.DataFrame()
    for i, row in good_reactions_AA.iterrows():
        reactants = get_list_from_reactants_str(row["Reactants"])
        r1 = reactants[0]
        r2 = reactants[1]
        if r1 in precursors:
            if r2 in precursors:
                if row["Inverse hull energy (eV/atom)"] < -0.010:

                    df = df.append(row,ignore_index=True)
    df.to_excel("reactions_list_using_%s_most_common_precursors.xlsx" % n)
    return precursors

def get_precursors_bar_plot():
    ctr1, all_reactions = count_precursors_frequency_and_get_all_reactions(good_reactions_AA)
    rank_list = ctr1.most_common()
    print(len(rank_list))
    els = []
    fres = []
    for el, fre in rank_list:
        els.append(el)
        fres.append(fre)
    fig = go.Figure([go.Bar(x=els, y=fres)])
    fig.update_traces(texttemplate='%{y}', textposition='outside')
    fig.update_xaxes(title_text="Precursors", ticks="inside")
    fig.update_yaxes(title_text="Frequency", ticks="inside")
    fig.update_layout(font=dict(size=30, family='Calibri', color='black'),
                      template="ggplot2",
                      title={'text': "Ocurrance of each precursor"},
                      hoverlabel = dict(font_size = 30, font_family='Calibri'))
    fig.show()
    return fig 

def how_many_materials_precursors_can_form():
    df1 = pandas.read_excel("gagahaha" +".xlsx",
                       engine='openpyxl',
                       index_col=0)
    pr1 = pandas.read_excel("ggbro" +".xlsx",
                       engine='openpyxl',
                       index_col=0)
    D = {}
    for j,rowpr in pr1.iterrows():
        D[rowpr["Reactants"]] = []
    for i, row in df1.iterrows():
        reactants = get_list_from_reactants_str(row["Reactants"])
        D[reactants[0]].append(row["Target"])
        D[reactants[1]].append(row["Target"])
    howmany = []
    whichm = []
    for j,rowpr in pr1.iterrows():
        whichm.append(D[rowpr["Reactants"]])
        howmany.append(len(D[rowpr["Reactants"]]))
    pr1["how_many_materials_can_form"] = howmany
    pr1["which_materials_can_form"] = whichm
    pr1.to_excel("tryyitry.xlsx")
    
if __name__ == "__main__":
#     fig1 = get_pareto_frontier_plot()
#     fig1.write_html("pareto_frontier_for_samsung.html")
# #     fig1.write_image("pareto_frontier_for_samsung.png")
    fig2 = get_precursors_bar_plot()
# #     fig2.write_image("frequency_of_precursors.png")
# #     precursors = make_reactions_list_with_most_common_precursors(43)
    how_many_materials_precursors_can_form()
 
 
 







