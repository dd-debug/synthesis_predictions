'''
Created on Jan 6, 2021

@author: jiadongc
'''
import pandas
import collections
from itertools import permutations
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
filename = "good_reaction_for_samsung_n"
# filename = "tryyitry"
good_reactions_AA = pandas.read_excel(filename +".xlsx",
                       engine='openpyxl',
                       index_col=0)
def count_precursors_frequency_and_get_all_reactions(good_reactions_AA):
    reactants = []
    all_reactions = []
    for i, row in good_reactions_AA.iterrows():
        if float(row["Inverse hull energy (eV/atom)"]) < -0.015:
            reactant = get_list_from_reactants_str(row["Reactants"])
            reactants += reactant
            all_reactions.append(reactant)
    
    print(len(all_reactions), "out of %s reactions with reactE < -0.1" % (len(good_reactions_AA["Reactants"])))
    ctr1 = collections.Counter(reactants)
    print(len(ctr1))
    return ctr1, all_reactions

def get_precursors_chord_diagram():
    ctr1, all_reactions = count_precursors_frequency_and_get_all_reactions(good_reactions_AA)
    precursors = list(ctr1.keys())
    print(len(precursors), "precursors")
    
    board = []
    for i in range(len(precursors)):
        board.append([0 for j in range(len(precursors))])
        
    for reactants in all_reactions:
        permus = list(permutations(reactants,2))
        for permu in permus:
            board[precursors.index(permu[0])][precursors.index(permu[1])] += 1
    
    '''check if the transpose is equal of the board'''
    # rez = [[board[j][i] for j in range(len(board))] for i in range(len(board[0]))]
    # if rez == board:
    #     print("haha") 
    title = "Chord Diagram for precursors with reactE < -0.1 eV/atom"
    cd = Chord(board, precursors, opacity = 0.6,
            font_size_large="10px",title=title)
    cd.to_html("chord_diagram_for_precursors.html")
    return cd



def get_most_common_precursors_chord_diagram(n):
    ctr1, all_reactions = count_precursors_frequency_and_get_all_reactions(good_reactions_AA)
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
    
    
    '''check if the transpose is equal of the board'''
#     rez = [[board[j][i] for j in range(len(board))] for i in range(len(board[0]))]
#     if rez == board:
#         print("haha") 

    title = "Chord Diagram for most common %s precursors with reactE < -0.1 eV/atom" % n
    cd = Chord(board, precursors, opacity = 0.6,
            font_size_large="20px",title=title)
    cd.to_html("chord_diagram_for_most_common_%s_precursors.html" % n)
    return cd

def get_precursors_bar_plot():
    ctr1, all_reactions = count_precursors_frequency_and_get_all_reactions(good_reactions_AA)
    rank_list = ctr1.most_common()
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


if __name__ == "__main__":
#     get_most_common_precursors_chord_diagram(13)
#     get_precursors_chord_diagram()
#     get_precursors_bar_plot()
    ctr = count_precursors_frequency_and_get_all_reactions(good_reactions_AA)[0]
    ctrlist = ctr.most_common()
    print(ctrlist)
    for tu in ctrlist:
        print(tu)
    precursors = [tu1[0] for tu1 in ctrlist]
#     allprs = ['Li2ZrO3', 'LiScO2', 'LiVO2', 'LiMoO2', 'Li2TiO3', 'LiNbO2', 'LiCrO2', 'LiYO2', 'LiFeO2', 'ZnMoO4', 'LiMnO2', 'LiPO3', 'Li(MoO2)2', 'MnMoO4', 'YMnO3', 'MnSO4', 'ZnSO4', 'Li2RhO3', 'CdMoO4', 'LiCoO2', 'CuSO4', 'CdSO4', 'LiTiO2', 'VSO4', 'LiTiVO4', 'NiSO4', 'NiMoO4', 'Li2FeO2', 'Li2MoO4', 'VSO5', 'LiVCrO4', 'ZnCrO4', 'AgSO4', 'YCrO4', 'CoSO4', 'Li2FeO3', 'FeSO4', 'LiNbO3', 'CrMoO4', 'PdSO4', 'CoPO4', 'LiRhO2', 'Li2CrO4', 'YCoO3', 'LiTiCrO4', 'LiTiPO5', 'LiVO3', 'TiNiO3', 'LiVPO4', 'NbMoO4', 'LiTi2O4', 'Li2NiO3', 'Li2MnO3', 'YCrO3', 'LiVCoO4', 'CdO', 'LiO8', 'LiMnCrO4', 'LiNiO2', 'TiVO3', 'LiMn2O4', 'TiFeO3', 'TiSO5', 'Li2PdO3', 'MnPO4', 'LiCuO2', 'LiMnPO4', 'YVO3', 'Li3NbO4', 'Li3VO4', 'LiMnNbO4', 'LiMnVO4', 'CrCdO4', 'LiFePO4', 'CrO2', 'CrNiO4', 'ZnO', 'LiVPO5', 'TiCdO3', 'Ti2FeO5', 'CuPO4', 'Li(NiO2)2', 'Li3CrO4', 'YNiO3', 'Li2FeP2O7', 'Li2VOF4', 'YMn2O5', 'LiV2O5', 'AgO', 'NiO', 'MnCoO3', 'MnNiO3', 'PRhO4', 'MnO', 'FeO', 'Li3Mn3(PO4)4', 'LiAgO2', 'CrPO4', 'Li3PO4', 'FePO4', 'YRhO3', 'AgPO3', 'MoPO5', 'Y2MoO6', 'Li2PdO2', 'CoO', 'Zn(FeO2)2', 'Li2MnF6', 'TiO', 'ScVO4', 'LiZnPO4', 'LiCoPO4', 'LiVZnO4', 'ScPO4', 'VAgO3', 'Li2CoO3', 'Li(CoO2)2', 'VO2', 'YVO4', 'Y(PdO2)2', 'Ti2CoO5', 'CoO2', 'SO2', 'Cu2O', 'YPO4', 'NbO2', 'LiScP2O7', 'Mn(FeO2)2', 'Zn(PO3)2', 'Mn2PO5', 'Y2O3', 'Cd(FeO2)2', 'Y(CuO2)2', 'NbO', 'TiMnO3', 'LiNiPO4', 'Y2TiO5', 'Cu2O3', 'Ag2O', 'LiNi2P3O10', 'MoO3', 'Li2CrF6', 'Zn3(PO4)2', 'Li3Fe3(PO4)4', 'Cd(PO3)2', 'Mn(PO3)2', 'Mn2CoO4', 'Li2VF6', 'Ni3(PO4)2', 'Sc2O3', 'CrCoO4', 'LiP', 'NbZnRu2', 'Zr2CoP', 'MoO2', 'YAgS2', 'LiMnF4', 'CuF2', 'CrF2', 'ScAgS2', 'LiCuO', 'Mn', 'MnO2', 'MnCuO2', 'Cr2O3', 'Zn(CoO2)2', 'Mn2ZnO4', 'Li3Mn4(FeO6)2', 'AgF2', 'AgRuO4', 'CrAgS2', 'Cd(CoO2)2', 'VOF', 'Li3CuF6', 'FeCuO2', 'CoNiO3', 'YCuS2', 'VPO4', 'Co(NiO2)2', 'Nb2ZnO6', 'CuO', 'PdO2', 'Li2Mn3NiO8', 'Fe2O3', 'TiO2', 'Sc2TiO5', 'Nb2FeO6', 'Li2Ti3NiO8', 'Li3CuO3', 'ScAgO2', 'CrNiPO5', 'VCdO3', 'ZnCr2O4', 'YMnFeO5', 'Cu3(PO4)2', 'Nb2CoO6', 'Ti2O', 'Ag3(PO4)2', 'MnNbRu2', 'CuAgO2', 'CuRhO2', 'Zn', 'Li4P2O7', 'V2Cd2O7', 'Li2ZrF6', 'TiAgPO5', 'FeF2', 'MnAgO2', 'MnCr2O4']
    print(precursors)
    print(len(precursors))
    dd = ["Reactants","Melting_temperture","decomposition temperature","temperature for structural transition","temperature for peritectic formation","temperature for eutectoid decomposition"]
    
    D = {e:[] for e in dd}
    dfnew = pandas.DataFrame(D)
    df1 = pandas.read_excel("precursors_for_samsung" +".xlsx",
                           engine='openpyxl',
                           index_col=0)
    for prs in ctrlist:
        for i,row in df1.iterrows():
            if row["Reactants"] == prs[0]:
                dfnew = dfnew.append(row)
    for i,row in df1.iterrows():
        if row["Reactants"] not in precursors:
            dfnew = dfnew.append(row)
    dfnew.to_excel("precursors_for_samsung_rank_based_on_frequency.xlsx")
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    