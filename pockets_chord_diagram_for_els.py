'''
Created on Jan 7, 2021

@author: jiadongc
'''
import pandas
import collections
from itertools import permutations
from pymatgen.core.composition import Composition
from chord import Chord
import plotly.graph_objects as go

Chord.user = "jdchen1997@outlook.com"
Chord.key = "CP-c4615f36-9908-43f4-bdbe-97ff1b9f3fbd"

# filename = "good_reactions_AA_purify"
# good_reactions_AA = pandas.read_csv(filename +".csv") #dataframe
filename = "reactions_target_is_deepest_in_IRPD"
good_reactions_AA = pandas.read_excel(filename +".xlsx",
                       engine='openpyxl',
                       index_col=0)
def count_els_frequency(good_reactions_AA):
    all_els = []
    for i, row in good_reactions_AA.iterrows():
        if row["Energy"] < -0.1: # only consider reaction with formE < -0.1 eV/atom
            els = [str(el) for el in Composition(row["Target"]).elements]
            all_els += els
    ctr = collections.Counter(all_els)
    print(ctr)
    return ctr


def display_board(board):
    for ii in board:
        print(ii)

def get_els_chord_diagram(good_reactions_AA, filename):
    ctr = count_els_frequency(good_reactions_AA)
    all_els = sorted(list(ctr.keys()))
#     all_els.remove("O")
#     all_els.remove("S")
    print(len(all_els), "els out of 32 els")
    # build a 29*29 matrix, with value is 0
    board = []
    for i in range(len(all_els)):
        board.append([0 for j in range(len(all_els))])
    
    for i,row in good_reactions_AA.iterrows():
        if row["Energy"] < -0.1: # only consider reaction with reactE < -0.1 eV/atom
            els = [str(el) for el in Composition(row["Target"]).elements]
            permus = list(permutations(els,2))
            for permu in permus:
                board[all_els.index(permu[0])][all_els.index(permu[1])] += 1
    
    display_board(board)
    
    # Chord(board, all_els).show()
    title = "Chord Diagram for " + filename + " with reactE < -0.1 eV/atom"
    cd = Chord(board, all_els, opacity = 0.6, font_size="12px",
            font_size_large="20px",title=title)
    cd.to_html("chord_diagram_for_" + filename + ".html")
    return cd

def get_bar_plot(good_reactions_AA):
    ctr = count_els_frequency(good_reactions_AA)
    rank_list = ctr.most_common()
    els = []
    fres = []
    for el, fre in rank_list:
        els.append(el)
        fres.append(fre)
    fig = go.Figure([go.Bar(x=els, y=fres)])
    fig.update_traces(texttemplate='%{y}', textposition='outside')
    fig.update_xaxes(title_text="Elements", ticks="inside")
    fig.update_yaxes(title_text="Frequency", ticks="inside")
    fig.update_layout(font=dict(size=40, family='Calibri', color='black'),
                      template="ggplot2",
                      title={'text': "Ocurrance of each element"},
                      hoverlabel = dict(font_size = 30, font_family='Calibri'))
    fig.show()
    return fig
    
if __name__ == "__main__":
    get_bar_plot(good_reactions_AA)
    get_els_chord_diagram(good_reactions_AA, filename)
    
    
    
    
    
    
    