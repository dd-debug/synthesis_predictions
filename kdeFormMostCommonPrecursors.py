import json
import plotly.figure_factory as ff
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import glob
import pandas
from pymatgen.analysis.phase_diagram import PhaseDiagram
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pockets_chord_diagram_for_precursorsAA import count_precursors_frequency_and_get_all_reactions
from pymatgen.core.composition import Composition
import plotly
def get_formE_kde_volcano_background():
    with open('formEdis','r') as f:
        D = json.load(f)
    Binary = D['comp2']
    Ternary = D['comp3']
    FourconvE = D['comp4']
    FiveconvE = D['comp5']
    SixconvE = D['comp6']
        
    fig2 = ff.create_distplot([Binary], ['allBinary'], bin_size=.2,show_hist=False,show_rug=False)
    trace2 = fig2['data'][0]
    trace2['showlegend']=False
    trace2['opacity'] = 0.2
    trace2['line'] = dict(color ="rgb(13, 109, 247)" ,width = 3)
    print(fig2['data'][0])
    
    fig3 = ff.create_distplot([Ternary], ['allTernary'], bin_size=.2,show_hist=False,show_rug=False)
    trace3 = fig3['data'][0]
    trace3['showlegend']=False
    trace3['opacity'] = 0.2
    trace3['line'] = dict(color ="rgb(247, 36, 13)",width = 3)
    print(fig3['data'][0])
    
    fig4 = ff.create_distplot([FourconvE], ['allQuanternary'], bin_size=.2,show_hist=False,show_rug=False)
    trace4 = fig4['data'][0]
    trace4['showlegend']=False
    trace4['opacity'] = 0.2
    trace4['line'] = dict(color ="rgb(59, 202, 21)",width = 3)
    print(fig4['data'][0])
    
    
    fig = make_subplots(rows=4, cols=1,shared_xaxes=True, vertical_spacing=0.02)#subplot_titles=("Plot 1", "Plot 2", "Plot 3", "Plot 4",'Plot 5')
    
    fig.append_trace(trace2, 1, 1)
    fig.append_trace(trace3, 2, 1)
    fig.append_trace(trace4, 3, 1)
    
    
    with open("allMatFormEFor5compMat", 'r') as f:
        Dict = json.load(f)
    for e in Dict:
        print(e,len(Dict[e]),Dict[e])
    
    Binary = list(Dict['comp2'].values())
    Ternary = list(Dict['comp3'].values())
    FourinverseE = list(Dict['comp4'].values())
    FiveinverseE = list(Dict['comp5'].values())
    
    
    fig2 = ff.create_distplot([Binary], ['Binary'], bin_size=.2,show_hist=False,show_rug=False)
    trace2 = fig2['data'][0]
    # trace2['opacity'] = 0.2
    # trace2['showlegend']=False
    trace2['line'] = dict(color ="rgb(13, 109, 247)" ,width = 3)
    print(fig2['data'][0])
    
    fig3 = ff.create_distplot([Ternary], ['Ternary'], bin_size=.2,show_hist=False,show_rug=False)
    trace3 = fig3['data'][0]
    # trace3['opacity'] = 0.2
    trace3['line'] = dict(color ="rgb(247, 36, 13)",width = 3)
    print(fig3['data'][0])
    
    fig4 = ff.create_distplot([FourinverseE], ['Quaternary'], bin_size=.2,show_hist=False,show_rug=False)
    trace4 = fig4['data'][0]
    # trace4['opacity'] = 0.2
    trace4['line'] = dict(color ="rgb(59, 202, 21)",width = 3)
    print(fig4['data'][0])
    
    fig5 = ff.create_distplot([FiveconvE], ['Pentanary'], bin_size=.2,show_hist=False,show_rug=False)
    trace5 = fig5['data'][0]
    trace5['opacity'] = 1
    trace5['line'] = dict(color ="rgb(172, 24, 247)",width = 3)
    print(fig5['data'][0])
    
    
    
    fig.append_trace(trace2, 1, 1)
    fig.append_trace(trace3, 2, 1)
    fig.append_trace(trace4, 3, 1)
    fig.append_trace(trace5, 4, 1)
    return fig


def get_most_common_precursors_vertical_lines(rank_list, fig):
    D1 = {"Binaryx":[], "Binaryy":[],"Binaryname":[],
          "Ternaryx":[], "Ternaryy":[],"Ternaryname":[],
          "Quanternaryx":[], "Quanternaryy":[],"Quanternaryname":[]}
    D2 = {"Binaryx":[], "Binaryy":[],"Binaryname":[],
          "Ternaryx":[], "Ternaryy":[],"Ternaryname":[],
          "Quanternaryx":[], "Quanternaryy":[],"Quanternaryname":[]}
    for tu in rank_list:
        els = [str(el) for el in Composition(tu[0]).elements]
        entries = getOrigStableEntriesList(els)
        pd = PhaseDiagram(entries)
        for e in entries:
            if e.name == tu[0]:
                if len(e.composition.elements) == 2:
                    formE = pd.get_form_energy_per_atom(e)
                    D1["Binaryx"] += [formE,formE,None]
                    D1["Binaryy"] += [0,0.55,None]
                    D1["Binaryname"] += [str(tu), str(tu), None]
                    D2["Binaryx"] += [formE]
                    D2["Binaryy"] += [0.55]
                    D2["Binaryname"] += [str(tu)]
                if len(e.composition.elements) == 3:
                    formE = pd.get_form_energy_per_atom(e)
                    D1["Ternaryx"] += [formE,formE,None]
                    D1["Ternaryy"] += [0,0.42,None]
                    D1["Ternaryname"] += [str(tu), str(tu), None]
                    D2["Ternaryx"] += [formE]
                    D2["Ternaryy"] += [0.42]
                    D2["Ternaryname"] += [str(tu)]
                if len(e.composition.elements) == 4:
                    formE = pd.get_form_energy_per_atom(e)
                    D1["Quanternaryx"] += [formE,formE,None]
                    D1["Quanternaryy"] += [0,0.32,None]
                    D1["Quanternaryname"] += [str(tu), str(tu), None]
                    D2["Quanternaryx"] += [formE]
                    D2["Quanternaryy"] += [0.32]
                    D2["Quanternaryname"] += [str(tu)]
    trace22 = go.Scatter(x=D1["Binaryx"], y=D1["Binaryy"], mode="lines", line=dict(color='rgb(30,140,220)',width=2),
                         text = D1["Binaryname"], hovertemplate=
            "<b>%{text}</b><br>" +
            "formE: %{x}<br>",opacity=1,name="Binary")
    trace33 = go.Scatter(x=D1["Ternaryx"], y=D1["Ternaryy"], mode="lines", line=dict(color='rgb(228,65,38)',width=2),
                         text = D1["Ternaryname"],hovertemplate=
            "<b>%{text}</b><br>" +
            "formE: %{x}<br>",opacity=1,name="Ternary")
    trace44 = go.Scatter(x=D1["Quanternaryx"], y=D1["Quanternaryy"], mode="lines", line=dict(color='rgb(45,220,80)',width=2),
                         text = D1["Quanternaryname"], hovertemplate=
            "<b>%{text}</b><br>" +
            "formE: %{x}<br>",opacity=1,name="Quanternary")
    fig.append_trace(trace22, 1, 1)
    fig.append_trace(trace33, 2, 1)
    fig.append_trace(trace44, 3, 1)
    fsize = 10
    for i in range(len(D2["Binaryx"])):
        fig.add_annotation(go.layout.Annotation(
            x=D2["Binaryx"][i],
            y=D2["Binaryy"][i],
            text=D2["Binaryname"][i], align="left",
            showarrow=False,
            yanchor='bottom',
            textangle=-90,font_size = fsize),row = 1,col = 1)
    for i in range(len(D2["Ternaryx"])):
        fig.add_annotation(go.layout.Annotation(
            x = D2["Ternaryx"][i],
            y = D2["Ternaryy"][i],
            xref="x",
            yref="y",
            text = D2["Ternaryname"][i], align="left",
            showarrow=False,
            yanchor='bottom',
            textangle=-90,font_size = fsize),row = 2,col = 1)
    for i in range(len(D2["Quanternaryx"])):
        fig.add_annotation(go.layout.Annotation(
            x=D2["Quanternaryx"][i],
            y=D2["Quanternaryy"][i],
            text=D2["Quanternaryname"][i], align="left",
            showarrow=False,
            yanchor='bottom',
            textangle=-90,font_size = fsize),row = 3,col = 1)
    # fig.update_xaxes(title_text="Formation energy (eV/atom)",row = 4, col = 1)
    # fig.update_yaxes(title_text="Kernel Density",row = 2.5, col = 1)
    fig.update_xaxes(range = [-4.1,0])
    fig.update_layout(showlegend = False,legend_orientation="v",font=dict(size=16, family='Arial', color='black'),
                      paper_bgcolor='rgb(255,255,255)',
                      template="plotly_white")
    fig.update_xaxes(title_font=dict(size=20, family='Arial', color='black'),showgrid=False,zeroline=False,showline=True, linewidth=2, linecolor='black',mirror=True)
    fig.update_yaxes(title_font=dict(size=20, family='Arial', color='black'),showgrid=False, zeroline=False,showline=True, linewidth=2, linecolor='black',mirror=True)
    
    fig.update_xaxes(ticks="inside",row = 5,col = 1)
    fig.update_yaxes(ticks="inside")
    return fig
    # fig.show()
    # plotly.io.orca.config.executable = r'C:\Users\jiadongc\anaconda3\pkgs\plotly-orca-1.3.1-1\orca_app\orca.exe'
    # plotly.io.orca.config.save()
#     fig.write_image("kdeFormMostCommonPrecursors.png",width=900,height=900,scale=3)



if __name__ == "__main__":
    fig = get_formE_kde_volcano_background()
    good_reactions_AA = pandas.read_csv("good_reactions_AA.csv")
    ctr, all_reactions = count_precursors_frequency_and_get_all_reactions(good_reactions_AA)
    rank_list = ctr.most_common()[0:11]
    print(rank_list)
    fig = get_most_common_precursors_vertical_lines(rank_list, fig)
    fig.show()
#     fig.write_image("kdeFormMostCommonPrecursors.png",width=900,height=900,scale=3)
    
    
    
    
    
    
    
    
    
    
    
    

