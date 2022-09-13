'''
Created on Jan 6, 2021
 
@author: jiadongc
'''
import pandas
import collections
import numpy as np
from chord import Chord
import plotly.graph_objects as go
 
 
def get_list_from_reactants_str(reactant_list_str):
    reactants = reactant_list_str[1:-1].split(", ")
    reactant_list = []
    for name in reactants:
        reactant_list.append(name[1:-1])
    return reactant_list
  
def make_hyperlink(value,urlfoldername="mpds_mp_samsung"):
    url = urlfoldername+"/{}"
    return '=HYPERLINK("%s", "%s")' % (url.format(value+".html"), value+"_html")

def add_hyperlinks(filename, urlfoldername = "mpds_mp_samsung"): 
    df = pandas.read_excel(filename+".xlsx",
                           engine='openpyxl',
                           index_col=0)
    df['IRPD_hyperlink'] = df['Target'].apply(lambda x: make_hyperlink(x,urlfoldername))
    df.to_excel(filename+".xlsx")
filename = "mp_reaction_for_samsung_comp4_purify"
add_hyperlinks(filename,urlfoldername = "mpds_mp_samsung")





# # df = pandas.read_excel("reactions_target_is_deepest_in_IRPD.xlsx",
# #                        engine='openpyxl',
# #                        index_col=0)
# # print(df)

#   
#   
# df1 = pandas.read_excel("comp5_Reactants_and_Reactions.xlsx",
#                        engine='openpyxl',
#                        index_col=0)
# reactants_melting_temperature = []
# reactants_decomposition_temperature = []
# for a, rowdf in df.iterrows():
#     reactants = get_list_from_reactants_str(rowdf["Reactants"])
#     melt_str = ""
#     decomp_str = ""
#     for re in reactants:
#         for i, row in df1.iterrows():
#             if row["Reactant"] == re:
# #                 print(row)
# #                 print(type(row["Melting Temperature [C]"]))
#                 melt_str += f"{re}: " f"{row['Melting Temperature [C]']}; "
#                 decomp_str += f"{re}: " f"{row['Decomposition Temperature [C]']}; "
#     print(melt_str)
#     print(decomp_str)
#     print()
#     reactants_melting_temperature.append(melt_str)
#      
#     if decomp_str.count("Unknown") < 2:
#         reactants_decomposition_temperature.append(decomp_str)
#     else:
#         reactants_decomposition_temperature.append("Unknown")
# df['Phasediagram_hyperlink'] = df['Target'].apply(lambda x: make_hyperlink(x))
# df["Reactants_melting_temperature [C]"] = reactants_melting_temperature
# df["Reactants_decomposition_temperature [C]"] = reactants_decomposition_temperature
#   
# df.to_excel("34reactions_list_with_32precursors1.xlsx")
# 
# # df = pandas.read_excel("34reactions_list_with_32precursors_peter1.xlsx",
# #                        engine='openpyxl',
# #                        index_col=0)
# # df1 = pandas.read_excel("34reactions_list_with_32precursors.xlsx",
# #                        engine='openpyxl',
# #                        index_col=0)
# # temps = []
# # print(df1)
# # for a, rowdf in df1.iterrows():
# #     for b, row in df.iterrows():
# #         if row["Target"] == rowdf["Target"]:
# #             temps.append(row["Reactants_melting/decomposition_temperature [C]"])
# # df1["Reactants_melting/decomposition_temperature [C]"]=temps
# # df1.to_excel("34reactions_list_with_32precursors1.xlsx")
 
 
 
 

