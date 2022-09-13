'''
Created on Jan 14, 2021

@author: jiadongc
'''

import pandas

def get_list_from_reactants_str(reactant_list_str):
    reactants = reactant_list_str[1:-1].split(", ")
    reactant_list = []
    for name in reactants:
        reactant_list.append(name[1:-1])
    return reactant_list
# good_reactions_AA = pandas.read_excel("mp_reaction_for_samsung_comp4_nonAA_n2.xlsx",
#                        engine='openpyxl',
#                        index_col=0)
# print(len(good_reactions_AA["Target"].tolist()))
# df = pandas.DataFrame()
# targets = []
# all_reactions = []
# energy = []
# invEs = []
# for i, row in good_reactions_AA.iterrows():
# #     if float(row["Energy"]) < -0.1:
#     if row["Target"] not in targets:
#         reactant = get_list_from_reactants_str(row["Reactants"])
#         all_reactions.append(reactant)
#         targets.append(row["Target"])
#         energy.append(row["Energy"])
#         invEs.append(row["invE"])
#     elif energy[targets.index(row["Target"])] != row["Energy"]:
#         reactant = get_list_from_reactants_str(row["Reactants"])
#         all_reactions.append(reactant)
#         targets.append(row["Target"])
#         energy.append(row["Energy"])
#         invEs.append(row["invE"])
# df["Target"] = targets
# df["Reactants"] = all_reactions
# df["Energy"] = energy
# df["invE"]=invEs
# print(len(energy))
# df.to_excel("mp_reaction_for_samsung_comp4_purify.xlsx")

df = pandas.read_excel("ggbro.xlsx",
                       engine='openpyxl',
                       index_col=0)
df1 = pandas.DataFrame()
targets = []
for i, row in df.iterrows():
    if row["Reactants"] in targets:
        continue
    targets.append(row["Reactants"])
    df1 = df1.append(row)
df1.to_excel("gagahaha.xlsx")
print(len(df),len(df1))