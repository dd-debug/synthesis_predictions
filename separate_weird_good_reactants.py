'''
Created on Feb 13, 2021

@author: jiadongc
'''
import pandas
from pymatgen.core.composition import Composition
def get_list_from_reactants_str(reactant_list_str):
    reactants = reactant_list_str[1:-1].split(", ")
    reactant_list = []
    for name in reactants:
        reactant_list.append(name[1:-1])
    return reactant_list
 
filename = "good_reaction_for_samsung_n"
df = pandas.read_excel(filename +".xlsx",
                       engine='openpyxl',
                       index_col=0)

df1 = pandas.DataFrame()
df2 = pandas.DataFrame()

weird_reactants = ["LiO8","Mn","Zn","SO2","SO3"]
for i,row in df.iterrows():
    if row["Reactants1"] in weird_reactants or row["Reactants2"] in weird_reactants:
        df1 = df1.append(row,ignore_index=True)
    else:
        df2 = df2.append(row,ignore_index=True)
df2.to_excel("reactions_with_good_Reactants_for_samsung.xlsx")
df1.to_excel("reactions_with_weird_Reactants_for_samsung.xlsx")