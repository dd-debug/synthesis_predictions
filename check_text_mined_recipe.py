'''
Created on May 1, 2021

@author: jiadongc
'''
import json
import pandas
filename = "gagahaha"
df = pandas.read_excel(filename +".xlsx",
                       engine='openpyxl',
                       index_col=0)
targes = list(df["Target"])
targets = []
for i in targes:
    if i not in targets:
        targets.append(i)
print(len(targets))
  
with open("text_mined_recipes.json","r") as f:
    txt_recipes = json.load(f)
  
recipes_overlap = []
for i in txt_recipes:
#     print(i["targets_string"][0])
    if i["targets_string"][0] in targets:
        recipes_overlap.append(i)
  
print(len(recipes_overlap))
ol = []
for i in recipes_overlap:
    if i["targets_string"][0] not in ol:
        ol.append(i["targets_string"][0])
print(len(ol))
# for i in ol:
#     print(i)
print()
for i in targets:
    if i in ol:
        print(i)
#          
with open("overlap_recipes","w") as f:
    json.dump(recipes_overlap,f)


with open("overlap_recipes","r") as f:
    overlap_recipes = json.load(f)
for i in overlap_recipes:
    if i["targets_string"][0] == "Li2VSiO5":
        try:
            print(i)
        except:
            print("emm")
        print(i["reaction_string"])
        print()

'''precursors that have recipes in text-mined database'''
# st = "LiPO3,Li2O,LiAsO3,LiBO2,Li3BO3,GeO2,LiSO3F,Li2SiO3,Li4GeO4,Li4SiO4,LiOsO3,SiO2,SrO,Li2Si2O5,LiAg3O2,LiCrO2,LiF,LiGaO2,LiHO,LiP(HO2)2,LiVO2,MnO,ZnO"
# prs = st.split(",")
# print(len(prs))
# 
# with open("text_mined_recipes.json","r") as f:
#     txt_recipes = json.load(f)
#   
# recipes_overlap = []
# for i in txt_recipes:
# #     print(i["targets_string"][0])
#     if i["targets_string"][0] in prs:
#         recipes_overlap.append(i)
# print(len(recipes_overlap))
# ol = []
# for i in recipes_overlap:
#     if i["targets_string"][0] not in ol:
#         ol.append(i["targets_string"][0])
# print(len(ol))
# for i in ol:
#     print(i)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        