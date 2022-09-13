
from urllib.parse import urlencode
import httplib2
import json
import pandas
import time
import os 

API_KEY = "uZAscr1bzEJN4tYgIisXUlkKAHzYc7Ch2TVQ4DEbzTJ5uHHx"
def get_list_from_reactants_str(reactant_list_str):
    reactants = reactant_list_str[1:-1].split(", ")
    reactant_list = []
    for name in reactants:
        reactant_list.append(name[1:-1])
    return reactant_list
def construct_unique_key(params):
    connector = "_"
    params_list = []
    for e in params:
        params_list.append(f'{e}_{"".join(params[e].split())}')
    params_list.sort()
    key_str = connector.join(params_list)
#     print(key_str)
    return key_str

directory = os.path.join("C:/Users/jiadongc/dataset/melting_temp_for_precursors_mpds")

def get_data_from_mpds_api(search): # your key
    fn = construct_unique_key(search)
    cache = os.path.join(directory, fn)
    if os.path.exists(cache):
        print('loading from cache.')
        with open(cache,'r') as f:
            result = json.load(f)
        return result
    else:
        print("reading from db.")
        endpoint = "https://api.mpds.io/v0/download/facet"
        req = httplib2.Http()
        response, content = req.request(
            uri=endpoint + '?' + urlencode({
                'q': json.dumps(search),
                'pagesize': 10,
                'dtype': 7 # see parameters documentation above
            }),
            method='GET',
            headers={'Key': API_KEY}
        )
        schema = json.loads(content)
        with open(cache,"w") as f:
            json.dump(schema,f)
        return schema
    
def combine_precursor_melting_temperature_with_reactions(filename,filename1):
#     filename = "good_reaction_for_samsung_n"
    df = pandas.read_excel(filename +".xlsx",
                           engine='openpyxl',
                           index_col=0)
#     filename1 = "precursors_for_samsung"
    df1 = pandas.read_excel(filename1 +".xlsx",
                           engine='openpyxl',
                           index_col=0)
    Reactants1 = []
    Reactants2 = []
    melting_temps1 = []
    melting_temps2 = []
    decomptemps1=[]
    decomptemps2=[]
    for i,row in df.iterrows():
        reactants = get_list_from_reactants_str(row["Reactants"])
        Reactants1.append(reactants[0])
        Reactants2.append(reactants[1])
        for j,row2 in df1.iterrows():
            if row2["Reactants"] == reactants[0]:
                melting_temps1.append(row2["Melting_temperture"])
                decomptemps1.append(row2["decomposition temperature"])
            if row2["Reactants"] == reactants[1]:
                melting_temps2.append(row2["Melting_temperture"])
                decomptemps2.append(row2["decomposition temperature"])
    df["Reactants1"] = Reactants1
    df["Melting_temps1_mpds_K"] = melting_temps1
    df["decomposition_temps1_mpds_K"] = decomptemps1
    df["Reactants2"] = Reactants2
    df["Melting_temps2_mpds_K"] = melting_temps2
    df["decomposition_temps2_mpds_K"] = decomptemps2
    df.to_excel("reactions_list_for_samsung.xlsx")
    
def isNone(str1):
    if str1 == "None":
        result = 10000
    else:
        if ", " in str1:
            strlist = str1.split(", ")
        else:
            strlist = [str1]
        templist = []
        for i in strlist:
            if ": " in i:
                i = i.split(": ")[-1]
            if "+/-" in i:
                temp1 = float(i.split("+/-")[0])-280
                templist.append(temp1)
            else:
                templist.append(float(i))
        result = min(templist)
    return result 
                
'''min decomposition or melting temperature'''
filename = "mp_reaction_for_samsung_comp4_purify"
df = pandas.read_excel(filename +".xlsx",
                       engine='openpyxl',
                       index_col=0)
mintemp = []
for j,row in df.iterrows():
    # find the min of 4 temps
    a = []
    a.append(isNone(row["Melting_temps1_mpds_K"]))
    a.append(isNone(row["decomposition_temps1_mpds_K"]))
    a.append(isNone(row["Melting_temps2_mpds_K"]))
    a.append(isNone(row["decomposition_temps2_mpds_K"]))
    minn = min(a)
    if minn == 10000:
        mintemp.append("None")
    else:
        mintemp.append(minn)
df["mintemp"] = mintemp
df.to_excel("reactions_list_for_samsung.xlsx")



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    