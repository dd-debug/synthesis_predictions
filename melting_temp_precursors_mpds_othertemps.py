
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
    df.to_excel("reactions_list_for_samsung1.xlsx")
    
    
filename = "mp_reaction_for_samsung_comp4_purify"
df = pandas.read_excel(filename +".xlsx",
                       engine='openpyxl',
                       index_col=0)
all_reactants = []
for j,row in df.iterrows():
    
    reactant = get_list_from_reactants_str(row["Reactants"])
    for re in reactant:
        if re not in all_reactants:
            all_reactants.append(re)
print(len(all_reactants))
n=0
filename = "precursors_for_samsung_mp_comp4"
df1 = pandas.read_excel(filename +".xlsx",
                       engine='openpyxl',
                       index_col=0)
reqs = ["decomposition temperature","temperature for eutectoid decomposition","temperature for peritectic formation","temperature for peritectoid formation","temperature for structural transition","pressure derivative of transition temperature"]
for requirement in reqs:
    melting_temps = []
    for re in all_reactants:
        mts = []
        n +=1
        print(re)
        search = {
                    "formulae": re,
                    "props": "phase transitions"
                }
          
        schema = get_data_from_mpds_api(search)
          
        print(schema)
        for e in schema["out"]:
            print(e['sample']['material']['chemical_formula'])
            for m in e['sample']['measurement']:
                mt = "haha"
                if m['property']['name']==requirement:
                    mt = str(m['property']['scalar'])
                    unit = m['property']['units']
                if mt != "haha":
                    if re == e['sample']['material']['chemical_formula']:
                        mts.append(mt)
                    elif re == "".join(e['sample']['material']['chemical_formula'].split("]")[0].split("[")):
                        mts.append(mt)
                    else:
                        mts.append(e['sample']['material']['chemical_formula']+": "+mt)
        print(n)
#         if n % 10 == 0 and n > 155:
#             time.sleep(10)
        if len(mts) == 0:
            melting_temps.append("None")
        else:
            mts.sort()
            melting_temps.append(", ".join(mts))
  
    df1["Reactants"] = all_reactants
    df1[requirement] = melting_temps

df1.to_excel(filename+"1.xlsx")
filename = "mp_reaction_for_samsung_comp4_purify"
filename2 = "precursors_for_samsung_mp_comp41"
combine_precursor_melting_temperature_with_reactions(filename,filename2)









