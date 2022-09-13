# #This is a test file

# from MPDS_api.ext.mpds_rest import MPDSRester 
# 
# 
# a = MPDSRester('uZAscr1bzEJN4tYgIisXUlkKAHzYc7Ch2TVQ4DEbzTJ5uHHx')
# 
# search = {
#             "classes": "binary, oxide",
#             "props": "temperature for congruent melting"
#         }
# 
# entries=a.mpds_save_query(search)
# 
# for e in entries: 
#     print(e['sample']['material']['chemical_formula'])
#     for m in e['sample']['measurement']:
#         if m['property']['name']=='temperature for congruent melting':
#             print(m['property']['scalar'],m['property']['units']) 
#     
#     print()

# from mpds_client import MPDSDataRetrieval, MPDSDataTypes
# client = MPDSDataRetrieval('uZAscr1bzEJN4tYgIisXUlkKAHzYc7Ch2TVQ4DEbzTJ5uHHx')
# answer = client.get_data(
#     {"elements": "ternary", "props": "temperature for congruent melting"}
# )
# lengths = []
# 
# for item in answer:
#     print(item)

from urllib.parse import urlencode
import httplib2
import json


api_key = "uZAscr1bzEJN4tYgIisXUlkKAHzYc7Ch2TVQ4DEbzTJ5uHHx" # your key
endpoint = "https://api.mpds.io/v0/download/facet"

search = {
            "formulae": "ZnCrO4",
            "props": "temperature for congruent melting"
        }

req = httplib2.Http()
response, content = req.request(
    uri=endpoint + '?' + urlencode({
        'q': json.dumps(search),
        'pagesize': 10,
        'dtype': 1 # see parameters documentation above
    }),
    method='GET',
    headers={'Key': api_key}
)
schema = json.loads(content)
print(schema)
for e in schema["out"]:
    print(e['sample']['material']['chemical_formula'])
    for m in e['sample']['measurement']:
        if m['property']['name']=='temperature for congruent melting':
            print(m['property']['scalar'],m['property']['units']) 
     
    print()