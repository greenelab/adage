'''
This is used to get operon information from the online database DOOR.
'''

import json
import urllib2
import pickle 
import sys

url_h = urllib2.urlopen('http://csbl.bmb.uga.edu/DOOR/search_ajax.php?keyword=pao1&mode=DataTable&sEcho=4&iColumns=6&sColumns=&iDisplayStart=0&iDisplayLength=6000&mDataProp_0=0&mDataProp_1=1&mDataProp_2=2&mDataProp_3=3&mDataProp_4=4&mDataProp_5=5&_=1384798560566', timeout=5)
obj = json.loads(url_h.read())
operons = []
operon_file = open(sys.argv[1],'w')
constraint = int(sys.argv[2])

for item in obj['aaData']:
    operon = item[3].split('; ')
    if len(operon) < constraint: #Only operons that have more genes than constraint will be kept.
        continue
    operon = [str(x) for x in operon]
    operons.append(operon)
    operon_file.write('\t'.join([x for x in operon])+'\n')

