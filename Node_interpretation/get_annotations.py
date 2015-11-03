import requests
import sys

TRIBE_URL = 'http://tribe.greenelab.com'
limit = 2000

def get_terms(offset=None, target = None):
    parameters = {'show_tip': 'true', 'organism': 9, 'xrdb': 'Symbol', 'limit': limit, 'offset': offset} 
    r = requests.get(TRIBE_URL + '/api/v1/geneset/?title__startswith='+target, params=parameters) #for GO terms

    print(r.status_code)
    result = r.json()
    genesets = result['objects']
    simplified_genesets = []
    count = 0
    for geneset in genesets:
        count += 1
        if geneset['tip'] == None:
            continue
        if (len(geneset['tip']['genes']) > 100) or (len(geneset['tip']['genes']) < 5):
            continue
        else:
            new_geneset = {}
            new_geneset['title'] = geneset['title']
            new_geneset['genes'] = geneset['tip']['genes']
            simplified_genesets.append(new_geneset)        
    if count == limit:
        print 'Need to increase the limit.'
    return simplified_genesets


if __name__ == '__main__':
    out_file = sys.argv[1]
    target = sys.argv[2]
    out_fh = open(out_file,'w')
    anno_terms = get_terms(target=target)
    for term in anno_terms:
        out_fh.write(term['title']+'\t'+ str(len(term['genes']))+'\t'+';'.join([str(x) for x in term['genes']])+'\n')
    out_fh.close()
