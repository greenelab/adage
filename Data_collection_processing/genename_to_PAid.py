"""
Convert gene name to PA ID.
"""
import sys

id_symbol = {}
gi_fh = open('./Data_collection_processing/Pseudomonas_aeruginosa_PAO1.gene_info')
gi_fh.next()#skip header
for line in gi_fh:
    toks = line.strip().split('\t')
    gid = toks[2]
    symbol = toks[3]
    id_symbol[gid] = symbol
gi_fh.close()

def rename_pcl(in_file, out_file):

    in_fh = open(in_file,'r')
    out_fh = open(out_file,'w')
    out_fh.write(in_fh.readline())
    for line in in_fh:
        toks = line.strip().split('\t')
        gene1 = toks[0].strip()
        try:
            gid1 = id_symbol[gene1]
        except KeyError:
            gid1 = gene1
            print gene1+' not converted!'
        out_fh.write(gid1 +'\t'+ '\t'.join(map(lambda x: x.strip(),toks[1:])) + '\n')

    out_fh.close()
    in_fh.close()

rename_pcl(sys.argv[1],sys.argv[2])