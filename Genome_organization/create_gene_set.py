'''
Convert pickle file that stores operons into the gmt file for GSEA.
'''

import pickle
import sys

def create_geneset(operon_file, out_file):
    out_fh = open(out_file,'w')
    operon_fh = open(operon_file,'r')

    for line in operon_fh:
        operon = line.strip().split('\t')
        oper_name = ';'.join(operon)
        out_fh.write(oper_name+'\t'+'na'+'\t'+'\t'.join(operon)+'\n')

create_geneset(sys.argv[1],sys.argv[2])