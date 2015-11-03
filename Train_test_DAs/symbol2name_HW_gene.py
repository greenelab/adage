'''
For each high-weight gene file, change their gene symbols to gene names.
'''

import sys
import glob
import numpy

#Build a dictionary that stores gene symbol and its corresponding gene name.
symbol_id = {}
gi_fh = open('Data_collection_processing/Pseudomonas_aeruginosa_PAO1.gene_info')
gi_fh.next()#skip header
for line in gi_fh:
    toks = line.strip().split('\t')
    gid = toks[2]
    symbol = toks[3]
    symbol_id[symbol] = gid
gi_fh.close()

def symbol_to_name(geneList_folder, out_folder):
    '''
    geneList_folder: the folder stores high-weight gene files for each node
    out_folder: the folder stores renamed high-weight gene files for each node
    '''
    #Count the number of files in the folder
    node_num =  len(glob.glob(geneList_folder + '/Node*.txt'))

    for node in range(0,node_num):
        gene_fh = open(geneList_folder+'/Node'+str(node+1)+'.txt','r')    
        gene_fh.next()
        out_fh = open(out_folder + '/Node'+str(node+1)+'.txt','w')
        out_fh.write('gene_id\tgene_symbol\tweight\n')
        for line in gene_fh:
            toks = line.strip().split('\t')
            symbol = toks[0]
            #Change gene symbol to gene name
            try:
                gid = symbol_id[symbol]
                out_fh.write(gid + '\t' + '\t'.join(toks[:]) + '\n')
            except KeyError:
                #if the symbol cannot be found in the dictionary, keep the original symbol
                out_fh.write(symbol + '\t' + '\t'.join(toks[:]) + '\n') 
        node += 1
        out_fh.close()
        gene_fh.close()

symbol_to_name(geneList_folder = sys.argv[1], out_folder = sys.argv[2])


