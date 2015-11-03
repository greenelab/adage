'''
Combine all the significant associated operons with each node into one file.
'''

import os
import glob
import sys

FDR = float(sys.argv[1]) #significance cutoff 
net_size = 50
out_file = 'Genome_organization/GSEA'+'_FDR_'+str(FDR)+'.txt'
out_fh = open(out_file,'w')
out_fh.write('node\toperon\tq_value\n')

#loop over each node
for node in xrange(net_size):

    #operons enriched with positive weights
    input_file = glob.glob('Genome_organization/GSEA/Node'+str(node+1)+'.*/gsea_report_for_na_pos_*.xls')[0]
    input_fh = open(input_file,'r')
    input_fh.next()
    for line in input_fh:
        tok = line.strip().split('\t')
        if float(tok[6]) <= FDR: #lower than cutoff
            out_fh.write('node'+str(node+1)+'\t'+tok[0]+'\t'+tok[6]+'\n')
        else:
            break

    #operons enriched with negative weights
    input_file = glob.glob('Genome_organization/GSEA/Node'+str(node+1)+'.*/gsea_report_for_na_neg_*.xls')[0]
    input_fh = open(input_file,'r')
    input_fh.next()
    for line in input_fh:
        tok = line.strip().split('\t')
        if float(tok[6]) <= FDR:
            out_fh.write('node'+str(node+1)+'\t'+tok[0]+'\t'+tok[6]+'\n')
        else:
            break
            
out_fh.close()