'''
Run Gene Set Enrichment Analysis for each node.
Operons are used as gene sets and the weight 
vector connect genes to a node is used as ranked
gene list.
'''

import sys
sys.path.insert(0,'Data_collection_processing/')
from pcl import PCLfile
import numpy
import os

def read_weight_matrix(data_file, network_file):
    '''
    This function reads weight matrix from network structure file
    and the corresponding gene id from data file
    '''
    datasets = PCLfile(data_file, skip_col=0)
    gene_id = datasets.id_list
    network_fh = open(network_file,'r')     
    input_size = len(gene_id)   
    network_fh.next() # skip the layer count line
    network_fh.next() # skip 'weight matrix' line
    W = []
    input_count = 0
    for line in network_fh:
        line = line.strip().split('\t')
        W.append(line)
        input_count += 1
        if input_count == input_size:
            break
    W = numpy.array(W, dtype= float)
    return gene_id, W

data_file = sys.argv[1]
net_file = sys.argv[2]
gene_id, weight = read_weight_matrix(data_file, net_file)
for node in xrange(weight.shape[1]):
    #prepare a rnk file for each node
    out_file = 'Genome_organization/GSEA/Node'+ str(node+1)+'.rnk' 
    out_fh = open(out_file,'w')
    for gene in xrange(weight.shape[0]):
        out_fh.write(gene_id[gene]+'\t'+str(weight[gene,node])+'\n')
    out_fh.close()
    
    #Run GSEA using the weight vector of each node.
    label = 'Node'+str(node+1)
    os.system('java -Xmx1024m -cp Genome_organization/gsea2-2.0.14.jar xtools.gsea.GseaPreranked -gmx Genome_organization/operon_set_3.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk '+ out_file+' -scoring_scheme weighted -rpt_label '+ label+' -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 3 -zip_report false -out Genome_organization/GSEA/ -gui false')    