'''
For each hidden node, write its high-weight genes into a file.
"High-weight" is defined by higher than a threshold. Here we use
two or three standard deviations.
'''

import numpy
import sys
sys.path.insert(0,'Data_collection_processing/')
from pcl import PCLfile
import pickle

def read_weight_matrix(data_file, network_file):
    '''
    This function read weight matrix from network structure file
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
    return W, gene_id



thre = int(sys.argv[1]) #the number of standard deviations to use as cutoff for high-weight genes
data_file = sys.argv[2]
network_file = sys.argv[3]
output_folder = sys.argv[4]

weight, gene_id = read_weight_matrix(data_file,network_file)
for node in xrange(weight.shape[1]):
    aver_weight = numpy.mean(weight[:,node])    
    std_weight = numpy.std(weight[:,node])
    cutoff_up = aver_weight + thre*std_weight
    cutoff_down = aver_weight - thre*std_weight
    out_fh = open(output_folder+'/Node'+str(node+1)+'.txt','w')
    out_fh.write('gene\tweight\n')
    for gene in xrange(weight.shape[0]):
        if weight[gene,node] >= cutoff_up or weight[gene,node] <= cutoff_down: #gene weight is outside of the cutoffs
            out_fh.write(gene_id[gene]+'\t'+str(weight[gene,node])+'\n')
    out_fh.close()





