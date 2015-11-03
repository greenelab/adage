"""
Usage:
        SdA_test.py <data-file> <skip-col> <network-file> <net-structure>
        SdA_test.py -h | --help

Options:
        -h --help               Show this screen.
        <data-file>             File path of the test sample file
        <skip-col>              int, the number of column to be skipped between the first gene ID column and the first experimental column
        <network-file>          Existing network that are used in the test
        <net-structure>         The net structure used by the exisiting network. A list of ints separated by comma 
"""


'''
This code tests the test set on existing network that are previously trained 
It outputs activity values and raw activity values for each node of each sample in the test set.
'''

import numpy
import sys
sys.path.insert(0,'Data_collection_processing/')
from pcl import PCLfile
from docopt import docopt

def logit(x):
    '''logistic function'''
    y = 1/(1 + numpy.exp(-x))
    return y

def SdA_test(data_file, skip_col, network_file, net_structure):

    activity_file = data_file.replace('.pcl','_activity') + '_with_' + network_file.split('/')[-1] 
    raw_activity_file = data_file.replace('.pcl','_rawActivity') + '_with_' + network_file.split('/')[-1] 
    datasets = PCLfile(data_file, skip_col)

    input_data = datasets.get_sample()
    input_data = numpy.matrix(input_data)
    sample_id = datasets.sample_list    

    network_fh = open(network_file,'r') 
    input_size = input_data.shape[1] #input_size is the number of genes     
    layer_para = []
    for each_layer in xrange(len(net_structure)):
        network_fh.next() # skip the layer count line
        network_fh.next() # skip 'weight matrix' line
        
        #Get the weight matrix
        W = []
        input_count = 0
        for line in network_fh:
            line = line.strip().split('\t')
            W.append(line)
            input_count += 1
            if input_count == input_size: #when it reach the number of genes
                break
        W = numpy.matrix(W, dtype= float)
        network_fh.next() # skip 'hidden bias vector' line

        #Get the bias vector of hidden layer
        h_bias = []
        output_count = 0
        for line in network_fh:
            line = line.strip().split('\t')
            h_bias.append(line)
            output_count += 1
            if output_count == int(net_structure[each_layer]): #when it reach the number of nodes in hidden layer
                break
        h_bias = numpy.matrix(h_bias, dtype = float)
        network_fh.next() # skip 'visible bias vector' line

        #Get the bias vector of visible output layer
        v_bias = []
        input_count = 0
        for line in network_fh:
            v_bias.append(line)
            input_count += 1
            if input_count == input_size:
                break
        v_bias = numpy.matrix(v_bias, dtype = float)

        layer_para.append((W,h_bias)) #Weight matrix and hidden bias vector are enough for calculating activities
        input_size = net_structure[each_layer]  


    activity_fh = open(activity_file, 'w')
    raw_activity_fh = open(raw_activity_file,'w')
    for each_layer in xrange(len(net_structure)):
        activity_fh.write('layer %i \n' %(each_layer+1))
        raw_activity_fh.write('layer %i \n' %(each_layer+1))
        temp = input_data * layer_para[each_layer][0]  #calculate raw activity, before sigmoid transformation
        output = logit(temp+ layer_para[each_layer][1].T) #calculate activity, after sigmoid transformation
        for each_sample in xrange(output.shape[0]):
            activity_fh.write(sample_id[each_sample]+'\t')
            raw_activity_fh.write(sample_id[each_sample]+'\t')
            numpy.savetxt(activity_fh, output[each_sample,], fmt= '%.8f', delimiter= '\t')
            numpy.savetxt(raw_activity_fh, temp[each_sample,], fmt= '%.8f', delimiter= '\t')
        input_data = output 


arguments = docopt(__doc__, version=None)
data_file = arguments['<data-file>']
network_file = arguments['<network-file>']
network_stru = [int(x) for x in arguments['<net-structure>'].strip().split(',')]
skip_col= int(arguments['<skip-col>'])
SdA_test(data_file,skip_col, network_file, network_stru)
