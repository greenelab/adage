"""
This code is adapted from Deep Learning Tutorials http://deeplearning.net/tutorial/

Copyright (c) 2008-2013, Theano Development Team All rights reserved.
Copyright (c) 2015, Jie Tan All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
Neither the name of Theano nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Usage:
        SdA_train.py <input-file> <skip-col> <net-structure> <batch-size> <epoch-size> <corruption-level> <learning-rate> [--seed1=<SEED1>] [--seed2=<SEED2>]
        SdA_train.py -h | --help

Options:
        -h --help               Show this screen.
        <input-file>            File path to the microarray file that needed to be analyzed
        <skip-col>              int, the number of column to be skipped between the first gene ID column and the first experimental column
        <net-structure>         A list of ints separated by comma
        <batch-size>            The size of a minibatch
        <epoch-size>            The number of epoch to do training
        <corruption-level>      The corruption percentage in denoised autoencoders
        <learning-rate>         The learning rate to be used during training
        --seed1 = <SEED1>       Random seed for training  [default: 123]
        --seed2 = <SEED2>       Random seed for permuting order of samples [default:123]


"""

import cPickle
import gzip
import os
import sys
import time
import logging

from docopt import docopt
import numpy
import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams

from dA import dA
sys.path.insert(0,'Data_collection_processing/')
from pcl import PCLfile


class SdA(object):

    def __init__(self, numpy_rng, theano_rng=None, n_ins=100,
                 hidden_layers_sizes=[500, 500], n_outs=10,
                 corruption_levels=[0.1, 0.1]):
        """ This class is made to support a variable number of layers.
        :type numpy_rng: numpy.random.RandomState
        :param numpy_rng: numpy random number generator used to draw initial
                    weights

        :type theano_rng: theano.tensor.shared_randomstreams.RandomStreams
        :param theano_rng: Theano random generator; if None is given one is
                           generated based on a seed drawn from `rng`

        :type n_ins: int
        :param n_ins: dimension of the input to the sdA

        :type n_layers_sizes: list of ints
        :param n_layers_sizes: intermediate layers size, must contain
                               at least one value

        :type n_outs: int
        :param n_outs: dimension of the output of the network

        :type corruption_levels: list of float
        :param corruption_levels: amount of corruption to use for each
                                  layer
        """

        self.dA_layers = []
        self.params = []
        self.n_layers = len(hidden_layers_sizes)

        assert self.n_layers > 0

        if not theano_rng:
            theano_rng = RandomStreams(numpy_rng.randint(2 ** 30))
        # allocate symbolic variables for the data
        self.x = T.matrix('x')  

        for i in xrange(self.n_layers):
            # construct the DA layer

            # the size of the input is either the number of hidden units of
            # the layer below or the input size if we are on the first layer
            if i == 0:
                input_size = n_ins
            else:
                input_size = hidden_layers_sizes[i - 1]

            # the input to this layer is either the activation of the hidden
            # layer below or the input of the SdA if you are on the first
            # layer
            if i == 0:
                layer_input = self.x
            else:
                layer_input = self.dA_layers[-1].output

            # Construct a denoising autoencoder
            dA_layer = dA(numpy_rng=numpy_rng,
                          theano_rng=theano_rng,
                          input=layer_input,
                          n_visible=input_size,
                          n_hidden=hidden_layers_sizes[i])
            self.dA_layers.append(dA_layer)

    def training_functions(self, train_set_x, batch_size):
        ''' Generates a list of functions, each of them implementing one
        step in trainnig the dA corresponding to the layer with same index.
        The function will require as input the minibatch index, and to train
        a dA you just need to iterate, calling the corresponding function on
        all minibatch indexes.

        :type train_set_x: theano.tensor.TensorType
        :param train_set_x: Shared variable that contains all datapoints used
                            for training the dA 

        :type batch_size: int
        :param batch_size: size of a [mini]batch    
        '''

        index = T.lscalar('index')  # index to a minibatch
        corruption_level = T.scalar('corruption')  # % of corruption to use
        learning_rate = T.scalar('lr')  # learning rate to use
        # number of batches
        n_batches = train_set_x.get_value(borrow=True).shape[0] / batch_size
        # begining of a batch, given `index`
        batch_begin = index * batch_size
        # ending of a batch given `index`
        batch_end = batch_begin + batch_size

        train_fns = []
        for dA in self.dA_layers:
            # get the cost and the updates list
            cost, updates = dA.get_cost_updates(corruption_level,
                                                learning_rate)
            # compile the theano function
            fn = theano.function(inputs=[index,
                              theano.Param(corruption_level, default=0.2),
                              theano.Param(learning_rate, default=0.1)],
                                 outputs=cost,
                                 updates=updates,
                                 givens={self.x: train_set_x[batch_begin:
                                                             batch_end]})
            # append `fn` to the list of functions
            train_fns.append(fn)

        return train_fns

    def return_activity(self, train_set_x):
        '''Given an input, this function returns the activity
        value of all the nodes in each hidden layer.'''
        
        activity_each_layer = []
        index = T.lscalar('index')  # index to a sample
        
        for dA in self.dA_layers:
            activity_fn = theano.function(inputs=[index],outputs = dA.output,
                                          givens={self.x: train_set_x[index:(index+1)]})
            activity_each_layer.append(activity_fn)
        return activity_each_layer

    def return_raw_activity(self, train_set_x):
        '''Given an input, this function returns the raw activity
        value of all the nodes in each layer.'''
        
        raw_activity_each_layer = []
        index = T.lscalar('index')  # index to a sample
        
        for dA in self.dA_layers:
            raw_activity_fn = theano.function(inputs=[index],outputs = dA.raw_output,
                                              givens={self.x: train_set_x[index:(index+1)]})
            raw_activity_each_layer.append(raw_activity_fn)
        return raw_activity_each_layer

    def return_network(self):
        '''This function returns weight matrix and bias vectors of each hidden layer in the 
        final network after training.'''

        weights_all_layer = []
        bias_all_layer = []
        bias_prime_all_layer = []

        for dA_layer in self.dA_layers:
            weight = dA_layer.W.get_value(borrow = True)
            bias = dA_layer.b.get_value(borrow = True)
            bias_prime = dA_layer.b_prime.get_value(borrow = True)
            weights_all_layer.append(weight)
            bias_all_layer.append(bias)
            bias_prime_all_layer.append(bias_prime)

        return weights_all_layer, bias_all_layer, bias_prime_all_layer

def train_SdA(training_epochs=15, train_lr=0.001, data_file = None, skip_col = 2,
              batch_size=1, random_seed_1 = 89677, random_seed_2 = 123,net_structure = [1000,1000,1000], 
              corruption_levels = [.1, .2, .3], output_file = None, net_file = None):


    logging.basicConfig(filename = output_file.replace('activity_SdA.txt', 'SdA.log'), level= logging.INFO)
    logging.info('Training the dataset:' + data_file)
    logging.info('The structure of the networks:'+ str(net_structure))
    logging.info('Training epoches:'+str(training_epochs)+'\t'+'Batch size:'+str(batch_size)+'\t'+'Learning rate:'+str(train_lr)+'\t'+'Corruption levels:'+str(corruption_levels)+'\n'
        +'Random seed for training:'+str(random_seed_1)+'\t'+ 'Ramdom seed for permuting sample order:'+str(random_seed_2))
    
    datasets = PCLfile(data_file, skip_col)    
    train_set_x, sample_id = datasets.get_permuted_sample(seed = random_seed_2)#Permute the order of samples using random_seed_2
    print '... finish reading the data'

    train_set_x = theano.shared(train_set_x,borrow=True)

    # compute number of minibatches for training
    train_size = train_set_x.get_value(borrow=True).shape[0]
    n_train_batches = train_size / batch_size

    # numpy random generator
    numpy_rng = numpy.random.RandomState(random_seed_1)

    # the number of input nodes
    input_node = len(datasets.id_list)

    print '... building the model'
    # construct the stacked denoising autoencoder class
    sda = SdA(numpy_rng=numpy_rng, n_ins= input_node,
              hidden_layers_sizes= net_structure)

    #########################
    # TRAINING THE MODEL #
    #########################
    print '... getting the training functions'
    training_fns = sda.training_functions(train_set_x=train_set_x,
                                          batch_size=batch_size)

    print '... training the model'
    start_time = time.clock()
    ## Train layer-wise
    corruption_levels = corruption_levels
    for i in xrange(sda.n_layers):
        # go through training epochs
        for epoch in xrange(training_epochs):
            # go through the training set
            c = []
            for batch_index in xrange(n_train_batches):
                c.append(training_fns[i](index=batch_index,
                         corruption=corruption_levels[i],
                         lr=train_lr))
            print 'Training layer %i, epoch %d, cost ' % (i, epoch),
            print numpy.mean(c)
            logging.info('Training layer %i, epoch %d, cost %f ' % (i, epoch, numpy.mean(c) ))

    end_time = time.clock()

    logging.info('The training code for file ' + os.path.split(__file__)[1] + ' ran for %.2fm' % ((end_time - start_time) / 60.))
    print '... training finished.'

    ##############################################################
    # Return the final activity value and raw activity value
    # for each node of each input sample 
    ##############################################################
    output_fh = open(output_file,'w')
    raw_output_fh = open(output_file.replace('activity','rawActivity'),'w')
    each_layer_output = sda.return_activity(train_set_x=train_set_x)
    each_layer_raw_output = sda.return_raw_activity(train_set_x=train_set_x)
    for i in xrange(sda.n_layers):
        output_fh.write('layer %i \n' %(i+1))
        raw_output_fh.write('layer %i \n' %(i+1))
        for train_sample in xrange(train_size):
            node_activation = each_layer_output[i](train_sample)
            node_raw_activation = each_layer_raw_output[i](train_sample)
            output_fh.write(sample_id[train_sample]+'\t')
            raw_output_fh.write(sample_id[train_sample]+'\t')
            numpy.savetxt(output_fh, node_activation, fmt= '%.8f', delimiter= '\t') 
            numpy.savetxt(raw_output_fh, node_raw_activation, fmt= '%.8f', delimiter= '\t') 


    ##############################################################
    # Return weight matrix and bias vectors of the final network #
    ##############################################################
    net_file = open(net_file,'w')
    weight_output, bias_output, bias_prime_output = sda.return_network()
    for i in xrange(len(weight_output)):
        net_file.write('layer %i \n' %(i+1))
        net_file.write('weight matrix \n')
        numpy.savetxt(net_file, weight_output[i], fmt= '%.8f', delimiter = '\t') 
        net_file.write('hidden bias vector \n')
        numpy.savetxt(net_file, bias_output[i], fmt= '%.8f', delimiter = '\t')
        net_file.write('visible bias vector \n')
        numpy.savetxt(net_file, bias_prime_output[i], fmt= '%.8f', delimiter = '\t')


if __name__ == '__main__':

    arguments = docopt(__doc__, version=None)
    input_file = arguments['<input-file>']
    network_stru = [int(x) for x in arguments['<net-structure>'].strip().split(',')]
    corrupt_levels = [float(x) for x in arguments['<corruption-level>'].strip().split(',')]
    train_SdA(training_epochs=int(arguments['<epoch-size>']) ,train_lr=float(arguments['<learning-rate>']), data_file= input_file, skip_col= int(arguments['<skip-col>']),batch_size=int(arguments['<batch-size>']), random_seed_1 = int(arguments['--seed1']), random_seed_2 = int(arguments['--seed2']), net_structure = network_stru, corruption_levels = corrupt_levels, output_file = input_file.replace('.pcl', '')+'_'+ arguments['<net-structure>']+'_batch'+arguments['<batch-size>'] + '_epoch' + arguments['<epoch-size>'] + '_corrupt' + arguments['<corruption-level>'] + '_lr' +arguments['<learning-rate>']  + '_seed1_' + arguments['--seed1'] + '_seed2_' + arguments['--seed2']  + "_activity_SdA.txt", net_file = input_file.replace('.pcl', '')+'_'+ arguments['<net-structure>'] +'_batch'+ arguments['<batch-size>'] + '_epoch' + arguments['<epoch-size>'] + '_corrupt' + arguments['<corruption-level>'] +  '_lr' +arguments['<learning-rate>'] + '_seed1_' + arguments['--seed1']+ '_seed2_' + arguments['--seed2'] +"_network_SdA.txt")
