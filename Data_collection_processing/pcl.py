'''
This object reads PCL files and prepare the microarray data as training set to DAs.
The input training vector can either be a gene's expression value over all sampels, or one 
microarray sample with all genes' expression value. To feed into DAs, the 
standard input dataset is a two-dimensional array with each row as a training sample.
'''

import numpy

class PCLfile(object):

    def __init__(self, dataset, skip_col=2):
        '''
        type dataset: string
        param dataset: path to the pcl file
        type skip_col: int
        param skip_col: the number of colunms to skip between the first gene ID column and the first 
                        experimental column.
        '''
        
        try:
            dataset_fh = open(dataset,'r')
        except IOError:
            print "Error, file not found."
            
        self.data_matrix = []
        self.id_list = []
       
        line_count = 0
        for line in dataset_fh:
            if line_count == 0:
                self.sample_list = line.rstrip().split('\t')[(skip_col+1):] #This stores samples' names
                line_count +=1
                continue
                
            line_new = line.strip().split('\t')
            self.data_matrix.append(line_new[(skip_col+1):]) #This extract microarray data with gene in rows, sample in columns. 
            self.id_list.append(line_new[0]) #This stores each gene's ID
        
        self.data_matrix = numpy.array(self.data_matrix, dtype = numpy.float64) #Convert data_matrix to a numpy array
        
    #Normalize every row linearly so that the min is 0 and max is 1
    #This directly change the self.data_matrix
    def zero_one_normalization(self):

        for i in xrange(self.data_matrix.shape[0]): #'shape' return the dimension of the the matrix and shape[0] return the first dimension which is the row.
            row_minimum = self.data_matrix[i,:].min()
            row_maximum = self.data_matrix[i,:].max()
            row_range = row_maximum - row_minimum
            self.data_matrix[i,:] = (self.data_matrix[i,:] - row_minimum)/row_range

    def zero_one_normalization_sample(self):

        for i in xrange(self.data_matrix.shape[1]): #'shape' return the dimension of the the matrix and shape[0] return the first dimension which is the row.
            row_minimum = self.data_matrix[:,i].min()
            row_maximum = self.data_matrix[:,i].max()
            row_range = row_maximum - row_minimum
            self.data_matrix[:,i] = (self.data_matrix[:,i] - row_minimum)/row_range
            
    #Normalize the data, one row at a time, by converting each value to a row-relative Z score.
    #This directly change the self.data_matrix
    def z_normalization(self):
    
        for i in xrange(self.data_matrix.shape[0]):
            mean = numpy.mean(self.data_matrix[i,:])
            standev = numpy.std(self.data_matrix[i,:])
            self.data_matrix[i,:] = (self.data.matrix[i,:] - mean) / standev

    def z_normalization_sample(self):
    
        for i in xrange(self.data_matrix.shape[1]):
            mean = numpy.mean(self.data_matrix[:,i])
            standev = numpy.std(self.data_matrix[:,i])
            self.data_matrix[:,i] = (self.data.matrix[:,i] - mean) / standev

    def logistic_normalization(self):

        for i in xrange(self.data_matrix.shape[0]):
            self.data_matrix[i,:] = 1.0 / (1.0+ numpy.exp(-self.data_matrix[i,:]))

        
    #This function returns a matrix with each gene in a row
    def get_gene(self):
        
        return self.data_matrix
        
    #This function returns a matrix with each sample in a row
    def get_sample(self):
        
        return self.data_matrix.T

    #This function permutes samples. It returns a data matrix with the order of samples being permuted.
    def get_permuted_sample(self,seed=123):
        
        transposed = self.data_matrix.T #After matrix transpose, each row represents one sample from the microarray data.
        if seed == 0:
            return transposed, self.sample_list
        else:
            numpy.random.seed(seed)
            numpy.random.shuffle(transposed) #numpy.random.shuffle only shuffles the array along the 
                                        #first index of a multi-dimensional array, which is the row here.
            numpy.random.seed(seed)
            numpy.random.shuffle(self.sample_list)
            return transposed , self.sample_list

    def permuted_gene_order(self, seed=123):
        numpy.random.seed(seed)
        numpy.random.shuffle(self.data_matrix)
        numpy.random.seed(seed)
        numpy.random.shuffle(self.id_list)

    #This function writes a PCLfile object to a text file
    def write_pcl(self, outputPath):
        try:
            outputFileHandle = open(outputPath, 'w')
        except IOError:
            print "Was not able to open the output file"
            return False
            
        #First write the header
        outputFileHandle.write('Gene_symbol\t')
        header = '\t'.join(map(str,self.sample_list))
        outputFileHandle.write(header)
        outputFileHandle.write('\n')
    
        #Now write the gene values
        for i in range(self.data_matrix.shape[0]):
            geneID = self.id_list[i]
            geneValue = self.data_matrix[i,:]
            outputFileHandle.write(geneID + '\t' + '\t'.join(map(str, geneValue)))
            outputFileHandle.write('\n')
        
        outputFileHandle.close()
        return True
        