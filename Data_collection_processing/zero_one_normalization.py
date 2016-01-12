'''
Linearly scale the expression range of one gene to be between 0 and 1.
If a reference dataset is provided, then the scaling of one gene in the
target dataset in done using the minimun and range of that gene in the 
reference dataset.
'''

import sys
import argparse
sys.path.insert(0,'Data_collection_processing/')
from pcl import PCLfile

parser = argparse.ArgumentParser(description="Linearly scale the expression range of one gene to be between 0 and 1. If a reference dataset is provided, then the scaling of one gene in the target dataset in done using the minimun and range of that gene in the reference dataset.")
parser.add_argument('tar',help='the target file for zero one normalization')
parser.add_argument('out',help='the output file after zero one normalization')
parser.add_argument('ref',help='the reference file. If reference file is \'None\', then zero one normalization will be done based on target file itself.')
args=parser.parse_args()

def zero_one_normal(tar=None, out=None, ref=None):
    '''
    tar: the target file for zero one normalization
    out: the output file after zero one normalization
    ref: the reference file. If reference file is 'None', 
         then zero one normalization will be done based on 
         target file itself.
    '''

    if ref=='None':
        tar_data=PCLfile(tar, skip_col = 0)
        tar_data.zero_one_normalization()
        tar_data.write_pcl(out)
    else:
        ref_data = PCLfile(ref, skip_col=0)
        tar_data = PCLfile(tar, skip_col=0)
        for i in xrange(ref_data.data_matrix.shape[0]): 
            row_minimum = ref_data.data_matrix[i,:].min()
            row_maximum = ref_data.data_matrix[i,:].max()
            row_range = row_maximum - row_minimum
            tar_data.data_matrix[i,:] = (tar_data.data_matrix[i,:] - row_minimum)/row_range
            #bound the values to be between 0 and 1
            tar_data.data_matrix[i,:] = [0 if x< 0 else x for x in tar_data.data_matrix[i,:]] 
            tar_data.data_matrix[i,:] = [1 if x> 1 else x for x in tar_data.data_matrix[i,:]]
        tar_data.write_pcl(out)

zero_one_normal(tar=args.tar, out=args.out, ref=args.ref)
