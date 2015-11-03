'''
This script create a file that stores each dataset's name and the samples contained in the dataset.
'''

import os

pcl_folder = 'Data_collection_processing/data/pcls'
out_fh = open('Data_collection_processing/datasets_list.txt','w')

pcl_list = os.listdir(pcl_folder)
for file in pcl_list:
    pcl_fh = open(os.path.join(pcl_folder,file),'r')
    header = pcl_fh.readline()
    samples = header.strip().split('\t')
    out_fh.write(file.strip('.pcl')+'\t'+';'.join(samples)+'\n')