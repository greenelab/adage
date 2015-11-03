'''
This code is used to find node whose high-weight genes are over-representated 
in a gene list (a biological pathway or a TF's downstream targets).

'''

import glob
from scipy import stats
import sys
sys.path.insert(0,'Data_collection_processing/')
from pcl import PCLfile
from statsmodels.stats.multitest import *

def find_enriched_node(geneList_folder = None, data_file= None, gold_std = None, out_file = None):
    '''
    geneList_folder: the folder stores high-weight gene files for each node
    data_file: the microarray file with gene genes
    gold_std: gold standard file contains a list of genes
    out_file: the output file that stores each node and its corresponding q value
    '''

    datasets = PCLfile(data_file, skip_col=0)
    gene_id = datasets.id_list  

    gold_fh = open(gold_std,'r')
    gold_set = []
    for line in gold_fh:
        gene = line.strip().split('\t')[0]
        gold_set.append(gene)

    p_all_node = []
    geneList_files =  glob.glob(geneList_folder + '/Node*.txt') #Get all the high-weight gene files under the geneList_folder
    for i in xrange(len(geneList_files)):
        gene_fh = open(geneList_folder + '/Node'+str(i+1)+'.txt','r')
        gene_fh.next()
        geneset = [] #geneset stores the high-weight gene for a node
        for line in gene_fh:
            gene = line.strip().split('\t')[0]
            geneset.append(gene)
        #Build the contengency table
        all_overlap_genes = set(gold_set).intersection(set(gene_id))
        selected_overlap_genes = set(gold_set).intersection(set(geneset))
        a = len(selected_overlap_genes)
        b = len(all_overlap_genes) - len(selected_overlap_genes)
        c = len(geneset) - len(selected_overlap_genes)
        d = len(gene_id) -a -b -c
        table = [[a,b],[c,d]]
        #Calculate p-value using fisher exact test
        oddsratio, pvalue = stats.fisher_exact(table)
        p_all_node.append(pvalue)   

    #multiple hypothesis correction
    result_adj_pvalue = multipletests(p_all_node, alpha=0.05, method='fdr_bh')[1] 
    all_node = ['Node'+str(x+1) for x in xrange(50)]
    #find the node with lowest q value, write the output
    qvalue_small = 1 
    node_small = None
    out_fh = open(out_file, 'w')
    for node, qvalue in zip(all_node, result_adj_pvalue):
        out_fh.write(node+'\t'+str(qvalue)+'\n')
        if qvalue < qvalue_small:
            qvalue_small = qvalue
            node_small = node
    print node_small+' is most significantly associated with this gene set with a q value of '+str(qvalue_small)


find_enriched_node(geneList_folder = sys.argv[1], data_file= sys.argv[2], gold_std = sys.argv[3], out_file= sys.argv[4])
