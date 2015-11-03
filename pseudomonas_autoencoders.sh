#############################################
#Section one: collecting and processing data#
#############################################

    ############
    ##Path One##
    ############
    #If you want to build an up-to-date Pseudomonas expression compendium, you can download datasets and process them as follows:

    #Download pseudonomas datasets from ArrayExpress
    mkdir -p Data_collection_processing/data/zips
    #For mac users, you need to install wget first or use curl instead
    for x in `python Data_collection_processing/get_pseudo.py`; do wget -N -P Data_collection_processing/data/zips/ $x ; done   

    #Unzip samples in all datasets into one folder and also unzip samples in one dataset into individual folders.
    mkdir -p Data_collection_processing/data/cels/all-pseudomonas
    for x in Data_collection_processing/data/zips/*; do unzip -n $x -d Data_collection_processing/data/cels/all-pseudomonas; done
    for x in Data_collection_processing/data/zips/*; do mkdir -p Data_collection_processing/data/cels/`basename -s .raw.1.zip $x`; unzip -n $x -d Data_collection_processing/data/cels/`basename -s .raw.1.zip $x`; done    

    #Process each dataset into pcl file. Also process samples in all datasets into one expression compendium.
    mkdir -p Data_collection_processing/data/pcls/
    #The following code requires installing R packages from Bioconductor: affy, affyio, AnnotationDbi, paeg1acdf
    for x in Data_collection_processing/data/cels/*; do R --no-save --args $x Data_collection_processing/data/pcls/`basename $x`.pcl < Data_collection_processing/ProcessToPCL.R; done

    #Create a file to store all datasets names and their sample names
    python Data_collection_processing/create_dataset_list.py

    #Remove controls and only keep genes starting with PA. The first argument takes input file (combined expression compendium) and second argument specifies output file.
    python Data_collection_processing/remove_control.py Data_collection_processing/data/pcls/all-pseudomonas.pcl Data_collection_processing/Pa_compendium.pcl

    #Process test sets that are not included in the expression compendium. Let's use the genome hybridization test set and Anr test set as examples. You can also provide your own dataset here.
    R --no-save --args Data_collection_processing/test_sets/Genome-hybs < Data_collection_processing/process_test_set.R 
    R --no-save --args Data_collection_processing/test_sets/Anr < Data_collection_processing/process_test_set.R 
    
    #Remove controls and only keep genes in the test set
    python Data_collection_processing/remove_control.py Data_collection_processing/test_sets/Genome-hybs.pcl Data_collection_processing/Genome-hybs-gene.pcl
    python Data_collection_processing/remove_control.py Data_collection_processing/test_sets/Anr.pcl Data_collection_processing/Anr-gene.pcl

    ##########################
    # process RNAseq testset #
    ##########################
    #convert gene names to PA IDs
    python Data_collection_processing/genename_to_PAid.py Data_collection_processing/test_sets/Anr-RNAseq/GSE68534_Processed_PAO1.txt Data_collection_processing/test_sets/Anr-RNAseq/GSE68534_Processed_PAO1_renamed.txt
    python Data_collection_processing/genename_to_PAid.py Data_collection_processing/test_sets/Anr-RNAseq/GSE68534_Processed_J215.txt Data_collection_processing/test_sets/Anr-RNAseq/GSE68534_Processed_J215_renamed.txt
    #process RNAseq data using tdm r package
    Rscript Data_collection_processing/process_rnaseq_testset.R Data_collection_processing/test_sets/Anr-RNAseq/GSE68534_Processed_PAO1_renamed.txt Data_collection_processing/Pa_compendium.pcl Data_collection_processing/Anr_RNAseq_PAO1_gene.pcl
    Rscript Data_collection_processing/process_rnaseq_testset.R Data_collection_processing/test_sets/Anr-RNAseq/GSE68534_Processed_J215_renamed.txt Data_collection_processing/Pa_compendium.pcl Data_collection_processing/Anr_RNAseq_J215_gene.pcl

    ############
    ##Path Two##
    ############
    #If you want to reproduce our analysis, we only included datasets that are available before 02/22/2014 in the compendium. You can start from here using the already processed expression compendium provided as Data_collection_processing/Pa_compendium_02.22.2014.pcl and the processed test sets Data_collection_processing/Genome-hybs-gene.pcl and Data_collection_processing/Anr-gene.pcl.
    #Note that test set processing also depends on the compendium. If you process the test set with an up-to-date compendium, then the expression values in the resulting test set might be slightly different.


#If you are using an up-to-date compendium, replace Pa_compendium_02.22.2014.pcl with the up-to-date compendium in the following analyses.
#Note that the node order will change if you train a new DA on an updated compendium.

#To prepare for DA training, each gene expression vector is linearly normalized to the range between 0 and 1. 
python Data_collection_processing/zero_one_normalization.py Data_collection_processing/Pa_compendium_02.22.2014.pcl Train_test_DAs/train_set_normalized.pcl None

#Test set is normalzied using the same minimums and ranges used above.
python Data_collection_processing/zero_one_normalization.py Data_collection_processing/Genome-hybs-gene.pcl Train_test_DAs/Genome-hybs_normalized.pcl Data_collection_processing/Pa_compendium_02.22.2014.pcl
python Data_collection_processing/zero_one_normalization.py Data_collection_processing/Anr-gene.pcl Train_test_DAs/Anr_normalized.pcl Data_collection_processing/Pa_compendium_02.22.2014.pcl
python Data_collection_processing/zero_one_normalization.py Data_collection_processing/Anr_RNAseq_PAO1_gene.pcl Train_test_DAs/Anr_RNAseq_PAO1_gene_normalized.pcl Data_collection_processing/Pa_compendium_02.22.2014.pcl
python Data_collection_processing/zero_one_normalization.py Data_collection_processing/Anr_RNAseq_J215_gene.pcl Train_test_DAs/Anr_RNAseq_J215_gene_normalized.pcl Data_collection_processing/Pa_compendium_02.22.2014.pcl

#######################################
#Section two: training and testing DAs#
#######################################

#Running the following codes requires installing python packages Theano and docopt
#Instruction for Theano: http://deeplearning.net/software/theano/install.html
#Instruction for docopt: https://pypi.python.org/pypi/docopt

#Train a denoising autoencoder using network size 50, batch size 10, epoch size 500, corruption level 0.1, learning rate 0.01
python Train_test_DAs/SdA_train.py Train_test_DAs/train_set_normalized.pcl 0 50 10 500 0.1 0.01 --seed1 123 --seed2 123

#Test the test set on already trained DA. 
python Train_test_DAs/SdA_test.py Train_test_DAs/Genome-hybs_normalized.pcl 0 Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt 50
python Train_test_DAs/SdA_test.py Train_test_DAs/Anr_normalized.pcl 0 Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt 50
python Train_test_DAs/SdA_test.py Train_test_DAs/Anr_RNAseq_PAO1_gene_normalized.pcl 0 Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt 50
python Train_test_DAs/SdA_test.py Train_test_DAs/Anr_RNAseq_J215_gene_normalized.pcl 0 Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt 50

#Plot the distribution of activity values for each hidden node and the distribution of weight vector for each node
mkdir -p Train_test_DAs/activity_plot
mkdir -p Train_test_DAs/weight_plot
mkdir -p Train_test_DAs/raw_activity_plot
R --no-save < Train_test_DAs/plot_distribution.R

#Get high-weight genes for each nodes
mkdir -p Train_test_DAs/high_weight_gene/2std
python Train_test_DAs/high_weight_gene.py 2 Train_test_DAs/train_set_normalized.pcl Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt Train_test_DAs/high_weight_gene/2std

#Change gene symbol to gene name
mkdir -p Train_test_DAs/high_weight_gene/2std_symbol
python Train_test_DAs/symbol2name_HW_gene.py Train_test_DAs/high_weight_gene/2std Train_test_DAs/high_weight_gene/2std_symbol

#Combine high-weight genes for each node into one file
for x in Train_test_DAs/high_weight_gene/2std_symbol/Node*; do basename -s .txt $x > $x.titled; cat $x >>$x.titled;cut -f 1 $x.titled > $x.cut; done
paste Train_test_DAs/high_weight_gene/2std_symbol/*.cut > Train_test_DAs/high_weight_gene/2std_combined.txt
rm Train_test_DAs/high_weight_gene/2std_symbol/*.cut
rm Train_test_DAs/high_weight_gene/2std_symbol/*.titled

#Plot the relationship between the number of nodes a gene contributing high-weight to and the expression variance of the gene
R --no-save < Train_test_DAs/HWgene_variance.R


##################################
#Section three: positive controls#
##################################

#Download operon gold standard from DOOR database and store them in a file
python Genome_organization/download_operon.py Genome_organization/operon_all.txt 0
#only include operons with more than 3 genes
python Genome_organization/download_operon.py Genome_organization/operon_3.txt 3 
#In case the DOOR online database is down, the operon files are also provided in the repository.

#Build a logistic regression model to predict high-weight genes based on genome distance and cooperonicness
R --no-save < Genome_organization/logit_reg.R #this analysis will take about 1.5 hours, you can continue to the next step while it is running.

#Build GSEA gene sets using operons
python Genome_organization/create_gene_set.py Genome_organization/operon_3.txt Genome_organization/operon_set_3.gmt

#Run GSEA analysis to identify operons associated with each node
mkdir -p Genome_organization/GSEA
python Genome_organization/GSEA_allnode.py Train_test_DAs/train_set_normalized.pcl Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt

#Combine significant GSEA results for every node into one file, 0.05 is FDR cutoff
python Genome_organization/combine_GSEA_allnode.py 0.05

#Calculate operon coverage, plot number of significant operons per node
R --no-save --args Genome_organization/GSEA_FDR_0.05.txt < Genome_organization/GSEA_operon_node_fig.R

#Use the weight matrix to do gene function prediction using 1-nearest-neighbor classifier and compare it with the performance of 1000 random weight matrices
R --no-save < Gene-gene_function/function_prediction_kegg.R #this script has a permutation test, so it will be running for a while


###################################
#Section four: node interpretation#
###################################

#Download GO and KEGG terms from Tribe. Only download terms that have at least 5 genes and at most 100 genes.
python Node_interpretation/get_annotations.py Node_interpretation/pseudomonas_GO_terms.txt GO
python Node_interpretation/get_annotations.py Node_interpretation/pseudomonas_KEGG_terms.txt KEGG

#Plot the top 10 most enriched GO and KEGG pathways for each node.
R --no-save < Node_interpretation/plot_enriched_pathways.R

#Plot the heatmaps of activity values for each dataset
mkdir -p Node_interpretation/heatmaps
R --no-save --args 50 Node_interpretation/heatmaps/ Data_collection_processing/datasets_list_02.22.2014.txt < Node_interpretation/make_heatmaps.R

#Plot heatmaps for a test set
R --no-save --args 50 Node_interpretation/heatmaps/ Train_test_DAs/Genome-hybs_normalized.pcl PAO1_PA14 < Node_interpretation/make_heatmaps_testset.R
R --no-save --args 50 Node_interpretation/heatmaps/ Train_test_DAs/Anr_normalized.pcl Anr-Microarray < Node_interpretation/make_heatmaps_testset.R
R --no-save --args 50 Node_interpretation/heatmaps/ Train_test_DAs/Anr_RNAseq_PAO1_gene_normalized.pcl Anr-RNAseq-PAO1 < Node_interpretation/make_heatmaps_testset.R
R --no-save --args 50 Node_interpretation/heatmaps/ Train_test_DAs/Anr_RNAseq_J215_gene_normalized.pcl Anr-RNAseq-J215 < Node_interpretation/make_heatmaps_testset.R

#Do enrichment analysis on a curated gene set.
python Node_interpretation/find_enriched_node.py Train_test_DAs/high_weight_gene/2std Train_test_DAs/train_set_normalized.pcl Node_interpretation/gene_sets/anr_regulated_cooperonic_genes.txt Node_interpretation/anr_enrichment.txt
