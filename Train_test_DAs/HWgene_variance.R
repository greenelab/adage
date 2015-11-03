networkFile = './Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt'
dataFile = './Data_collection_processing/Pa_compendium_02.22.2014.pcl'
data = read.table(dataFile,header=T,row.names=1,sep='\t',check.names=F,quote='')
weight = read.table(networkFile,header=F,skip=2,fill=T,sep='\t', colClasses='character')
gene_num = nrow(data)
weight = weight[1:gene_num,]
weight = data.matrix(weight)
net_size = 50

#Get high-weight genes for each node and combine them into one gene list
high_weight_genes = c()
for (i in 1:net_size){
  HW_file = paste('./Train_test_DAs/high_weight_gene/2std/Node',i,'.txt',sep='')
  HW_genes = read.table(HW_file, header=T,sep='\t',stringsAsFactors =F)
  high_weight_genes = c(high_weight_genes, HW_genes$gene)
}

gene_var = apply(data, 1, function(x) var(x)) #calculate the expression variance of each gene
gene_var = as.data.frame(gene_var)
count = as.data.frame(table(high_weight_genes)) #count the times of appearance of each gene in the high-weight gene list
rownames(count) = count[,1]
merged = merge(count, gene_var,by=0,all=T)
merged$Freq[is.na(merged$Freq)] = 0 #for genes never appear in the high-weight gene list, its appearance is set to 0
print (cor(merged$Freq,merged$gene_var)) #calculate the correlation between the number of nodes a gene gives high weight to and its expression variance

library(ggplot2)
#plot the relationship between the number of nodes a gene gives high-weight to and its expression variance
plot_file_1 = './Train_test_DAs/NodeNum_GeneVar.pdf'
pdf(plot_file_1, height = 5, width = 8)
p = ggplot(merged,aes(factor(Freq),gene_var))
p + geom_boxplot()+labs( x= 'Number of nodes a gene gives high weight to',y='Variance of gene expression')
dev.off()

