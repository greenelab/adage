library(ggplot2)

kegg = read.delim('./Node_interpretation/pseudomonas_KEGG_terms.txt',header=F,stringsAsFactors=F)
data = read.delim('./Train_test_DAs/train_set_normalized.pcl',row.names=1,header=T)
all_genes = rownames(data)
weight = read.delim('./Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt',skip=2,fill=T,header=F)
gene_num = nrow(data)
weight = weight[1:gene_num,]
rownames(weight) = all_genes

#get a list of genes that are covered by at least one KEGG term and combine all terms a gene participates in
annotated_genes = c()
kegg_goldstd = matrix(nrow=0,ncol=2)
for (i in 1:nrow(kegg)){
  term = kegg[i,1]
  genes = unlist(strsplit(kegg[i,3],';'))
  for (gene in genes){
    #only include genes that are measured on the affy chip
    if (gene %in% all_genes){   
      if (!(gene %in% annotated_genes)){
        annotated_genes = c(annotated_genes,gene)
        kegg_goldstd = rbind(kegg_goldstd, c(term,gene))}
      else{
        idx = which(kegg_goldstd[,2]==gene)
        kegg_goldstd[idx,1] = paste(kegg_goldstd[idx,1],term,sep=';')}
    }
  }
}

weight_selected = weight[annotated_genes, ] #Only include genes that have been annotated.
weight_dist = as.matrix(dist(weight_selected)) #calculate the distance matrix

correct_count_all = 0
correct_count_overlap = 0
for (i in 1:length(annotated_genes)){
  gene = annotated_genes[i]
  j = order(weight_dist[,i])[2] #order the weight distance and get the index to the gene that have closest distance
  nearest_gene = annotated_genes[j] 
  predicted_pathway = unlist(strsplit(kegg_goldstd[j,1],';'))
  actual_pathway = unlist(strsplit(kegg_goldstd[i,1],';'))
  #all annotations need to be correct
  if (setequal(predicted_pathway,actual_pathway)){
    correct_count_all = correct_count_all + 1
  }
  #as long as there is an overlap between predicted pathways and actual pathways, count as correct
  if (length(intersect(predicted_pathway,actual_pathway))>=1){
    correct_count_overlap = correct_count_overlap + 1
  }
}
accuracy_all_final = correct_count_all / length(annotated_genes)
accuracy_overlap_final = correct_count_overlap / length(annotated_genes)

#################################################
#Remove operon structures in function prediction#
#################################################

# #read the operon file, get the operon organizations
# operon_file = file('./Genome_organization/operon_3.txt',open='r')
# lines = readLines(operon_file)
# operons = c()
# for (i in 1:length(lines)){
#   operon = strsplit(lines[i],'\t')
#   operons = c(operons,operon)
# }
# 
# correct_count_all = 0
# correct_count_overlap = 0
# for (i in 1:length(annotated_genes)){
#   gene = annotated_genes[i]
#   #if the target gene is in an operon
#   if (gene %in% unlist(operons)){
#     its_operon = which(sapply(1:length(operons), function(x) is.element(gene,unlist(operons[x])))==T)#get the index of operon which contains the target gene
#     its_operon_genes = unlist(operons[its_operon])#get the genes in this operon
#     weight_dist_nonOperon = weight_dist[-which(rownames(weight_dist) %in% its_operon_genes),]#remove genes in this operon
#     j = order(weight_dist_nonOperon[,i])[1] #find the index of the gene having the closest distance
#     nearest_gene = rownames(weight_dist_nonOperon)[j]#get the gene name of the nearest gene
#   }
#   #if the target gene is not in an operon
#   else{
#     j = order(weight_dist[,i])[2]
#     nearest_gene = annotated_genes[j]
#   }
#   predicted_pathway = unlist(strsplit(kegg_goldstd[which(kegg_goldstd[,2]==nearest_gene),1],';'))
#   actual_pathway = unlist(strsplit(kegg_goldstd[i,1],';'))
#   #all annotations need to be correct
#   if (setequal(predicted_pathway,actual_pathway)){
#     correct_count_all = correct_count_all + 1
#   }
#   #as long as there is an overlap between predicted pathways and actual pathways, count as correct
#   if (length(intersect(predicted_pathway,actual_pathway))>=1){
#     correct_count_overlap = correct_count_overlap + 1
#   }
# }
# accuracy_all_nonOperon = correct_count_all / length(annotated_genes)
# accuracy_overlap_nonOperon = correct_count_overlap / length(annotated_genes)
# print (accuracy_all_nonOperon)
# print (accuracy_overlap_nonOperon)

################
#Random network#
################

#Do function  prediction using a permuted weight matrix for 1000 times
permute_N = 1000
accuracy_all_list = c()
accuracy_overlap_list = c()

for (n in 1:permute_N){
  weight_permuted = weight[sample(nrow(weight)),]#permute the weight matrix
  rownames(weight_permuted) = all_genes
  weight_permuted_selected = weight_permuted[annotated_genes, ] #Only include genes that have been annotated.
  weight_permuted_dist = as.matrix(dist(weight_permuted_selected)) #calculate the distance matrix
  correct_count_all = 0
  correct_count_overlap = 0
  for (i in 1:length(annotated_genes)){
    gene = annotated_genes[i]
    j = order(weight_permuted_dist[,i])[2] #index to closest neighbor gene
    nearest_gene = annotated_genes[j]
    predicted_pathway = unlist(strsplit(kegg_goldstd[j,1],';'))
    actual_pathway = unlist(strsplit(kegg_goldstd[i,1],';'))
    #all annotations need to be correct
    if (setequal(predicted_pathway,actual_pathway)){
      correct_count_all = correct_count_all + 1
    }
    #as long as there is an overlap between predicted pathways and actual pathways, count as correct
    if (length(intersect(predicted_pathway,actual_pathway))>=1){
      correct_count_overlap = correct_count_overlap + 1
    }
  }
  accuracy_all = correct_count_all / length(annotated_genes)
  accuracy_overlap = correct_count_overlap / length(annotated_genes)
  accuracy_all_list = c(accuracy_all_list, accuracy_all)
  accuracy_overlap_list = c(accuracy_overlap_list, accuracy_overlap)
}
accuracy_all_list = as.data.frame(accuracy_all_list)
accuracy_overlap_list = as.data.frame(accuracy_overlap_list)

#plot the distribution of accuracy when using random weight matrices to do function prediction
plot_file_all = './Gene-gene_network/allCorrect_acc_dist.pdf'
pdf(plot_file_all, height=5, width = 8)
acc_all_frame = data.frame(density=1, accuracy= accuracy_all_final)
p = ggplot(accuracy_all_list,aes(x=accuracy_all_list*100))
p + geom_histogram(binwidth=0.5) + geom_point(data = acc_all_frame, aes(x=accuracy*100, y = density, colour='red'))+labs(x='Accuracy(100%)')+scale_color_manual("",labels=c('Observed'),values=c('red'))+theme_bw()
dev.off()

plot_file_overlap = './Gene-gene_network/overlapCorrect_acc_dist.pdf'
pdf(plot_file_overlap, height=5, width = 8)
acc_overlap_frame = data.frame(density=1, accuracy= accuracy_overlap_final)
p = ggplot(accuracy_overlap_list,aes(x=accuracy_overlap_list*100))
p + geom_histogram(binwidth=0.5) + geom_point(data = acc_overlap_frame, aes(x=accuracy*100, y = density, colour='red'))+labs(x='Accuracy(100%)')+scale_color_manual("",labels=c('Observed'),values=c('red'))+theme_bw()
dev.off()

print (paste('The accuracy when all predicted terms are correct is ',accuracy_all_final*100,'%'))
print (paste('The accuracy when predicted terms and real terms overlap is ',accuracy_overlap_final*100,'%'))
