#Plot the number of significant operons per node and calculate operon coverage

library(ggplot2)
sig_operon = commandArgs(trailingOnly = TRUE)[1]
operon = read.table('Genome_organization/operon_set_3.gmt', sep='\n')
total_operon= nrow(operon) #the total number of operons
data = read.table(sig_operon,header=T,sep='\t')
data$new_node = as.numeric(substr(data$node,5,6))-1 #extract the node number
data_ordered = data[order(data$new_node),] #order by the node number
operon_num = as.data.frame(table(data$new_node)) #count the number of significant operons per node
pdf('Genome_organization/GSEA_FDR_0.05.pdf')
p = ggplot(data_ordered,aes(x=new_node))
p + geom_bar(width = 0.1,fill='blue',colour = 'white',binwidth =1)+ theme_bw()+labs(x='Node',y='Number of significant operons')
dev.off()

coverage = length(table(data$operon))/total_operon #calculate the coverage of all operons
sprintf("Significant operons for all nodes cover %.3f%% of all operons being tested",coverage)
