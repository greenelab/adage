data_file = './Train_test_DAs/train_set_normalized.pcl'
net_file = './Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt'
max_dist = 10
plot_file = "./Genome_organization/distance_plot.pdf"

#read data file
data = read.table(data_file, header=F,skip =1,sep= '\t',stringsAsFactors=F)
gene_id = data[,1]
gene_id = data.frame(gene_id,stringsAsFactors = F)
gene_num = nrow(gene_id)
net_size = 50

#read weight matrix of denoising autoencoder network
weight = read.table(net_file,fill=T,sep='\t',header=F,skip=2)

#convert weight matrix with real values into a weight matrix that use 1 to indicate a high-weight gene and 0 to indicate a non high-weight gene.
weight_matrix = matrix(data=0,nrow=gene_num,ncol=net_size,dimnames=gene_id)
colnames(weight_matrix) = paste('Node',rep(1:50),sep='')
for (i in 1:net_size){
  node_weight = weight[c(1:gene_num),i]
  node_weight = data.frame(node_weight)
  node_weight$node_weight = as.numeric(as.character(node_weight$node_weight))
  #the cutoff of high-weight genes is +- 2 standard deviations.
  pos_cutoff = mean(node_weight$node_weight) + 2* sd(node_weight$node_weight)
  neg_cutoff = mean(node_weight$node_weight) - 2* sd(node_weight$node_weight)
  high_weight_gene = c(gene_id[node_weight$node_weight >= pos_cutoff,],gene_id[node_weight$node_weight <= neg_cutoff,])
  weight_matrix[high_weight_gene,i] =1
}

#Read operon file
operon_file = file('./Genome_organization/operon_all.txt',open='r')
lines = readLines(operon_file)
operons = c()
for (i in 1:length(lines)){
  operon = strsplit(lines[i],'\t')
  #Only keep operons measured on affy array
  if (sum(unlist(operon) %in% rownames(weight_matrix)) == length(unlist(operon))){
    operons = c(operons,operon)
  }
}

#a function to check whether two genes share an operon
CO_OP = function(target_gene, neighbor_gene){
  test = lapply(operons, function(x) sum(c(target_gene, neighbor_gene) %in% x) == 2)
  if (sum(unlist(test)) == 1){
    return(T)
  }
  else{
    return(F)
  }
}

data_table = matrix(nrow=0,ncol=3) #the table that stores all instances
colnames(data_table) = c('distance','cooperonic','HW')

#this loop will take about 75 mins
#start_time = proc.time()
#loop over each hidden node
for (i in 1:net_size){ 
  #loop over each gene in the node
  for (j in 1:gene_num){
    target_gene = rownames(weight_matrix)[j]
    #if the target gene is a high-weight gene
    if (weight_matrix[j,i] == 1){ 
      #loop over each distance
      for (distance in 1:max_dist){
        
        #if the gene locate at the the lower end of the genome
        if (j <= distance){          
          #left neighbor, the chromosome of P.a. is circular
          neighbor = j-distance+gene_num 
          if(weight_matrix[neighbor,i]==1){HW = 1} #if the neighboring gene is also a high-weight gene, set HW=1
          else{HW= 0} #otherwise, set HW=0
          #check whether the target gene j and its neighboring gene are in the same operon
          neighbor_gene = rownames(weight_matrix)[neighbor]
          co_op = CO_OP(target_gene, neighbor_gene)
          data_table = rbind(data_table, c(distance, co_op, HW))
        
          #right neighbor
          neighbor = j+distance  
          if(weight_matrix[neighbor,i]==1){HW = 1} #if the neighboring gene is also a high-weight gene, set HW=1
          else{HW= 0} #otherwise, set HW=0
          #check whether the target gene j and its neighboring gene are in the same operon
          neighbor_gene = rownames(weight_matrix)[neighbor]
          co_op = CO_OP(target_gene, neighbor_gene)
          data_table = rbind(data_table, c(distance, co_op, HW))  
        }
        
        #if the gene locate at the the higher end of the genome
        else if(j > gene_num-distance){
          #left neighbor
          neighbor = j-distance
          if(weight_matrix[neighbor,i]==1){HW =1}
          else{HW=0}
          neighbor_gene = rownames(weight_matrix)[neighbor]
          co_op = CO_OP(target_gene, neighbor_gene)
          data_table = rbind(data_table, c(distance, co_op, HW))
          
          #right neighbor
          neighbor = j+distance-gene_num
          if(weight_matrix[neighbor,i]==1){HW =1}
          else{HW=0}
          neighbor_gene = rownames(weight_matrix)[neighbor]
          co_op = CO_OP(target_gene, neighbor_gene)
          data_table = rbind(data_table, c(distance, co_op, HW))
        }
        
        #if the gene does not locate in the middle of the genome
        else{
          #left neighbor
          neighbor = j-distance
          if(weight_matrix[neighbor,i]==1){HW = 1}
          else{HW = 0}
          neighbor_gene = rownames(weight_matrix)[neighbor]
          co_op = CO_OP(target_gene, neighbor_gene)
          data_table = rbind(data_table, c(distance, co_op, HW))
          
          #right neighbor
          neighbor = j+distance
          if(weight_matrix[neighbor,i]==1){HW = 1}
          else{HW = 0}
          neighbor_gene = rownames(weight_matrix)[neighbor]
          co_op = CO_OP(target_gene, neighbor_gene)
          data_table = rbind(data_table, c(distance, co_op, HW))
        }
      }
    }
  }
}
#proc.time() - start_time 
write.table(data_table, './Genome_organization/dist_coop_HW.txt',quote=F, sep='\t',row.names=F,col.names=T)#write data_table into a file
data_table = as.data.frame(data_table)
#data_table = read.table('./Genome_organization/dist_coop_HW.txt',sep='\t',header=T) 

data_count = as.data.frame(table(data_table))
data_summary = data_count[which(data_count$HW==0),c('distance','cooperonic')]
data_summary$nonHW = data_count[which(data_count$HW==0),'Freq']
data_summary$HW = data_count[which(data_count$HW==1),'Freq']
data_summary$prop = data_summary$HW/ (data_summary$HW+data_summary$nonHW)
data_summary$distance = as.numeric(data_summary$distance)

#logistic regression
data_table = as.data.frame(data_table)
data_table$HW = factor(data_table$HW)
data_table$cooperonic = factor(data_table$cooperonic)
logit = glm(HW ~ distance + cooperonic + distance*cooperonic, data = data_table, family = "binomial")
summary(logit)
exp(coef(logit))

#predicting new data using the fitted logit model
newdata <- with(data_table, data.frame(distance = rep(seq(from = 1, to = 10, by=1),2), cooperonic = factor(rep(0:1, each = 10))))
newdata_predict = cbind(newdata, predict(logit, newdata, type='link', se = T))
#calculate probability and confident interval
newdata_predict <- within(newdata_predict, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})

background_rate = sum(weight_matrix)/length(weight_matrix)#calculate the background high-weight gene rate

library(ggplot2)
pdf('./Genome_organization/dist_coop_plot.pdf', height=4, width=7)
ggplot(data_summary,aes(x=distance,y = prop,fill=cooperonic))+ geom_bar(stat='identity',position='dodge',alpha=0.9)+scale_fill_brewer('Observed frequency',palette="Set1",labels=c('Noncooperonic','Cooperonic')) + geom_line(data=newdata_predict, aes(x = distance, y = PredictedProb,colour=cooperonic),size=1 ) +geom_ribbon(data=newdata_predict,aes(x= distance,y = PredictedProb,ymin = LL,ymax = UL, fill=cooperonic), alpha = 0.2)+ scale_x_continuous(breaks = seq(1,10,by=1))+ scale_y_continuous(breaks = seq(0.0,0.5,by=0.1),limit = c(0,0.5))+theme_bw()+scale_colour_brewer('Predicted probability',palette='Set1',labels=c('Noncooperonic','Cooperonic'))+labs(x='Distance',y='H|W| Frequency/Probability')+ geom_hline(aes(yintercept=background_rate),linetype='dotted')
dev.off()
