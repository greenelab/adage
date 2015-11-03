##########################
#Plot weight distribution, node activity distribution and node raw activity distribution.
##########################

library(ggplot2)

data = read.table('Train_test_DAs/train_set_normalized.pcl', header=F,skip =1,sep= '\t',stringsAsFactors=F)
gene_id = data[,1]
net_size = 50
weight = read.table('Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt',fill=T,sep='\t',header=F,skip=2,stringsAsFactors=F)
activity = read.table('Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_activity_SdA.txt',header=F, sep= '\t', fill=T,skip=1,row.names=1)
raw_activity = read.table('Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_rawActivity_SdA.txt',header=F, sep= '\t', fill=T,skip=1,row.names=1)

weight_plotpath='Train_test_DAs/weight_plot'
activity_plotpath = 'Train_test_DAs/activity_plot'
raw_activity_plotpath = 'Train_test_DAs/raw_activity_plot'

for (i in 1:net_size){
  node_weight = weight[c(1:length(gene_id)),i]
  node_weight = data.frame(node_weight)
  node_weight$node_weight = as.numeric(as.character(node_weight$node_weight))
  
  #plot distribution of node weight    
  path = file.path(weight_plotpath, paste('Node_', (i),'.eps', sep= ''))
  setEPS()
  postscript(file = path)
  title = paste('Weight of Node', i)
  p = ggplot(node_weight, aes(x= node_weight))
  p = p + geom_density() + labs(title=title,x='Weight')
  print(p)
  dev.off()
  
  #plot node activity of each node
  path = file.path(activity_plotpath, paste('Node_', (i),'.eps', sep= ''))
  setEPS()
  postscript(file = path)
  title = paste('Activity of Node', i)
  p = ggplot(activity, aes(x= activity[,i]))
  p = p + geom_density() + labs(title=title,x='Activity')
  print(p)
  dev.off()
  
  #plot node raw activity of each node
  path = file.path(raw_activity_plotpath, paste('Node_', (i),'.eps', sep= ''))
  setEPS()
  postscript(file = path)
  title = paste('Raw activity of Node', i)
  p = ggplot(raw_activity, aes(x= raw_activity[,i]))
  p = p + geom_density() + labs(title=title,x='Raw activity')
  print(p)
  dev.off()
}

