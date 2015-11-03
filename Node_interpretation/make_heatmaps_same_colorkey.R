library("gplots")
####################################################################
#Load data
####################################################################

node = commandArgs(trailingOnly = TRUE)[1]
out_folder = commandArgs(trailingOnly = TRUE)[2]
dataset_list =commandArgs(trailingOnly = TRUE)[3]
testFile = commandArgs(trailingOnly = TRUE)[4]
test_file_name = commandArgs(trailingOnly = TRUE)[5]
key_range = commandArgs(trailingOnly = TRUE)[6]

node = as.numeric(node)
key_min = as.numeric(unlist(strsplit(key_range,','))[1])
key_max = as.numeric(unlist(strsplit(key_range,','))[2])
key_gap = as.numeric(unlist(strsplit(key_range,','))[3])
key_range = seq(key_min,key_max,by=key_gap)
dir.create(out_folder)
networkFile = './Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt' 
dataFile = './Train_test_DAs/train_set_normalized.pcl' 
data = read.table(dataFile,header=T,row.names=1,sep='\t',check.names=F)
weight = read.table(networkFile,header=F,skip=2,fill=T,sep='\t', colClasses='character')
gene_num = nrow(data)
weight = weight[1:gene_num,]
weight = data.matrix(weight)
dataset_list = read.table(dataset_list,header=F, sep='\t',stringsAsFactors=F )
test = read.table(testFile,header=T,row.names=1,sep='\t',check.names=F)

###################################################################
#one heatmap for each node in each experiment
###################################################################
raw_activity = t(data) %*% weight[,node] #get the raw activity values of this hidden node across all datasets
#loop over each dataset
for (row in 1:nrow(dataset_list)){
  dataset = dataset_list[row,1]
  sampleList = unlist(strsplit(dataset_list[row,2],';')) #get the sample names 
  #only consider datasets that are in the train set
  if (sampleList[1] %in% colnames(data)){      
    chosen = raw_activity[sampleList,] #get the raw activity values of this hidden node in this dataset
    scaled = data.frame(chosen - mean(chosen)) #scale it by subtracting the mean
    colnames(scaled) = paste('Node',node,sep='')
    scaled = cbind(scaled,scaled) #to draw a heatmap, it needs at least two column. So two same columns are combined.
    outputFile = paste(out_folder,'/Node',node,'_',dataset,'.pdf',sep='')
    height = 2 + length(sampleList)/6
    pdf(outputFile, width=5, height=height)
    par(cex.main=0.75)
    dddcolors = colorpanel(round((key_max-key_min)/key_gap), rgb(0/255, 176/255, 240/255), rgb(230/255, 230/255, 230/255), rgb(255/255, 255/255, 0/255))
    heatmap.2(as.matrix(scaled),trace='none',breaks = key_range,col=dddcolors, symbreaks=F, symkey=F, Colv=F, Rowv=F,main= paste(dataset,'\n','Node',node,sep=''),cexRow=0.5,labCol = c('','') ,cexCol=0.1,keysize=3,lhei=c(1,1),margin=c(1,10),density.info='none')
    dev.off()
  }
}

#Plot test dataset
raw_activity = t(test) %*% weight[,node]
scaled = data.frame(raw_activity - mean(raw_activity))
colnames(scaled) = paste('Node',node,sep='')
scaled = cbind(scaled,scaled)
outputFile = paste(out_folder,'/Node',node,'_',test_file_name,'.pdf',sep='')
height = 2 + nrow(scaled)/6
pdf(outputFile, width=5, height=height)
par(cex.main=0.75) 
dddcolors = colorpanel(round((key_max-key_min)/key_gap), rgb(0/255, 176/255, 240/255), rgb(230/255, 230/255, 230/255), rgb(255/255, 255/255, 0/255))
heatmap.2(as.matrix(scaled),trace='none',breaks = key_range,col=dddcolors, symbreaks=F, symkey=F, Colv=F, Rowv=F,main= paste(test_file_name,'\n','Node',node,sep=''),cexRow=0.5,labCol = c('','') ,cexCol=0.1,keysize=3,lhei=c(1,1),margin=c(1,10),density.info='none')
dev.off()  