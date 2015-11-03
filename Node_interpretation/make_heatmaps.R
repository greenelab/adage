####################################################################
#Load data
####################################################################
net_size = commandArgs(trailingOnly = TRUE)[1]
out_folder = commandArgs(trailingOnly = TRUE)[2]
dataset_list = commandArgs(trailingOnly = TRUE)[3]
dir.create(out_folder)

library("gplots")
networkFile = './Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt' 
dataFile = './Train_test_DAs/train_set_normalized.pcl' 
data = read.table(dataFile,header=T,row.names=1,sep='\t',check.names=F)
weight = read.table(networkFile,header=F,skip=2,fill=T,sep='\t', colClasses='character')
gene_num = nrow(data)
weight = weight[1:gene_num,]
weight = data.matrix(weight)
dataset_list = read.table(dataset_list,header=F, sep='\t',stringsAsFactors=F )

####################################################################
#Plot one heatmap for each dataset
####################################################################
#Loop over all datasets
for (row in 1:nrow(dataset_list)){
  dataset = dataset_list[row,1]
  sampleList = unlist(strsplit(dataset_list[row,2],';')) #get the sample names 
  #only consider datasets that are in the train set
  if (sampleList[1] %in% colnames(data)){
    chosen = data[,sampleList] #extract gene expression values of this dataset with columns names
    raw_activity = t(chosen) %*% weight #multiply expression values with weight and get raw activity values of this dataset
    min_val = apply(raw_activity,2,min) #get the minimun value of each column
    raw_activity_scaled = scale(raw_activity,center=min_val,scale=F) #scale the values by substracting minimun of each column
    colnames(raw_activity_scaled) = paste('Node',rep(1:as.numeric(net_size)),sep='') 
    height = 2 + length(sampleList)/6
    outputFile = paste(out_folder,dataset,'.pdf',sep='')
    pdf(outputFile, width=15, height=height)
    dddcolors = colorpanel(255, rgb(0/255, 176/255, 240/255), rgb(230/255, 230/255, 230/255), rgb(255/255, 255/255, 0/255))
    heatmap.2(raw_activity_scaled, trace='none',margins=c(3, 10),main=dataset,cexRow=0.6, cexCol=0.6, lhei = c(1,1),col=dddcolors)
    dev.off()
  }
}

# ###################################################################
# #one heatmap for each node in each experiment
# ###################################################################
#loop over each hidden node
for (i in 1:net_size){
  subDir = paste('Node',i,sep='')
  dir.create(file.path(out_folder, subDir)) #creat a folder with name 'Node#'
  raw_activity = t(data) %*% weight[,i] #get the raw activity values of this hidden node across all datasets
  #loop over each dataset
  for (row in 1:nrow(dataset_list)){
    dataset = dataset_list[row,1]
    sampleList = unlist(strsplit(dataset_list[row,2],';')) #get the sample names 
    #only consider datasets that are in the train set
    if (sampleList[1] %in% colnames(data)){      
      chosen = raw_activity[sampleList,] #get the raw activity values of this hidden node in this dataset
      scaled = data.frame(chosen - mean(chosen)) #scale it by subtracting the mean
      colnames(scaled) = paste('Node',i,sep='')
      scaled = cbind(scaled,scaled) #to draw a heatmap, it needs at least two column. So two same columns are combined.
      outputFile = paste(file.path(out_folder, subDir),'/Node',i,'_',dataset,'.pdf',sep='')
      height = 2 + length(sampleList)/6
      pdf(outputFile, width=5, height=height)
      par(cex.main=0.75)
      dddcolors = colorpanel(255, rgb(0/255, 176/255, 240/255), rgb(230/255, 230/255, 230/255), rgb(255/255, 255/255, 0/255))
      heatmap.2(as.matrix(scaled),trace='none',symbreaks=F, symkey=F, Colv=F, Rowv=F,main= paste(dataset,'\n','Node',i,sep=''),cexRow=0.5,labCol = c('','') ,cexCol=0.1,keysize=3,lhei=c(1,1),margin=c(1,10),density.info='none',col=dddcolors)
      dev.off()
    }
  }
}