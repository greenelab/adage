####################################################################
#Load data
####################################################################
net_size = commandArgs(trailingOnly = TRUE)[1]
out_folder = commandArgs(trailingOnly = TRUE)[2]
testFile = commandArgs(trailingOnly = TRUE)[3]
test_file_name = commandArgs(trailingOnly = TRUE)[4]
dir.create(out_folder)

library("gplots")
networkFile = './Train_test_DAs/train_set_normalized_50_batch10_epoch500_corrupt0.1_lr0.01_seed1_123_seed2_123_network_SdA.txt' 
test = read.table(testFile,header=T,row.names=1,sep='\t',check.names=F)
weight = read.table(networkFile,header=F,skip=2,fill=T,sep='\t', colClasses='character')
gene_num = nrow(test)
weight = weight[1:gene_num,]
weight = data.matrix(weight)

####################################################################
#Plot a heatmap for the test set
####################################################################
raw_activity = t(test) %*% weight #calculate the raw activity of test dataset
min_val = apply(raw_activity,2,min)
raw_activity_scaled = scale(raw_activity,center=min_val,scale=F)
colnames(raw_activity_scaled) = paste('Node',rep(1:as.numeric(net_size)),sep='') 
height = 2 + nrow(raw_activity)/6
outputFile = paste(out_folder,test_file_name,'.pdf',sep='')
pdf(outputFile, width=15, height=height)
dddcolors = colorpanel(255, rgb(0/255, 176/255, 240/255), rgb(230/255, 230/255, 230/255), rgb(255/255, 255/255, 0/255))
heatmap.2(raw_activity_scaled, trace='none',margins=c(3, 10),main=test_file_name,cexRow=0.6, cexCol=0.6, lhei = c(1,1),col=dddcolors)
dev.off()

# ###################################################################
# #one-node heatmap for each node in the test set
# ###################################################################
#loop over each hidden node
for (i in 1:net_size){
  subDir = test_file_name
  dir.create(file.path(out_folder, subDir)) #creat a folder with name 'Node#'
  #Plot test dataset
  raw_activity = t(test) %*% weight[,i]
  scaled = data.frame(raw_activity - mean(raw_activity))
  colnames(scaled) = paste('Node',i,sep='')
  scaled = cbind(scaled,scaled)
  outputFile = paste(file.path(out_folder, subDir),'/Node',i,'_',test_file_name,'.pdf',sep='')
  height = 2 + nrow(scaled)/6
  pdf(outputFile, width=5, height=height)
  par(cex.main=0.75)
  dddcolors = colorpanel(255, rgb(0/255, 176/255, 240/255), rgb(230/255, 230/255, 230/255), rgb(255/255, 255/255, 0/255))
  heatmap.2(as.matrix(scaled),trace='none',symbreaks=F, symkey=F, Colv=F, Rowv=F,main= paste(test_file_name,'\n','Node',i,sep=''),cexRow=0.5,labCol = c('','') ,cexCol=0.1,keysize=3,lhei=c(1,1),margin=c(1,10),density.info='none',col=dddcolors)
  dev.off()  
}