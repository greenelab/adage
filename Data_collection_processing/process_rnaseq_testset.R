#To install TDM package, please refer to https://github.com/greenelab/TDM

library(TDM)
library(data.table)

testFile = commandArgs(trailingOnly = TRUE)[1]
refFile = commandArgs(trailingOnly = TRUE)[2]
outFile = commandArgs(trailingOnly = TRUE)[3]

target_data = read.table(testFile,header=T,sep='\t')
ref_data = read.table(refFile,header=T,sep='\t')
colnames(target_data)[1] = 'gene'
colnames(ref_data)[1] = 'gene'

target_data = data.table(target_data)
setkey(target_data,gene)
ref_data = data.table(ref_data)
setkey(ref_data,gene)

data_tdm=tdm_transform(target_data,ref_data)
write.table(data_tdm,outFile,row.names=F,col.names=T,quote=F,sep='\t')

