source("http://bioconductor.org/biocLite.R")
library("affy")
library("affyio")

test_dir = commandArgs( trailingOnly = TRUE)[1]
test_files = list.celfiles(test_dir)
ptype = sapply(test_files, function(f) read.celfile.header(paste(test_dir, f, sep="/"))[1])
test_pfiles = paste(test_dir, subset(test_files, ptype=='Pae_G1a'), sep="/")

compendium_dir = './Data_collection_processing/data/cels/all-pseudomonas'
comp_files = list.celfiles(compendium_dir)
ptype = sapply(comp_files, function(f) read.celfile.header(paste(compendium_dir, f, sep="/"))[1])
comp_pfiles = paste(compendium_dir, subset(comp_files, ptype=='Pae_G1a'), sep="/")
combined_pfiles = c(test_pfiles, comp_pfiles)

#ReadAffy loads the array data using the custom CDF
Data = ReadAffy(filenames = combined_pfiles)
#rma processes the data by multi-array average expression measure
express = rma(Data)
express_val = exprs(express)
test_set = express_val[,1:length(test_files)]
write.table(test_set, paste(test_dir,'.pcl',sep=''),sep='\t',quote=F,row.names=T,col.names=NA)

