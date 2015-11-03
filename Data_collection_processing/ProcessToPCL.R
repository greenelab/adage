#The first argument is the directory to process
processDir = commandArgs( trailingOnly = TRUE)[1]
outputFile = commandArgs( trailingOnly = TRUE)[2]

source("http://bioconductor.org/biocLite.R")

#the affy libraries are for working with CEL files
library("affy")
library("affyio")

#files now holds a list of all celfiles in the directory to be processed
files = list.celfiles(processDir)

#ptype is now a vector with the type of each array
ptype = sapply(files, function(f) read.celfile.header(paste(processDir, f, sep="/"))[1])

#ptype levels are the levels of that vector (i.e. the types of arrays present)
#this is important because only arrays of the same type can be processed together
ptype_levels = levels(as.factor(unlist(ptype)))

for (level in ptype_levels) {
    if (level == 'Pae_G1a') {
        #pfiles vector only contains files of this type
        pfiles = paste(processDir, subset(files, ptype==level), sep="/")
        #ReadAffy loads the array data using the custom CDF
        Data = ReadAffy(filenames = pfiles)
        #rma processes the data by multi-array average expression measure
        express = rma(Data)
        #this line writes out the PCL file. 
        write.exprs(express, file=outputFile)
    }
}


