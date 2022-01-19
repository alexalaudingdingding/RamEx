
#################################################################
# Function:  Ramanome IRCA clean_data function
# Author: Shi Huang   
# Last update: 2022-01-10, Yuehui He, Shi Huang
#################################################################
# install necessary libraries

## Clean R environment
rm(list=ls())
setwd('./')
sourcedir <- Sys.getenv("RamEx")
source(sprintf('%s/Rscript/util_clean.R',sourcedir))


args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
	make_option(c("-i", "--raman_data"),type="character", help="Input_the_raman_data"),
	make_option(c("-o", "--out_dir"), type="character", default='clean_data', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$raman_data)) stop('Please input a test raman file')

matrixfile <- opts$raman_data # "approxfun_Alldata.txt"
outpath <- opts$out_dir #"outputpath" 


#outputpath creation
dir.create(outpath)


options(warn=-1)
#-------------------------------
# Spectral data input
#-------------------------------
mat <- read.table(matrixfile, header = T, row.names = 1, sep="\t");
#-------------------------------
# scale by dividing sum area
#-------------------------------
scale_sum<-function(data){
  sum<-sum(data)
  nor<-data/sum
}
mat<-apply(mat, 1, scale_sum)
mat<-data.frame(t(mat))
write.table(mat,paste(outpath,"/Nor_data.txt",sep=""),row.names=T,quote=F, sep="\t")
mat<-CleanData(mat)
write.table(mat,paste(outpath,"/Clean_data.txt",sep=""),row.names=T,quote=F, sep="\t")
cat("The number of Raman shifts: ", ncol(mat) , "\n")

