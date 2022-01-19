
#################################################################
# Function:  Ramanome IRCA wave number trimming function  
# Author: Yuehui He, Shi Huang
# Last update: 2021-09-22, Shi Huang
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
	make_option(c("-i", "--clean_data"),type="character", help="Input_the_clean_data"),
	make_option(c("-o", "--out_dir"), type="character", default='wave_number_trimming', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$clean_data)) stop('Please input the clean data')

matrixfile <- opts$clean_data # "approxfun_Alldata.txt"
outpath <- opts$out_dir #"outputpath" 


#outputpath creation
dir.create(outpath)


options(warn=-1)
#-------------------------------
# Spectral data input
#-------------------------------
mat <- read.table(matrixfile, header = T, row.names = 1, sep="\t");
mat <- data.frame(mat)
#-------------------------------
# scale by dividing sum area
#-------------------------------
raw_wn<-as.numeric(gsub("[A-Z]", "", colnames(mat)))
wn0<-raw_wn
#wn0<-formatC(raw_wn, digits=0, format="f")
if(length(unique(wn0)) == length(raw_wn)){
  colnames(mat)<-paste("B", wn0,sep="")
}else{
  wn1<-formatC(raw_wn, digits=1, format="f")
  colnames(mat)<-paste("B", wn1, sep="")
}
write.table(mat,paste(outpath,"/Wnt_result",".txt",sep=""),row.names=T,quote=F,sep="\t")
cat("Finish the wave number trimming function. \n")

