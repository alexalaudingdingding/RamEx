
#################################################################
# Function:  Ramanome IRCA PCA function  
# Author: Yuehui He
# Last update: 2021-09-24, Yuehui He
#################################################################
# install necessary libraries
library(factoextra)
## Clean R environment
rm(list=ls())
setwd('./')
sourcedir <- Sys.getenv("RamEx")
source(sprintf('%s/Rscript/util_clean.R',sourcedir))


args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
	make_option(c("-i", "--clean_data"),type="character", help="Input_the_raman_clean_data"),
	make_option(c("-m", "--meta_data"),type="character", help="Input_the_meta_data"),
	make_option(c("-o", "--out_dir"), type="character", default='PCA', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$clean_data)) stop('Please input the clean raman data')
if(is.null(opts$meta_data)) stop('Please input the meta data')
matrixfile<-opts$clean_data
metadatafile<-opts$meta_data
plot_local_IRCN<-TRUE
plot_global_IRCN<-TRUE

outpath <- opts$out_dir#"outputpath"
category<-c("group_A", "group_B", "group_C")



#outputpath creation
dir.create(outpath)

options(warn=-1)
#-------------------------------
# Metadata input
#-------------------------------
allmetadata<-read.table(metadatafile,header=T,sep="\t",row.names=1); 
metadata<-data.frame(allmetadata[order(rownames(allmetadata)), ]) 
colnames(metadata)<-colnames(allmetadata)
all_group<-colnames(metadata)
all_group_f<-colnames(metadata)[sapply(metadata,class)=="factor"]; try(metadata_f<-metadata[, all_group_f])
all_group_n<-colnames(metadata)[sapply(metadata,class)!="factor"]; try(metadata_n<-metadata[, all_group_n])

if(length(category)>1){
  metadata$group_c <- as.factor(apply( metadata[ , category ], 1 , paste, collapse = "__" ))
  ori_category<-category
  category<-paste(ori_category, collapse="__")
  names(metadata)[length(metadata)]<-category
  all_group_f[length(all_group_f)+1]<-category
}else{
  ori_category<-category
}
Group<-metadata[, category]
#-------------------------------
# cleandata input
#-------------------------------
mat <- read.table(matrixfile, header = T, row.names = 1, sep="\t"); 
raw_wn<-as.numeric(gsub("[A-Z]", "", colnames(mat)))
wn0<-formatC(raw_wn, digits=0, format="f")

cor_mats_rpvalue<-by(as.matrix(mat),Group, rcorr_df)
v_cor_mats_rpvalue<-lapply(cor_mats_rpvalue, function(x) vectorize_dm_rcorr(x, group=NULL, duplicate=FALSE)) 
v_cor_mats_r_df<-v_cor_mats_rpvalue[[1]][,1:2]
v_cor_mats_pvalue_df<-v_cor_mats_rpvalue[[1]][,1:2]
for(name in names(v_cor_mats_rpvalue)){
  #name<-"000h"
  v_cor_mats_r_df_tem<-data.frame(v_cor_mats_rpvalue[[name]][,"value"])
  v_cor_mats_pvalue_df_tem<-data.frame(v_cor_mats_rpvalue[[name]][,"p_value"])
  colnames(v_cor_mats_r_df_tem)<-colnames(v_cor_mats_pvalue_df_tem)<-name
  v_cor_mats_r_df<-cbind(v_cor_mats_r_df,v_cor_mats_r_df_tem)
  v_cor_mats_pvalue_df<-cbind(v_cor_mats_pvalue_df,v_cor_mats_pvalue_df_tem)
}
v_cor_mats_rpvalue_df<-as.matrix(v_cor_mats_r_df[,3:ncol(v_cor_mats_r_df)])
v_cor_mats_rpvalue_df_neg6<-as.matrix(v_cor_mats_r_df[,3:ncol(v_cor_mats_r_df)])
v_cor_mats_pvalue_df_neg6<-as.matrix(v_cor_mats_pvalue_df[,3:ncol(v_cor_mats_pvalue_df)])
#--------------------------------------------------------
# Generate correlation matrix by one category in metadata
# calculate both correlation and p-value
# vectorize both correlation and p-value
#--------------------------------------------------------
cor_mats_rpvalue<-by(as.matrix(mat),Group, rcorr_df)
v_cor_mats_rpvalue<-lapply(cor_mats_rpvalue, function(x) vectorize_dm_rcorr(x, group=NULL, duplicate=FALSE))


#-------------------------------
#PCA based on MP/MI/MC
#-------------------------------
mean_spectra<-aggregate(mat,list(Group), mean)
rownames(mean_spectra)<-mean_spectra[,1]
mean_spectra<-mean_spectra[,-1]
write.csv(data.frame(mean_spectra),paste(outpath,"/meanspectra_mats.csv",sep=""),row.names=T,quote=F)
mean_spectra_mat<-as.matrix(mean_spectra)
mat<-mean_spectra
tem<-mat[1:2,1:5]
MP.pca <- prcomp(mat, scale = F)
#--------------
# 2. MC PCA (connectedness<-0.6 & p.value<=0.05)
#--------------
mat<-t(v_cor_mats_rpvalue_df_neg6)
tem<-mat[1:2,1:5]
MC.pca <- prcomp(mat, scale = F)
#-------------
# 3. MI PCA (connectedness all & p.value<=0.05)
#--------------
mat<-t(v_cor_mats_rpvalue_df)
tem<-mat[1:2,1:5]
tem<-mat[1:2,ncol(mat)]
tem<-mat[1:2,(ncol(mat)-2):ncol(mat)]
MI.pca <- prcomp(mat, scale = F)

output_PCA <- paste(outpath,"/Results_PCA_MPMIMC_20201117/",sep="")
dir.create(output_PCA)
# 1. MP PCA 
A<-fviz_eig(MP.pca, addlabels = TRUE, ylim = c(0, 50))
pdf(paste(output_PCA,"/PCA_MP.pdf",sep=""),height = 8,width = 8)
fviz_pca_ind(MP.pca,
             xlab="PC1 (41.84%)",
             ylab="PC2 (28.55%)",
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             repel = TRUE     # Avoid text overlapping
             )
dev.off()
#--------------
# 2. MC PCA 
#                              "Cluster_6_dist_cor_mats_rpvalue_neg6.csv",sep = ""))
A<-fviz_eig(MC.pca, addlabels = TRUE, ylim = c(0, 50))
pdf(paste(output_PCA,"PCA_MC.pdf",sep=""),height = 8,width = 8)
fviz_pca_ind(MC.pca,
             xlab="PC1 (22.07%)",
             ylab="PC2 (12.26%)",
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             repel = TRUE     # Avoid text overlapping
             )
dev.off()
#--------------
# 3. MI PCA 
#                              "Cluster_6_dist_cor_mats_rpvalue0.csv",sep = ""))
A<-fviz_eig(MI.pca, addlabels = TRUE, ylim = c(0, 50))
pdf(paste(output_PCA,"PCA_MI.pdf",sep=""),height = 8,width = 8)
fviz_pca_ind(MI.pca,
             xlab="PC1 (25.38%)",
             ylab="PC2 (15.57%)",
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             repel = TRUE     # Avoid text overlapping
             )
dev.off()
