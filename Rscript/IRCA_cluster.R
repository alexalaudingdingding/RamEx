
#################################################################
# Function:  Ramanome IRCA Correlation matrix transfer cluster function
# Author: Yuehui He 
# Last update: 2021-09-24, Yuehui He
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
	make_option(c("-i", "--clean_data"),type="character", help="Input_the_raman_clean_data"),
	make_option(c("-m", "--meta_data"),type="character", help="Input_the_meta_data"),
	make_option(c("-o", "--out_dir"), type="character", default='cluster', help="outpath file name[default %default]")
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
#--------------------------------------------------------
# Generate correlation matrix by one category in metadata
# calculate both correlation and p-value
# vectorize both correlation and p-value
#--------------------------------------------------------
cor_mats_rpvalue<-by(as.matrix(mat),Group, rcorr_df)
v_cor_mats_rpvalue<-lapply(cor_mats_rpvalue, function(x) vectorize_dm_rcorr(x, group=NULL, duplicate=FALSE))
#-------------------------------
# Graph cluster based on connectedness (correlation and p-value)
#-------------------------------
# 1. cor_mats_rpvalue format change
v_cor_mats_r_df<-v_cor_mats_rpvalue[[1]][,1:2]
v_cor_mats_pvalue_df<-v_cor_mats_rpvalue[[1]][,1:2]
for(name in names(v_cor_mats_rpvalue)){
  v_cor_mats_r_df_tem<-data.frame(v_cor_mats_rpvalue[[name]][,"value"])
  v_cor_mats_pvalue_df_tem<-data.frame(v_cor_mats_rpvalue[[name]][,"p_value"])
  colnames(v_cor_mats_r_df_tem)<-colnames(v_cor_mats_pvalue_df_tem)<-name
  v_cor_mats_r_df<-cbind(v_cor_mats_r_df,v_cor_mats_r_df_tem)
  v_cor_mats_pvalue_df<-cbind(v_cor_mats_pvalue_df,v_cor_mats_pvalue_df_tem)
}
write.csv(file=paste(outpath,'/global_connectness_v_cor_mats_r_df.csv',sep=""),as.data.frame(v_cor_mats_r_df), 
          quote=F, row.names=T)
write.csv(file=paste(outpath,'/global_connectness_v_cor_mats_pvalue_df.csv',sep=""),as.data.frame(v_cor_mats_pvalue_df), 
          quote=F, row.names=T)
# 2. v_cor_mats_df ramanomes distance (connectedness all & p.value no cutoff)
require("proxy")
dist_cor_mats_r<-as.matrix(dist(t(v_cor_mats_r_df[,3:ncol(v_cor_mats_r_df)]), method="Euclidean")) #Euclidean dist
write.csv(file=paste(outpath,'/global_connectness_euclidean_dist_cor_mats_r_hclust.csv',sep=""),as.data.frame(dist_cor_mats_r), 
          quote=F, row.names=T)
mtrx2cols_1col = function(m1,val1){
  lt = lower.tri(m1)  
  res = data.frame(row = row(m1,as.factor = T)[lt],  
                   col = col(m1,as.factor = T)[lt],  
                   val1 = m1[lt]) 
  names(res)[3] = c(val1) 
  return(res)
}
dist_cor_mats_r_df<-mtrx2cols_1col(dist_cor_mats_r,"Euclidean_dist")
write.csv(file=paste(outpath,'/global_connectness_euclidean_dist_cor_r_df_hclust.csv',sep=""),as.data.frame(dist_cor_mats_r_df), 
          quote=F, row.names=T)
# Graph cluster based on dist_cor_mats
CairoPDF(file = paste(outpath, "/global_connectness_euclidean_dist_cor_r_df_hclust.pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");
plot(hclust(as.dist(dist_cor_mats_r), "ward.D"))
dev.off()

# 3.v_cor_mats_rpvalue_df ramanomes distance(connectedness all & p.value<=0.05)
v_cor_mats_rpvalue_df<-as.matrix(v_cor_mats_r_df[,3:ncol(v_cor_mats_r_df)])
dist_cor_mats_rpvalue<-as.matrix(dist(t(v_cor_mats_rpvalue_df), method="Euclidean")) #Euclidean dist
write.csv(file=paste(outpath,'/global_connectness_euclidean_dist_cor_mats_rpvalue_hclust.csv',sep=""),as.data.frame(dist_cor_mats_rpvalue), 
          quote=F, row.names=T)
device="pdf"
CairoPDF(file = paste(outpath, "/global_connectness_euclidean_dist_cor_rpvalue_df_hclust.pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");
plot(hclust(as.dist(dist_cor_mats_rpvalue), "ward.D"))
dev.off()


# 4.v_cor_mats_rpvalue_df ramanomes distance(connectedness<-0.6 & p.value<=0.05)
v_cor_mats_rpvalue_df_neg6<-as.matrix(v_cor_mats_r_df[,3:ncol(v_cor_mats_r_df)])
v_cor_mats_pvalue_df_neg6<-as.matrix(v_cor_mats_pvalue_df[,3:ncol(v_cor_mats_pvalue_df)])
v_cor_mats_rpvalue_df_neg6[which(v_cor_mats_rpvalue_df_neg6>(-0.6))]<-0 #rho>-0.6 
v_cor_mats_rpvalue_df_neg6[which(v_cor_mats_pvalue_df_neg6>(0.05))]<-0 #p.value>0.05


dist_cor_mats_rpvalue_neg6<-as.matrix(dist(t(v_cor_mats_rpvalue_df_neg6), method="Euclidean")) #Euclidean dist
write.csv(file=paste(outpath,'/global_connectness_euclidean_dist_cor_mats_rpvalue_neg6_hclust.csv',sep=""),as.data.frame(dist_cor_mats_rpvalue_neg6), 
          quote=F, row.names=T)
device="pdf"
CairoPDF(file = paste(outpath, "/global_connectness_euclidean_dist_cor_rpvalue_df_neg0.6_hclust.pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");

plot(hclust(as.dist(dist_cor_mats_rpvalue_neg6), "ward.D"))
dev.off()

# 5.v_cor_mats_rpvalue_df ramanomes distance(connectedness<-0.8 & p.value<=0.05)
v_cor_mats_rpvalue_df_neg8<-as.matrix(v_cor_mats_r_df[,3:ncol(v_cor_mats_r_df)])
v_cor_mats_pvalue_df_neg8<-as.matrix(v_cor_mats_pvalue_df[,3:ncol(v_cor_mats_pvalue_df)])
v_cor_mats_rpvalue_df_neg8[which(v_cor_mats_rpvalue_df_neg8>(-0.8))]<-0 #rho>-0.6 
v_cor_mats_rpvalue_df_neg8[which(v_cor_mats_pvalue_df_neg8>(0.05))]<-0 #p.value>0.05
dist_cor_mats_rpvalue_neg8<-as.matrix(dist(t(v_cor_mats_rpvalue_df_neg8), method="Euclidean")) #Euclidean dist
write.csv(file=paste(outpath,'/global_connectness_euclidean_dist_cor_mats_rpvalue_neg8_hclust.csv',sep=""),as.data.frame(dist_cor_mats_rpvalue_neg8), 
          quote=F, row.names=T)
device="pdf"
CairoPDF(file = paste(outpath, "/global_connectness_euclidean_dist_cor_rpvalue_df_neg0.8_hclust.pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");
plot(hclust(as.dist(dist_cor_mats_rpvalue_neg8), "ward.D"))
dev.off()

#-------------------------------
# Graph cluster based on mean spectra
#-------------------------------
mean_spectra<-aggregate(mat,list(Group), mean)
rownames(mean_spectra)<-mean_spectra[,1]
mean_spectra<-mean_spectra[,-1]
write.csv(data.frame(mean_spectra),paste(outpath,"/meanspectra_mats.csv",sep=""),row.names=T,quote=F)
mean_spectra_mat<-as.matrix(mean_spectra)
# 1. jsd distance
jsd_meanspectra<-JSD(mean_spectra_mat)
jsd_meanspectra[is.na(jsd_meanspectra)]<-0
write.csv(data.frame(jsd_meanspectra),paste(outpath,"/jsd_meanspectra.csv",sep=""),row.names=T,quote=F)
CairoPDF(file = paste(outpath, "/global_jsd_meanspectra_hclust.", category, ".pdf" , sep=""),width=3+nlevels(Group)*1.5, height=70, bg="white");
plot(hclust(as.dist(jsd_meanspectra), "ward.D"))
dev.off()
cosine_meanspectra<-as.matrix(proxy::dist(mean_spectra_mat, method="cosine")) #cosine dist
cosine_meanspectra[is.na(cosine_meanspectra)]<-0
write.csv(data.frame(cosine_meanspectra),paste(outpath,"/cosine_meanspectra.csv",sep=""),row.names=T,quote=F)
CairoPDF(file = paste(outpath, "/global_cosine_meanspectra_hclust.", category, ".pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");
plot(hclust(as.dist(cosine_meanspectra), "ward.D"))
dev.off()
# 3. euclidean distance
euclidean_meanspectra<-as.matrix(proxy::dist(mean_spectra_mat, method="Euclidean")) # euclidean dist
euclidean_meanspectra[is.na(euclidean_meanspectra)]<-0
write.csv(data.frame(euclidean_meanspectra),paste(outpath,"/euclidean_meanspectra.csv",sep=""),row.names=T,quote=F)
CairoPDF(file = paste(outpath, "/global_euclidean_meanspectra_hclust.", category, ".pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");
plot(hclust(as.dist(euclidean_meanspectra), "ward.D"))
dev.off()
