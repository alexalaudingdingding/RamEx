
#################################################################
# Function:  Ramanome IRCA Generate correlation matrix
# Author: Yuehui He, Shi Huang  
# Last update: 2021-09-22, Yuehui He, Shi Huang
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
	make_option(c("-g", "--global_bands_annotation_file"),type="character", help="Input_global_bands_annotation_file"),
	make_option(c("-l", "--local_bands_annotation_file"),type="character", help="Input_local_bands_annotation_file"),
	make_option(c("-p", "--pos_edge"),type="logical", default=FALSE, help="Input_the_value_of_the_pos_edge[default %default]"),
	make_option(c("-n", "--neg_edge"),type="logical", default=TRUE, help="Input_the_value_of_the_neg_edge[default %default]"),
	make_option(c("-c", "--cor_cutoff"),type="double", default=0.6, help="The_cutoff_value_of_the_correlation[default %default]"),
	make_option(c("-o", "--out_dir"), type="character", default='correlation_matrix', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$clean_data)) stop('Please input the clean raman data')
if(is.null(opts$meta_data)) stop('Please input the meta data')
matrixfile<-opts$clean_data
metadatafile<-opts$meta_data
GlobalBandsAnnfile <- opts$global_bands_annotation_file #"approxfun_Alldata.txt"
LocalBandsAnnfile <- opts$local_bands_annotation_file
cor_cutoff<-opts$cor_cutoff
outpath <- opts$out_dir#"outputpath" 
category<-c("group_A", "group_B", "group_C")
Pos_Edge<-opts$pos_edge
Neg_Edge<-opts$neg_edge
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

#-------------------------------
# cleandata input
#-------------------------------
mat <- read.table(matrixfile, header = T, row.names = 1, sep="\t"); 
raw_wn<-as.numeric(gsub("[A-Z]", "", colnames(mat)))
wn0<-formatC(raw_wn, digits=0, format="f")
#-------------------------------
# Bands annotation file input
#-------------------------------
if(!is.null(GlobalBandsAnnfile)){ global_bands_ann<-read.table(GlobalBandsAnnfile,header=T, sep="\t") 
Wave_num <- data.frame(Wave_num = wn0)
global_bands_ann <- merge(Wave_num, global_bands_ann, by="Wave_num", all.x=TRUE, sort = TRUE)
global_bands_ann[is.na(global_bands_ann)] <- "Unknown"
global_bands_ann$Wave_num<-as.numeric(global_bands_ann$Wave_num)
global_bands_ann<-global_bands_ann[with(global_bands_ann, order(Wave_num)),]

}else{
  stop("Please provide global Bands annotation file!")
}
if(!is.null(LocalBandsAnnfile)){ local_bands_ann<-read.table(LocalBandsAnnfile,header=T, sep="\t") }

#-------------------------------
# If multiple categories specified
#-------------------------------
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
write.csv(Group,paste(outpath,"/group.csv",sep=""),row.names=T,quote=F)

#-------------------------------
# Generate correlation matrix by one category in metadata
#-------------------------------
cor_mats<-by(mat, Group, cor)
if(!is.null(GlobalBandsAnnfile)){
                                 v_cor_mats<-lapply(cor_mats, function(x) vectorize_dm(x, group=global_bands_ann$Group, duplicate=FALSE) )
                                 }else{
                                 v_cor_mats<-lapply(cor_mats, function(x) vectorize_dm(x, group=NULL, duplicate=FALSE)) 
                                 }

g_array<-sapply(names(cor_mats), function(x) Plot_network_graph(cor_mats[[x]], Pos_Edge=Pos_Edge, Neg_Edge=Neg_Edge, Threshold=cor_cutoff, node_size=1, node_color="Cluster_membership", layout="layout_components", outdir=NULL)) # paste(outpath, x, "_global_network_igraph.pdf",sep="")
g<-g_array["g",]
#-------------------------------
# Global network statistics by igraph
#-------------------------------
out_array<-t(g_array[1:16,])
net_stats<-data.frame(matrix(unlist(out_array), nrow=nrow(out_array), byrow=F))
dimnames(net_stats)<-dimnames(out_array)
sink(paste(outpath,"/global_network_stats_by_", category,".xls",sep="")); cat("\t"); write.table(net_stats,sep="\t",quote=FALSE); sink()
#-------------------------------
# Visualization: Global network statistics by ggplot2
#-------------------------------
net_stats_df <-melt(data.frame( Group=rownames(net_stats), net_stats))
p<-ggplot(net_stats_df, aes(x=Group, y=value))+ geom_bar(stat = "identity") + 
   xlab( category ) + ylab("Network-level topological feature") + 
   coord_flip()+
   facet_grid(. ~ variable, scales="free") 
ggsave(filename=paste(outpath, "/global_network_stats_barplot_by_", category,".png", sep=""), plot=p, limitsize=FALSE, width=24, height=1 + nlevels(Group)*0.2, device='png') 

#-------------------------------
# Cluster membership distribution
#-------------------------------
Cluster_membership<-sapply(g, function(x) clusters(x, mode = "strong")$membership)
Cluster_membership_df<-data.frame(Cluster_membership, global_bands_ann)
#Cluster_membership_df<-data.frame(Cluster_membership)
sink(paste(outpath,"/global_network_Cluster_membership_by_", category,".xls",sep="")); cat("\t"); write.table(Cluster_membership_df,sep="\t",quote=FALSE); sink()

Cluster_walktrap_membership<-sapply(g, function(x) cluster_walktrap(x)$membership)
Cluster_walktrap_membership_df<-data.frame(Cluster_walktrap_membership, global_bands_ann)
sink(paste(outpath,"/global_network_Cluster_walktrap_membership_by_", category,".xls",sep="")); cat("\t"); write.table(Cluster_walktrap_membership_df,sep="\t",quote=FALSE); sink()

#-------------------------------
# Degree distribution
#-------------------------------
Degree<-sapply(g, function(x) igraph::degree(x)); 
Degree_c<-data.frame(Wave_num=as.numeric(gsub("[A-Z]", "", rownames(Degree))), Degree)
Degree_df<-merge(global_bands_ann, Degree_c, by=1); 
rownames(Degree_df)<-paste("B", Degree_df$Wave_num, sep="")
sink(paste(outpath,"/global_network_Degree_by_", category,".xls",sep=""));
cat("\t"); write.table(Degree_df,sep="\t",quote=FALSE); sink()
#-------------------------------
    d<-Degree_df[, -(1:3)]
    ann = data.frame(Degree_df[, 2])
    rownames(d)<-rownames(ann)<-paste(Degree_df[, 1], Degree_df[, 3], sep="__")
    ann_colors = list(Group = ann[, 1])
