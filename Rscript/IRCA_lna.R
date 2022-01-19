
#################################################################
# Function:  Ramanome IRCA Local network by Arcdiagram
# Author: Yuehui He, Shi Huang
# Last update: 2021-09-24, Yuehui He; Shi Huang
#################################################################
# install necessary libraries
library(factoextra)
library(dplyr)
library(tidyverse)
library(circlize)
options(knitr.table.format = "html")
library(viridis)
library(igraph)

library(colormap)
library(hrbrthemes)
library(kableExtra)
library(ggraph)
## Clean R environment
rm(list=ls())
setwd('./')
sourcedir <- Sys.getenv("RamEx")
source(sprintf('%s/Rscript/util_clean.R', sourcedir))


args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
	make_option(c("-i", "--clean_data"),type="character", help="Input_the_raman_clean_data"),
	make_option(c("-m", "--meta_data"),type="character", help="Input_the_meta_data"),
	make_option(c("-l", "--localBandsAnnfile"),type="character", help="Input_the_local_bands_annotation_file"),
	make_option(c("-o", "--out_dir"), type="character", default='local_network_by_arcdiagram', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$clean_data)) stop('Please input the clean raman data')
if(is.null(opts$meta_data)) stop('Please input the meta data')
matrixfile<-opts$clean_data
metadatafile<-opts$meta_data
plot_local_IRCN<-TRUE
LocalBandsAnnfile<- opts$localBandsAnnfile
if(!is.null(LocalBandsAnnfile)){ local_bands_ann<-read.table(LocalBandsAnnfile,header=T, sep="\t") }
Peak_sub<-as.character(local_bands_ann$Wave_num)
bands<-Peak_sub
outpath <- opts$out_dir#"outputpath"
category<-c("group_A", "group_B", "group_C")



#outputpath creation
dir.create(outpath)
outpath_Arc <- paste(outpath,"1.Arcdiagram_Peaks_Neg0/",sep="/")
dir.create(outpath_Arc)
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
cor_mats_rpvalue_sig005<-NULL
for (list_name in names(cor_mats_rpvalue)) {
  #list_name<-"000h"
  cor_mats_r_tem<-cor_mats_rpvalue[[list_name]]$r
  cor_mats_pvalue_tem<-cor_mats_rpvalue[[list_name]]$P
  cor_mats_r_tem[which(cor_mats_pvalue_tem>0.05)]<-0 # p.value>0.05
  #cor_mats_r_tem[which(cor_mats_r_tem>(-0.6))]<-0 # p.value>0.05
  cor_mats_rpvalue_tem<-list(cor_mats_r_tem)
  names(cor_mats_rpvalue_tem)<-list_name
  cor_mats_rpvalue_sig005<-c(cor_mats_rpvalue_sig005,cor_mats_rpvalue_tem)
}
Peaks_fixed_cor_mats_rpvalue <- lapply(cor_mats_rpvalue, function(x) submat <- list(x$r[bands, bands],x$P[bands, bands],x$n[bands, bands]))
names(Peaks_fixed_cor_mats_rpvalue)<-names(cor_mats_rpvalue)
bands_cor_mats_rpvalue<-Peaks_fixed_cor_mats_rpvalue
for (x in names(bands_cor_mats_rpvalue)){
  #--------------------
  # 0.loading dataset
  #--------------------
  #x<-"00h"
  #x<-"CC124__144h"
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  #Corr_mat[Corr_mat>(-0.6)]<-NA # p.value<0.05
  Corr_mat[Corr_mat>(0)]<-NA # p.value<0.05
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  
  #--------------------
  # 1. Make the graph
  # no grouped
  #--------------------
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    dplyr::summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  coauth
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  coauth
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    #scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath_Arc,"/1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 15,width = 20)
  
  #----------------------
  #add group
  # local_bands_ann<-read.table(paste("E:/RWAS/IRCN_CC124_DU_20201022/",
  #                                   "Local_bands_annotation_DU.txt",sep = ""),
  #                             header = T,sep="\t")
  # local_bands_ann$Wave_num<-paste("B",local_bands_ann$Wave_num,sep = "")
  coauth$grp_mannul<-local_bands_ann[which(local_bands_ann$Wave_num==coauth$name),"Group"]
  coauth$grp_mannul<-factor(coauth$grp_mannul)
  coauth<-coauth[order(coauth$grp_mannul),]
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  p23<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], 
                    fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath_Arc,"/2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
         p23,height = 4,width = 2)
}
outpath_Arc <- paste(outpath,"2.Arcdiagram_Peaks_Pos0/",sep="")
dir.create(outpath_Arc)

for (x in names(bands_cor_mats_rpvalue)){
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  Corr_mat[Corr_mat<(0)]<-NA # r<0
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    dplyr::summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  #dim(coauth)
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#FF6666", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    #scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath_Arc,"/1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 15,width = 20)
  
  coauth$grp_mannul<-local_bands_ann[which(local_bands_ann$Wave_num==coauth$name),"Group"]
  coauth$grp_mannul<-factor(coauth$grp_mannul)
  coauth<-coauth[order(coauth$grp_mannul),]
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  p23<-ggraph(mygraph2, layout="linear") + 
    geom_edge_arc(edge_colour="#FF6666", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], 
                    fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath_Arc,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
         p23,height = 4,width = 2)
}
  
################
#3.Neg 6 & p.value<0.05 all peaks
################
#---------------------
# output pathway
#---------------------

outpath_Arc <- paste(outpath,"/3.Arcdiagram_Peaks_Neg6/",sep="")
dir.create(outpath_Arc)

for (x in names(bands_cor_mats_rpvalue)){
  #--------------------
  # 0.loading dataset
  #--------------------
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  Corr_mat[Corr_mat>(-0.6)]<-NA # p.value<0.05
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  
  #--------------------
  # 1. Make the graph
  # no grouped
  #--------------------
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    dplyr::summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath_Arc,"/1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 15,width = 20)
  
  #----------------------
  #add group
  coauth$grp_mannul<-local_bands_ann[which(local_bands_ann$Wave_num==coauth$name),"Group"]
  coauth$grp_mannul<-factor(coauth$grp_mannul)
  coauth<-coauth[order(coauth$grp_mannul),]
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  p23<-ggraph(mygraph2, layout="linear") + 
	geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], 
                    fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath_Arc,"/2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
         p23,height = 4,width = 2)
}  
################
#4.Neg 0.6 & p.value<0.05 DU only
################
#---------------------
# output pathway
#---------------------
outpath_Arc <- paste(outpath,"/4.Arcdiagram_Peaks_Pos6/",sep="")
dir.create(outpath_Arc)

for (x in names(bands_cor_mats_rpvalue)){
  #--------------------
  # 0.loading dataset
  #--------------------
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  Corr_mat[Corr_mat>(-0.6)]<-NA # r<0
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  
  #--------------------
  # 1. Make the graph
  # no grouped
  #--------------------
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    dplyr::summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  #dim(coauth)
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    geom_edge_arc(edge_colour="#FF6666", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath_Arc,"/1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 15,width = 20)
  
  #----------------------
  #add group
  coauth$grp_mannul<-local_bands_ann[which(local_bands_ann$Wave_num==coauth$name),"Group"]
  coauth$grp_mannul<-factor(coauth$grp_mannul)
  coauth<-coauth[order(coauth$grp_mannul),]
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  p23<-ggraph(mygraph2, layout="linear") + 
    geom_edge_arc(edge_colour="#FF6666", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], 
                    fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath_Arc,"/2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),p23,height = 40,width = 20)
}
