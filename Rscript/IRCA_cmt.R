
#################################################################
# Function:  Ramanome IRCA Correlation matrix transfer function
# Author: Yuehui He, Shi Huang
# Last update: 2021-09-24, Yuehui He, Shi Huang
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
	make_option(c("-o", "--out_dir"), type="character", default='correlation_matrix_transfer', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$clean_data)) stop('Please input the clean raman data')
if(is.null(opts$meta_data)) stop('Please input the meta data')
matrixfile<-opts$clean_data
metadatafile<-opts$meta_data
plot_local_IRCN<-FALSE
plot_global_IRCN<-TRUE

outpath <- opts$out_dir#"outputpath"
category<-c("group_A")
cor_cutoff<-0.6


#outputpath creation
dir.create(outpath)
outpath_png <- paste(outpath,"Global-Local-IRCA_png/",sep="/")
dir.create(outpath_png)

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
#---------------------------------------------------------------------------------
#cor_msts transfer to cor_mats_rpvalue for Plot_global_chordDiagram_rpvalue
#---------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------------------------------
# calculate connectedness P<0.05 (i.e.MeansNegCorr / MeansPosCorr / SumsNegCorr / SumsPosCorr/ Degree_Neg / Degree_Pos)
# cor_mats_rpvalue_new (P<0.05) cor_cutoff<-0
#--------------------------------------------------------------------------------------------------------
outpath_global_stat<-paste(outpath,"/global_stat/",sep="")
dir.create(outpath_global_stat)
cor_cutoff_new<-0
SumsNegCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
SumsPosCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
MeansNegCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
MeansPosCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
CountsNegCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
CountsPosCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
for (group in names(cor_mats_rpvalue_sig005)) {
  #group="000h"
  a1<-cor_mats_rpvalue_sig005[[group]]
  a2<-a1
  a2[a2<cor_cutoff_new]<-0 #postive corr
  a3<-a1
  a3[a3>(-cor_cutoff_new)]<-0 #negative corr
  s2<-colSums(a2)#postive corr
  s3<-colSums(a3)#negative corr
  SumsNegCorr_cutoff0<-cbind(SumsNegCorr_cutoff0,data.frame(s3))
  colnames(SumsNegCorr_cutoff0)[ncol(SumsNegCorr_cutoff0)]<-group
  SumsPosCorr_cutoff0<-cbind(SumsPosCorr_cutoff0,data.frame(s2))
  colnames(SumsPosCorr_cutoff0)[ncol(SumsPosCorr_cutoff0)]<-group
  
  s4<-colMeans(a2)#postive corr
  s5<-colMeans(a3)#negative corr
  MeansNegCorr_cutoff0<-cbind(MeansNegCorr_cutoff0,data.frame(s5))
  colnames(MeansNegCorr_cutoff0)[ncol(MeansNegCorr_cutoff0)]<-group
  MeansPosCorr_cutoff0<-cbind(MeansPosCorr_cutoff0,data.frame(s4))
  colnames(MeansPosCorr_cutoff0)[ncol(MeansPosCorr_cutoff0)]<-group
  
  a2[a2<=0]<-0 
  a2[a2>0]<-1 #count numbers of postive corr 
  a3[a3>=0]<-0 
  a3[a3<0]<-1 #count numbers of negative corr
  count2<-colSums(a2)
  count3<-colSums(a3)
  CountsNegCorr_cutoff0<-cbind(CountsNegCorr_cutoff0,data.frame(count3))
  colnames(CountsNegCorr_cutoff0)[ncol(CountsNegCorr_cutoff0)]<-group
  CountsPosCorr_cutoff0<-cbind(CountsPosCorr_cutoff0,data.frame(count2))
  colnames(CountsPosCorr_cutoff0)[ncol(CountsPosCorr_cutoff0)]<-group
}

SumsNegCorr_cutoff0<-SumsNegCorr_cutoff0[-1]
SumsPosCorr_cutoff0<-SumsPosCorr_cutoff0[-1]
MeansNegCorr_cutoff0<-MeansNegCorr_cutoff0[-1]
MeansPosCorr_cutoff0<-MeansPosCorr_cutoff0[-1]
CountsNegCorr_cutoff0<-CountsNegCorr_cutoff0[-1]
CountsPosCorr_cutoff0<-CountsPosCorr_cutoff0[-1]
write.csv(SumsNegCorr_cutoff0,paste(outpath_global_stat,"rpvalue_SumsNegCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(SumsPosCorr_cutoff0,paste(outpath_global_stat,"rpvalue_SumsPosCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(MeansNegCorr_cutoff0,paste(outpath_global_stat,"rpvalue_MeansNegCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(MeansPosCorr_cutoff0,paste(outpath_global_stat,"rpvalue_MeansPosCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(CountsNegCorr_cutoff0,paste(outpath_global_stat,"rpvalue_Degree_Neg_cutoff_",cor_cutoff_new,"_CountsNegCorr.csv",sep=""))
write.csv(CountsPosCorr_cutoff0,paste(outpath_global_stat,"rpvalue_Degree_Pos_cutoff_",cor_cutoff_new,"_CountsPosCorr.csv",sep=""))
#-------------------------------------------------------------------------------------------------------
# calculate connectedness P<0.05 (i.e.MeansNegCorr / MeansPosCorr / SumsNegCorr / SumsPosCorr/ Degree_Neg / Degree_Pos)
# cor_mats_rpvalue_new (P<0.05)
#--------------------------------------------------------------------------------------------------------
cor_cutoff_new<-0.6
SumsNegCorr<-data.frame(A=rep("A",ncol(mat)))
SumsPosCorr<-data.frame(A=rep("A",ncol(mat)))
MeansNegCorr<-data.frame(A=rep("A",ncol(mat)))
MeansPosCorr<-data.frame(A=rep("A",ncol(mat)))
CountsNegCorr<-data.frame(A=rep("A",ncol(mat)))
CountsPosCorr<-data.frame(A=rep("A",ncol(mat)))
for (group in names(cor_mats_rpvalue_sig005)) {
  #group="000h"
  a1<-cor_mats_rpvalue_sig005[[group]]
  a2<-a1
  a2[a2<cor_cutoff_new]<-0 #postive corr
  a3<-a1
  a3[a3>(-cor_cutoff_new)]<-0 #negative corr
  s2<-colSums(a2)#postive corr
  s3<-colSums(a3)#negative corr
  SumsNegCorr<-cbind(SumsNegCorr,data.frame(s3))
  colnames(SumsNegCorr)[ncol(SumsNegCorr)]<-group
  SumsPosCorr<-cbind(SumsPosCorr,data.frame(s2))
  colnames(SumsPosCorr)[ncol(SumsPosCorr)]<-group
  
  s4<-colMeans(a2)#postive corr
  s5<-colMeans(a3)#negative corr
  MeansNegCorr<-cbind(MeansNegCorr,data.frame(s5))
  colnames(MeansNegCorr)[ncol(MeansNegCorr)]<-group
  MeansPosCorr<-cbind(MeansPosCorr,data.frame(s4))
  colnames(MeansPosCorr)[ncol(MeansPosCorr)]<-group
  
  a2[a2<cor_cutoff_new]<-0 
  a2[a2>=cor_cutoff_new]<-1 #count numbers of postive corr 
  a3[a3>(-cor_cutoff_new)]<-0 
  a3[a3<=(-cor_cutoff_new)]<-1 #count numbers of negative corr
  count2<-colSums(a2)
  count3<-colSums(a3)
  CountsNegCorr<-cbind(CountsNegCorr,data.frame(count3))
  colnames(CountsNegCorr)[ncol(CountsNegCorr)]<-group
  CountsPosCorr<-cbind(CountsPosCorr,data.frame(count2))
  colnames(CountsPosCorr)[ncol(CountsPosCorr)]<-group
}
SumsNegCorr<-SumsNegCorr[-1]
SumsPosCorr<-SumsPosCorr[-1]
MeansNegCorr<-MeansNegCorr[-1]
MeansPosCorr<-MeansPosCorr[-1]
CountsNegCorr<-CountsNegCorr[-1]
CountsPosCorr<-CountsPosCorr[-1]
write.csv(SumsNegCorr,paste(outpath_global_stat,"rpvalue_SumsNegCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(SumsPosCorr,paste(outpath_global_stat,"rpvalue_SumsPosCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(MeansNegCorr,paste(outpath_global_stat,"rpvalue_MeansNegCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(MeansPosCorr,paste(outpath_global_stat,"rpvalue_MeansPosCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(CountsNegCorr,paste(outpath_global_stat,"rpvalue_Degree_Neg_cutoff_",cor_cutoff_new,"_CountsNegCorr.csv",sep=""))
write.csv(CountsPosCorr,paste(outpath_global_stat,"rpvalue_Degree_Pos_cutoff_",cor_cutoff_new,"_CountsPosCorr.csv",sep=""))
#--------------------------------------------------------------------------------------------------------
if(plot_global_IRCN){
  device="png"
  outpath_global_neg6_png_rpvalue <-paste(outpath_png,"/global_neg6_png_rpvalue_degree/",sep="")
  dir.create(outpath_global_neg6_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
  Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
                         x=as.numeric(gsub("[A-Z]", "", row.names(CountsNegCorr))),
                         y=CountsNegCorr[,group_name]/max(CountsNegCorr))
    Cairo(file = paste(outpath_global_neg6_png_rpvalue, "Gobal_chordDiagram1_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=30, height=30, type=device, bg="white",warning=FALSE)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=a_degree,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=NULL,
                                     main_title=group_name,
                                     Pos_Edge=FALSE, 
                                     Neg_Edge=TRUE, 
                                     Threshold=0.6)
    dev.off()  
  }

  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 2.1 # Pos_Edge=T, Neg_Edge=F, Threshold=0.6, degree
  #---------------------------------------------------
  device="png"
  outpath_global_pos6_png_rpvalue <-paste(outpath_png,"/global_pos6_png_rpvalue_degree/",sep="")
  dir.create(outpath_global_pos6_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
  Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
                         x=as.numeric(gsub("[A-Z]", "", row.names(CountsPosCorr))),
                         y=CountsPosCorr[,group_name]/max(CountsPosCorr))
    Cairo(file = paste(outpath_global_pos6_png_rpvalue, "Gobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=30, height=30, type=device, bg="white", warning=FALSE)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=a_degree,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=NULL,
                                     main_title=group_name,
                                     Pos_Edge=TRUE, 
                                     Neg_Edge=FALSE, 
                                     Threshold=0.6)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 3.# Pos_Edge=T, Neg_Edge=T, Threshold=0.6,
  #---------------------------------------------------
  device="png"
  outpath_global_neg6pos6_png_rpvalue <- paste(outpath_png,"/global_neg6pos6_png_rpvalue/",sep="")
  dir.create(outpath_global_neg6pos6_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
  Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    Cairo(file = paste(outpath_global_neg6pos6_png_rpvalue, "Gobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=30, height=30, type=device, bg="white", warning=FALSE)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=NULL,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=NULL,
                                     main_title=group_name,
                                     Pos_Edge=TRUE, 
                                     Neg_Edge=TRUE,
                                     Threshold=0.6)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 4.1 # Pos_Edge=F, Neg_Edge=T, Threshold=0,degree
  #---------------------------------------------------
  device="png"
  outpath_global_neg0_png_rpvalue <- paste(outpath_png,"/global_neg0_png_rpvalue_degree/",sep="")
  dir.create(outpath_global_neg0_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
  Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
                         x=as.numeric(gsub("[A-Z]", "", row.names(CountsNegCorr_cutoff0))),
                         y=CountsNegCorr_cutoff0[,group_name]/max(CountsNegCorr_cutoff0))
    Cairo(file = paste(outpath_global_neg0_png_rpvalue, "Gobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=30, height=30, type=device, bg="white", warning=FALSE)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=a_degree,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=NULL,
                                     main_title=group_name,
                                     Pos_Edge=F, 
                                     Neg_Edge=T,
                                     Threshold=0)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  
  #---------------------------------------------------
  # 5.1 # Pos_Edge=T, Neg_Edge=F, Threshold=0,degree(CountsPosCorr_cutoff0)
  #---------------------------------------------------
  device="png"
  outpath_global_pos0_png_rpvalue <- paste(outpath_png,"/global_pos0_png_rpvalue_degree/",sep="")
  dir.create(outpath_global_pos0_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
  Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
                         x=as.numeric(gsub("[A-Z]", "", row.names(CountsPosCorr_cutoff0))),
                         y=CountsPosCorr_cutoff0[,group_name]/max(CountsPosCorr_cutoff0))
    Cairo(file = paste(outpath_global_pos0_png_rpvalue, "Gobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=30, height=30, type=device, bg="white",warning=FALSE)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=a_degree,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=NULL,
                                     main_title=group_name,
                                     Pos_Edge=T, 
                                     Neg_Edge=F, 
                                     Threshold=0)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 6.# Pos_Edge=T, Neg_Edge=T, Threshold=0,
  #---------------------------------------------------
  device="png"
  outpath_global_neg0pos0_png_rpvalue <- paste(outpath_png,"/global_neg0pos0_png_rpvalue/",sep="")
  dir.create(outpath_global_neg0pos0_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
    Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    Cairo(file = paste(outpath_global_neg0pos0_png_rpvalue, "Gobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=30, height=30, type=device, bg="white",warning=FALSE)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=NULL,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=NULL,
                                     main_title=group_name,
                                     Pos_Edge=T,
                                     Neg_Edge=T,
                                     Threshold=0)
  }
  #-----------------------------------------------------------------------------------------------
}
