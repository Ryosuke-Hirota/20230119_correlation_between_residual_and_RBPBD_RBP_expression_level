# This script is to investigate correlation between miRNA biogenesis efficiency and RBP expression level and draw plot about it
# made 2023/01/19

library(stringr)
library(ggplot2)

# set function for calculating outlier
cal.outlier <-function(x){
  q <-as.numeric(quantile(x))
  iqr <-IQR(x)
  outlier1 <-q[2]-iqr*1.5
  outlier2 <-q[4]+iqr*1.5
  outliers <-append(outlier1,outlier2)
  return(outliers)
}

# set function for calculating minimum value of expression level
find.min <-function(x,y){
  for (i in 1:y) {
    m <-x[x[,i]!=0,i]
    mv <-min(m)
    if(i==1){
      mvs <-mv
    }else{
      mvs <-append(mvs,mv)
    }}
  return(min(mvs))
}

# make new directory
setwd("C:/Rdata")
dir.create("20230119_correlation_between_residual_and_RBPBD_RBP_expression_level")

# import table of gene counts
# this table is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\TCGA_plot\20230117_TCGA_colon_gene_counts"
setwd("C:/Rdata/20230117_TCGA_colon_gene_counts")
gene.counts.table <-read.table("table_of_TCGA_colon_gene_counts.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)

# calculate minimum value of gene expression level
table <-gene.counts.table[,c(-1,-2)]
min.gene <-find.min(table,293)

# import list of RBPDB RBP
# this list is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\refereces"
setwd("C:/Rdata/references")
RBPDB.RBP.list <-read.csv("RBPDB_v1.3.1_proteins_human_2012-11-21.csv",sep=",",header = F,stringsAsFactors = F)
RBPDB.RBP <-unique(RBPDB.RBP.list[RBPDB.RBP.list[,5]!="",5])

# import list of RBP2GO RBP
# this list is located at "https://www.dropbox.com/home/Okamura%20Lab%20share%20folder/Hirota/results_and_matterials/RBP2GO's_RBPs_contained_to_CCLE"
setwd("C:/Rdata/RBP2GO's_RBPs_contained_to_CCLE")
RBP2GO.RBP.list <-read.table("all_RBP2GO_RBPs_which_merged_with_CCLE.txt",sep="\t",header = T,stringsAsFactors = F)

# import CCLE RBP list
# this list is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\CCLE_Data"
setwd("C:/Rdata/CCLE_data")
df <-read.table("RNAseq_of_RBP.txt",sep="\t",header = T,stringsAsFactors = F)
ccle.RBP <-colnames(df)
ccle.RBP <-ccle.RBP[c(-1,-412)]
m1 <-match(ccle.RBP,RBP2GO.RBP.list[,2])
RBP <-RBP2GO.RBP.list[m1,4]
m2 <-match(RBP,gene.counts.table[,2])
RBP <-gene.counts.table[m2,2]
RBP <-sort(as.character(unique(na.omit(RBP))))

# list file of residual
# these files are located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\TCGA_plot\20230116_TCGA_colon_transcriptome_bam_calculate_residual\table"
residual.file <-list.files(path="C:/Rdata/20230116_TCGA_colon_transcriptome_bam_calculate_residual/table",pattern = ".txt")

# make table about residual file
miRNA.comb <-gsub("(table_of_residual_about_|.txt)","",residual.file)
miRNA.comb.file.table <-as.data.frame(str_split(miRNA.comb,pattern = "_vs_",simplify = T))
miRNA.comb.file.table[,3] <-residual.file


# import correspondence table of TCGA file
# this table is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\TCGA_plot\20230117_TCGA_colon_gene_counts"
setwd("C:/Rdata/20230116_TCGA_colon_make_correspond_gene_counts_files")
cor.table <-read.table("TCGA_colon_gene_counts_correspondence_table.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)

# import list of CCLE pri-miRNA
# this list is located at "https://github.com/Ryosuke-Hirota/20230110_TCGA_transcript_that_intersect_with_miRNA"
setwd("C:/Rdata/20230110_TCGA_transcript_that_intersect_with_miRNA")
mir.list <-read.table("list_of_TCGA_pri-miRNA_that_exsit_in_CCLE.txt",sep="\t",header = T,stringsAsFactors = F)

# make empty summary
sm <-as.data.frame(matrix(nrow = 298*402,ncol = 7))
colnames(sm) <-c("miRNA","transcript","RBP","r","p.value","number_of_sample","number_of_sample_without_expression")

# calculate correlation coefficient and draw plot
for (i in 1:nrow(mir.list)) {
  for (k in 1:length(RBP)) {
  
  if(k==1){
   setwd("C:/Rdata/20230119_correlation_between_residual_and_RBPBD_RBP_expression_level")
   combination <-paste0(mir.list[i,1],"_",mir.list[i,2])
   dir.create(combination)
   }
  
  # import a file about residual between expression levels of miRNA and transcript
  r.file <-miRNA.comb.file.table[miRNA.comb.file.table[,1]==mir.list[i,1]&miRNA.comb.file.table[,2]==mir.list[i,2],3]
  setwd("C:/Rdata/20230116_TCGA_colon_transcriptome_bam_calculate_residual/table")
  residual.df <-read.table(r.file,sep="\t",header = T,stringsAsFactors = F)
  
  # extract expression level of a certain gene from summarized table  
  gene.counts.df <-gene.counts.table[gene.counts.table[,2]==RBP[k],]
  gene.counts.df <-as.data.frame(t(gene.counts.df))
  gene.counts.df[,2] <-rownames(gene.counts.df)
  gene.counts.df <-gene.counts.df[c(-1:-2),]
  rownames(gene.counts.df) <-NULL
  colnames(gene.counts.df) <-c(RBP[k],"filename")
  
  # dataframe
  f <-match(gene.counts.df[,2],cor.table[,2])
  gene.counts.df[,3] <-cor.table[f,4]
  colnames(gene.counts.df)[3] <-"transcriptome_bam"
  
  # merge residual dataframe and gene expression dataframe
  merged.df <-merge(residual.df,gene.counts.df,by="transcriptome_bam")
  merged.df <-merged.df[,c(2,4,1,3,5)]
  
  # count number of sample without RBP expression level
  zero.exp <-nrow(merged.df[merged.df[,2]==0,])
  
  # if a sample hasn't RBP expression level, substitute minmum value 
  merged.df[,2] <-ifelse(merged.df[,2]==0,min.gene,merged.df[,2])

  # normalize RBP expression
  merged.df[,2] <-as.numeric(merged.df[,2])
  merged.df[,2] <-log2(merged.df[,2])
  
  # calculate outliers
  r.outlier <-cal.outlier(merged.df[,1])
  p.outlier <-cal.outlier(merged.df[,2])
  
  # remove outliers
  merged.df <-merged.df[merged.df[,1]>r.outlier[1]&merged.df[,1]<r.outlier[2],]
  merged.df <-merged.df[merged.df[,2]>p.outlier[1]&merged.df[,2]<p.outlier[2],]
  
  # calculate correlation coefficient
  r <-try(cor.test(merged.df[,1],merged.df[,2],method = "spearman"),silent = T)
  
  if(class(r)!="try-error"){
  # draw plot
  p <-ggplot(data = merged.df,aes(x=merged.df[,1],y=merged.df[,2]))+
      geom_point(color="blue")+
      geom_smooth(data=merged.df,mapping = aes(x=merged.df[,1],y=merged.df[,2]),method="lm",formula='y~x',se=FALSE,colour="black",size=0.5)+
      labs(title=paste0("R =",signif(r$estimate,3),", p = ",signif(r$p.value,3),", n = ",nrow(merged.df)),x=paste0(mir.list[i,1],"_",mir.list[i,2]),
           y=RBP[k])+ 
      theme_bw()+
      theme(legend.background = element_rect(fill = "white", colour = "black"))
    
  setwd(paste0("C:/Rdata/20230119_correlation_between_residual_and_RBPBD_RBP_expression_level/",combination))
  ggsave(filename=paste0("plot_of_correlation_between_",mir.list[i,1],"_",mir.list[i,2],"_and_",RBP[k],".pdf"),plot = p)
  
  # write summary
  sm[402*(i-1)+k,1] <-mir.list[i,1]
  sm[402*(i-1)+k,2] <-mir.list[i,2]
  sm[402*(i-1)+k,3] <-RBP[k]
  sm[402*(i-1)+k,4] <-r$estimate
  sm[402*(i-1)+k,5] <-r$p.value
  sm[402*(i-1)+k,6] <-nrow(merged.df)
  sm[402*(i-1)+k,7] <-zero.exp
  }else{
    # write summary
    sm[402*(i-1)+k,1] <-mir.list[i,1]
    sm[402*(i-1)+k,2] <-mir.list[i,2]
    sm[402*(i-1)+k,3] <-RBP[k]
    sm[402*(i-1)+k,4] <-NA
    sm[402*(i-1)+k,5] <-NA
    sm[402*(i-1)+k,6] <-NA
    sm[402*(i-1)+k,7] <-zero.exp
  }
  }}

# output summary
setwd("C:/Rdata/20230119_correlation_between_residual_and_RBPBD_RBP_expression_level")
sm <-sm[order(sm[,4],decreasing = T),]
write.table(sm,"summary_of_TCGA_colon_correlation_between_residual_and_RBP_expression.txt",sep="\t",row.names = F,quote = F)


# do different process by number of sample without RBP expression level
# less than , remove sample
# larger than , substitute minimum value of expression level
#if(zero.exp<100){
#merged.df <-merged.df[merged.df[,2]!=0,]#
#}else{
#merged.df[,2] <-ifelse(merged.df[,2]==0,min.gene,merged.df[,2])
#}




