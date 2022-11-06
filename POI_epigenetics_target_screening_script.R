#POI epigenetics

library(GEOquery)
library(AnnoProbe)
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(vioplot)
library(pheatmap)
library(ggVennDiagram)
library(ggplot2)
library(tidyr)
library(car)
library(biomaRt)
library(readr)
library(tibble)

setwd("C:/Users/86186/Desktop/POI_epigenetics/GSE135697")

#get data
GSE135697_exp<-read.csv("GSE135697_exp.csv",sep=",")
rownames(GSE135697_exp)<-GSE135697_exp$probe_id
GSE135697_exp$probe_id<-NULL
head(GSE135697_exp)

#We will fix the mistakes made by the data provider (2022/11/06)
#The data provider confused the labels of the control group and the experimental group
#GI.is control group, GC.is experimental(bPOI) group
GSE135697_exp<-GSE135697_exp[,c(11:20,1:10)]
head(GSE135697_exp)

#check data
max(GSE135697_exp)
min(GSE135697_exp)
all(is.na(GSE135697_exp)==F)
boxplot(GSE135697_exp)
#need standardization
GSE135697_exp_norm<-as.data.frame(normalizeBetweenArrays(GSE135697_exp))
boxplot(GSE135697_exp_norm)

#get chip
GSE135697_gpl<-'GPL16956'
GSE135697_gpl<-idmap(GSE135697_gpl,type="pipe")
head(GSE135697_gpl)

#annotation
GSE135697_exp_norm_note<-GSE135697_exp_norm
GSE135697_exp_norm_note$probe_id<-rownames(GSE135697_exp_norm_note)
GSE135697_exp_norm_note<-merge(GSE135697_exp_norm_note,GSE135697_gpl,by="probe_id")
GSE135697_exp_norm_note$probe_id<-NULL
head(GSE135697_exp_norm_note)

#Average the data corresponding to the same gene
GPL16956_gene<-unique(GSE135697_exp_norm_note$symbol)
head(GPL16956_gene)

GSE135697_exp_avg<-matrix(NA,nrow=length(GPL16956_gene),ncol=length(colnames(GSE135697_exp_norm_note))-1)
GSE135697_exp_avg<-as.data.frame(GSE135697_exp_avg)
colnames(GSE135697_exp_avg)<-colnames(GSE135697_exp_norm_note)[1:20]
head(GSE135697_exp_avg)

for(i in 1:length(GPL16956_gene)){
  gene<-GPL16956_gene[i]
  gene_probe_frame<-GSE135697_exp_norm_note[GSE135697_exp_norm_note$symbol==gene,]
  exp_avg<-colMeans(gene_probe_frame[,-c(length(colnames(gene_probe_frame)))])
  GSE135697_exp_avg[i,]<- exp_avg
  rownames(GSE135697_exp_avg)[i]<-gene
}
#test1
GSE135697_exp_avg["SALRNA1",]
GSE135697_exp_norm_note[GSE135697_exp_norm_note$symbol=="SALRNA1",]
GSE135697_gpl[GSE135697_gpl$symbol=="SALRNA1",]
GSE135697_exp["ASHGA5P040285",]
GSE135697_exp_norm["ASHGA5P040285",]
#"ASHGA5P029592" non-existent
GSE135697_exp["ASHGA5P029592",]
GSE135697_exp_norm["ASHGA5P029592",]
#test2
GSE135697_exp_avg["GACAT2",]
GSE135697_exp_norm_note[GSE135697_exp_norm_note$symbol=="GACAT2",]

#Retain genes annotated as ncRNA
keytypes(org.Hs.eg.db)
GSE135697_gene_map<-bitr(rownames(GSE135697_exp_avg),fromType = 'SYMBOL', toType = "REFSEQ",OrgDb = org.Hs.eg.db)
head(GSE135697_gene_map)
GSE135697_gene_lnc<-GSE135697_gene_map[c(grep("NR",GSE135697_gene_map$REFSEQ),grep("XR",GSE135697_gene_map$REFSEQ)),]
GSE135697_gene_lnc<-unique(GSE135697_gene_lnc$SYMBOL)
head(GSE135697_gene_lnc)

GSE135697_for_limma<-GSE135697_exp_avg[GSE135697_gene_lnc,]
head(GSE135697_for_limma)

#POI for limma
colnames(GSE135697_for_limma)<-c("t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10")
head(GSE135697_for_limma)

GSE135697_group<-c("t","t","t","t","t","t","t","t","t","t","c","c","c","c","c","c","c","c","c","c")
GSE135697_design<-model.matrix(~0+factor(GSE135697_group))
colnames(GSE135697_design)<-levels(factor(GSE135697_group))
rownames(GSE135697_design)<-colnames(GSE135697_for_limma)
GSE135697_design

GSE135697_contrasts_fit<-makeContrasts(t-c,levels = GSE135697_design)
GSE135697_fit<-lmFit(GSE135697_for_limma,GSE135697_design)
GSE135697_fit2<-contrasts.fit(GSE135697_fit,GSE135697_contrasts_fit)
GSE135697_fit2<-eBayes(GSE135697_fit2)
GSE135697_limma_result<-topTable(GSE135697_fit2,coef = 1,n=Inf,adjust="BH")
head(GSE135697_limma_result)

#test
GSE135697_limma_result["SALRNA1",]
GSE135697_for_limma["SALRNA1",]

#Select DEG
limma_GSE135697_sig<-GSE135697_limma_result[(GSE135697_limma_result$P.Val<0.05&(GSE135697_limma_result$logFC>=1|GSE135697_limma_result$logFC<=-1)),]
head(limma_GSE135697_sig)
length(limma_GSE135697_sig$logFC)

#Obtain POI related pathogenic genes from OMIM (https://www.omim.org/entry/300510?search=POI&highlight=poi)
#https://www.omim.org/phenotypicSeries/PS233300
#https://www.omim.org/phenotypicSeries/PS311360
setwd("C:/Users/86186/Desktop/POI_epigenetics")
omim_PS233300<-read.csv("OMIM_PS233300.csv",sep = ",")
omim_PS311360<-read.csv("OMIM_PS311360.csv",sep = ",")
omim_PS233300
omim_PS311360

#get GFF3 data from ensembl
#use mirror : USA (NC) [https] - Duke University, Durham, NC
mart<-useMart('ensembl')
list_db<-listDatasets(mart)
hsp<-biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl",
                        host = "www.ensembl.org")
attributes<-c(
  "ensembl_gene_id",
  "hgnc_symbol",
  "chromosome_name",
  "start_position",
  "end_position"
)
hsp_info<-getBM(attributes = attributes, mart = hsp)
length(hsp_info$ensembl_gene_id)
head(hsp_info)

#Confirm the location of POI genes
omim_POI<-unique(rbind(omim_PS233300,omim_PS311360))
omim_POI_gene<-unique(omim_POI$Gene.Locus)
omim_POI_gene
omim_POI_gene_position<-data.frame()
for(i in 1:length(omim_POI_gene)){
  a<-hsp_info[hsp_info$hgnc_symbol==omim_POI_gene[i],]
  omim_POI_gene_position<-rbind(omim_POI_gene_position,a)
}
omim_POI_gene_position

#Confirm the location of POI DEG lnc
sig_lnc_gene<-rownames(limma_GSE135697_sig)
head(sig_lnc_gene)
sig_lnc_gene_position<-data.frame()
for(i in 1:length(sig_lnc_gene)){
  a<-hsp_info[hsp_info$hgnc_symbol==sig_lnc_gene[i],]
  sig_lnc_gene_position<-rbind(sig_lnc_gene_position,a)
}
head(sig_lnc_gene_position)

#Obtain adjacent lncRNA of POI pathogenic gene
#Take 1Mb as threshold
clum<-data.frame()
for(i in 1:nrow(omim_POI_gene_position)){
  a<-omim_POI_gene_position[i,]
  b<-sig_lnc_gene_position[sig_lnc_gene_position$chromosome_name==a$chromosome_name,]
  c<-a$start_position-1000000
  d<-a$end_position+1000000
  e<-b[b$start_position>c&b$end_position<d,]
  e<-rbind(e,c(0,0,0,0,0))
  colnames(e)<-colnames(a)
  e$target<-a$hgnc_symbol
  e$target_strat<-a$start_position
  e$target_end<-a$end_position
  e$distance<-abs(((e$target_end+e$target_strat)/2)-((e$start_position+e$end_position)/2))
  clum<-rbind(clum,e)
}
clum<-clum[clum$ensembl_gene_id!=0,]
clum<-clum[-grep("_",clum$chromosome_name),]
cis_bind<-clum
cis_bind

write.csv(cis_bind,"cis_bind.csv")

#Predicting the binding degree of lncRNA and RBP
catRAPID_result<-read.csv("catRAPID_result.csv",sep=",")
catRAPID_result<-separate(catRAPID_result,c(2),into = c("RNA","note"),sep="-")
catRAPID_result$note<-NULL
catRAPID_result$bind<-paste(catRAPID_result$Protein,catRAPID_result$RNA,sep="-")

catRAPID_result_PS<-aggregate(catRAPID_result$Prediction.Score,by=list(type=catRAPID_result$bind),mean)
colnames(catRAPID_result_PS)[2]<-"Prediction.Score"
catRAPID_result_PS

catRAPID_result_PzS<-aggregate(catRAPID_result$Prediction.z.Score,by=list(type=catRAPID_result$bind),mean)
colnames(catRAPID_result_PzS)[2]<-"Prediction.z.Score"
catRAPID_result_PzS

catRAPID_result_mean<-merge(catRAPID_result_PS,catRAPID_result_PzS,by="type")
catRAPID_result_mean$bind<-catRAPID_result_mean$type
catRAPID_result_mean<-separate(catRAPID_result_mean,c(1),into = c("TF","RNA"),sep="-")
catRAPID_result_mean
#test
catRAPID_result[catRAPID_result$bind=="SIRT1-CYP3A5",]
catRAPID_result_mean[catRAPID_result_mean$bind=="SIRT1-CYP3A5",]

#HCP5-YB1,rediction.Score:8.23,Prediction.z.Score:-1.09
#get HCP5 data
HCP5_for_limma<-na.omit(GSE135697_exp[GSE135697_gpl[GSE135697_gpl$symbol=="HCP5",][,1],])
HCP5_for_limma_t<-HCP5_for_limma[,1:10]
HCP5_for_limma_c<-HCP5_for_limma[,11:20]
HCP5_for_limma_t$t_mean<-apply(HCP5_for_limma_t,1,mean)
HCP5_for_limma_c$c_mean<-apply(HCP5_for_limma_c,1,mean)
HCP5_for_limma_t
HCP5_for_limma_c
HCP5_for_limma<-cbind(HCP5_for_limma_t,HCP5_for_limma_c)
HCP5_for_limma$t_c_diff<-HCP5_for_limma$t_mean-HCP5_for_limma$c_mean
HCP5_for_limma
#According to the original method, the amount of differential expression should be -1.7922488(ASHGA5P052651)

#LIU-ZhanAo completed all modeling and bioinformatics parts
#All codes are written by LIU-ZhanAo

#HCP5 limma test
HCP5_for_limma<-na.omit(GSE135697_exp_norm[GSE135697_gpl[GSE135697_gpl$symbol=="HCP5",][,1],])
colnames(HCP5_for_limma)<-c("t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10")
head(HCP5_for_limma)

HCP5_group<-c("t","t","t","t","t","t","t","t","t","t","c","c","c","c","c","c","c","c","c","c")
HCP5_design<-model.matrix(~0+factor(HCP5_group))
colnames(HCP5_design)<-levels(factor(HCP5_group))
rownames(HCP5_design)<-colnames(HCP5_for_limma)
HCP5_design

HCP5_contrasts_fit<-makeContrasts(t-c,levels = HCP5_design)
HCP5_fit<-lmFit(HCP5_for_limma,HCP5_design)
HCP5_fit2<-contrasts.fit(HCP5_fit,HCP5_contrasts_fit)
HCP5_fit2<-eBayes(HCP5_fit2)
HCP5_limma_result<-topTable(HCP5_fit2,coef = 1,n=Inf,adjust="BH")
head(HCP5_limma_result)
#Using our method, the differential expression multiple should be -1.2993616(ASHGA5P052651)
#get location
hsp_info[hsp_info$hgnc_symbol=="HCP5",]
hsp_info[hsp_info$hgnc_symbol=="MSH5",]
#distance
abs(((31463170+31478936)/2)-((31739677+31762676)/2))
#HCP5-YB1-MSH5
#Prediction.Score:8.23
#Prediction.z.Score:-1.09
#logFC -1.2993616
#distance 280123.5

#data merge
head(catRAPID_result_mean)
head(cis_bind)
head(limma_GSE135697_sig)

colnames(cis_bind)[2]<-"RNA"
lnc_RBP_target_bind<-merge(cis_bind,catRAPID_result_mean,by="RNA")
cis_lnc_logFC<-limma_GSE135697_sig[cis_bind$RNA,][,1,drop=F]
cis_lnc_logFC$RNA<-rownames(cis_lnc_logFC)
lnc_RBP_target_bind<-merge(lnc_RBP_target_bind,cis_lnc_logFC,by="RNA")
lnc_RBP_target_bind$abslogFC<-abs(lnc_RBP_target_bind$logFC)
lnc_RBP_target_bind

#use Topsis
POI_for_topsis<-lnc_RBP_target_bind[,c(1,6,9,10,11,15)]
POI_for_topsis$node<-paste(POI_for_topsis$RNA,POI_for_topsis$TF,POI_for_topsis$target,sep="-")
POI_for_topsis<-POI_for_topsis[,c(7,3,5,6)]
rownames(POI_for_topsis)<-POI_for_topsis$node
POI_for_topsis$node<-NULL
POI_for_topsis[nrow(POI_for_topsis)+1,]<-c(280123.5,8.23,1.2993616)
rownames(POI_for_topsis)[45]<-c("HCP5-YB1-MSH5")
POI_for_topsis$distance_Multiply_minus_one<-POI_for_topsis$distance*(-1)
POI_for_topsis<-POI_for_topsis[,c(2,3,4)]
POI_for_topsis

z_value <- function(x){
  x / sqrt(sum(x^2))
}
POI_for_topsis_z <- POI_for_topsis %>% mutate(across(c(1:3), z_value))
POI_for_topsis_z
#Get the best combination
z_max <- POI_for_topsis_z %>% summarise(across(c(1:3), max)) %>% unlist
z_max
#Get the worse combination
z_min <- POI_for_topsis_z %>% summarise(across(c(1:3), min)) %>% unlist
z_min
#Calculate Euclid distance
dist <-function(x, std){
  res <- c()
  for ( i in 1 : nrow(x)) {
    res[i] = sqrt(sum((unlist(x[i,])-std)^2))
  }
  
  return(res)
}
#get D+
du <- dist(POI_for_topsis_z, z_max)
#get D-
dn <- dist(POI_for_topsis_z, z_min)
#get CI
POI_for_topsis_result<-POI_for_topsis_z %>% add_column(du = du, dn = dn) %>% 
  mutate(ci= dn/(du+dn)) %>%
  arrange(-ci)
POI_for_topsis_result
#select ci>HCP5-YB1-MSH5 is reserved
POI_for_topsis_result["HCP5-YB1-MSH5",]$ci
POI_for_topsis_sel<-POI_for_topsis_result[POI_for_topsis_result$ci>=POI_for_topsis_result["HCP5-YB1-MSH5",]$ci,]
POI_for_topsis_sel

#merge
lnc_RBP_target_bind$node<-paste(lnc_RBP_target_bind$RNA,lnc_RBP_target_bind$TF,lnc_RBP_target_bind$target,sep="-")
rownames(lnc_RBP_target_bind)<-lnc_RBP_target_bind$node
lnc_RBP_target_bind_sel<-merge(lnc_RBP_target_bind,POI_for_topsis_sel,by="row.names")
lnc_RBP_target_bind_sel<-lnc_RBP_target_bind_sel[,c(1,4,2,5,6,11,7,8,9,10,12,15,23)]
colnames(lnc_RBP_target_bind_sel)[1]<-"node"
lnc_RBP_target_bind_sel

#'lnc_RBP_target_bind_sel' is the potential pathogenic combination screened by our original method
#We decided to carry out validation experiments on SALRNA1-TF-C14orf39



