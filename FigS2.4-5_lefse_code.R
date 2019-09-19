# XG Manuscript
# LEfSe Analysis
# Last revised by Matt Schnizlein, 6/26/19

setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/experiment_1_20171127_youngmice/mothur_analysis/lefse/lefse_prepost_cef_20190419/")

library(ggplot2)
library(reshape2)
library(dplyr)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n){
  substr(x, 1, n)
}

#####
# Create design file for Pre/Post Cef analysis of both diet groups
# xg.dn14<-subset(genus.2p, group=="xanthan" & timepoint=="-14")
# xg.dn12<-subset(genus.2p, group=="xanthan" & timepoint=="-12")
# xg.dp00<-subset(genus.2p, group=="xanthan" & timepoint=="0")

genus.2p<-read.table(file="../../xanthan_genfrac2p.all_w.meta.txt")

# The contaminated sample was still in the dataframe and needs to be removed...
genus.2p<-subset(genus.2p, sampleID!="17401p01Amis")
genus.2p<-subset(genus.2p, cecum_result!="cecum")
genus.2p<-subset(genus.2p, timepoint!="-3")

genus.2p$timepoint<-factor(genus.2p$timepoint,levels=c(0,2,14,15,17,18,20,21,22,23,24))

# Make the dataframe smaller so we only have the data we need to plot the families
genus.2p<-genus.2p[,c(1,4,6,235:251)]
rownames(genus.2p)<-genus.2p$sampleID

genus.2p$sampleID<-as.character(genus.2p$sampleID)
genus.2p$timepoint<-substrRight(genus.2p$sampleID,3)
xg.diet<-genus.2p[genus.2p$group %in% c("xanthan","xanthan_cdiff") & genus.2p$timepoint %in% c("n14","n12"),]
xg.cef<-genus.2p[genus.2p$group %in% c("xanthan","xanthan_cdiff") & genus.2p$timepoint %in% c("n12","p00"),]

std.diet<-genus.2p[genus.2p$group %in% ("standard_cdiff") & genus.2p$timepoint %in% c("n14","n12"),]
std.cef<-genus.2p[genus.2p$group %in% ("standard_cdiff") & genus.2p$timepoint %in% c("n12","p00"),]

colnames(xg.diet)<-c("Group", "value0")
colnames(xg.cef)<-c("Group", "value0")
colnames(std.diet)<-c("Group", "value0")
colnames(std.cef)<-c("Group", "value0")

# write.table(xg.diet[,1:2], file="lefse/lefse_prepost_cef_20190419/xg_diet_lefse.design",sep="\t",row.names = FALSE,quote=FALSE)
# write.table(xg.cef[,1:2], file="lefse/lefse_prepost_cef_20190419/xg_cef_lefse.design",sep="\t",row.names = FALSE,quote=FALSE)
# write.table(std.diet[,1:2], file="lefse/lefse_prepost_cef_20190419/std_diet_lefse.design",sep="\t",row.names = FALSE,quote=FALSE)
# write.table(std.cef[,1:2], file="lefse/lefse_prepost_cef_20190419/std_cef_lefse.design",sep="\t",row.names = FALSE,quote=FALSE)

# Create the subset the shared (OTU) file based on what conditions are in the design file
shared.lefse<-read.table(file="xanthan_name.final.shared",header=TRUE,sep = "\t")

shared.xg.diet<-shared.lefse[shared.lefse$Group %in% xg.diet$Group,]
shared.xg.cef<-shared.lefse[shared.lefse$Group %in% xg.cef$Group,]
shared.std.diet<-shared.lefse[shared.lefse$Group %in% std.diet$Group,]
shared.std.cef<-shared.lefse[shared.lefse$Group %in% std.cef$Group,]

# write.table(shared.xg.diet, file = "lefse/lefse_prepost_cef_20190419/xg_diet.shared",sep="\t",row.names = FALSE,quote=FALSE)
# write.table(shared.xg.cef, file = "lefse/lefse_prepost_cef_20190419/xg_cef.shared",sep="\t",row.names = FALSE,quote=FALSE)
# write.table(shared.std.diet, file = "lefse/lefse_prepost_cef_20190419/std_diet.shared",sep="\t",row.names = FALSE,quote=FALSE)
# write.table(shared.std.cef, file = "lefse/lefse_prepost_cef_20190419/std_cef.shared",sep="\t",row.names = FALSE,quote=FALSE)

# Put shared and design files into mothur on flux using these commands in mothur:

# Comparing pre/post cef in std chow group
## lefse(shared= std_cef.shared, design= std_cef_lefse.design,class=value0)

# Comparing pre/post diet in std chow group (shouldn't be much)
## lefse(shared= std_diet.shared, design= std_diet_lefse.design,class=value0)

# Comparing pre/post cef in xg chow group
## lefse(shared= xg_cef.shared, design= xg_cef_lefse.design,class=value0)

# Comparing pre/post diet in xg chow group
## lefse(shared= xg_diet.shared, design= xg_diet_lefse.design,class=value0)

# After analyzing with mothur, retrieved the files and brought them back into R to analyze which OTUs were significantly correlated with each group

#####
# Analyzing the mothur output files

setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/experiment_1_20171127_youngmice/mothur_analysis/")

std_cef_lefse<-read.table(file="lefse/lefse_prepost_cef_20190419/std_cef.0.03.lefse_summary",fill=NA,header=TRUE)
std_diet_lefse<-read.table(file="lefse/lefse_prepost_cef_20190419/std_diet.0.03.lefse_summary",fill=NA,header=TRUE)
xg_cef_lefse<-read.table(file="lefse/lefse_prepost_cef_20190419/xg_cef.0.03.lefse_summary",fill=NA,header=TRUE)
xg_diet_lefse<-read.table(file="lefse/lefse_prepost_cef_20190419/xg_diet.0.03.lefse_summary",fill=NA,header=TRUE)

tax_file<-read.table(file="xanthan.taxonomy.names.txt",header=TRUE)
tax_file<-tax_file[,c(1,5,6)]

std_cef_lefse<-std_cef_lefse[!is.na(std_cef_lefse$LDA),]
std_diet_lefse<-std_diet_lefse[!is.na(std_diet_lefse$LDA),]
xg_cef_lefse<-xg_cef_lefse[!is.na(xg_cef_lefse$LDA),]
xg_diet_lefse<-xg_diet_lefse[!is.na(xg_diet_lefse$LDA),]

std_cef_lefse<-merge(std_cef_lefse, tax_file, by.all="OTU")
std_diet_lefse<-merge(std_diet_lefse, tax_file, by.all="OTU")
xg_cef_lefse<-merge(xg_cef_lefse, tax_file, by.all="OTU")
xg_diet_lefse<-merge(xg_diet_lefse, tax_file, by.all="OTU")

std_cef_lefse$LDA.abs<-ifelse(std_cef_lefse$Class=="n12", -std_cef_lefse$LDA, std_cef_lefse$LDA)
xg_cef_lefse$LDA.abs<-ifelse(xg_cef_lefse$Class=="n12", -std_cef_lefse$LDA, std_cef_lefse$LDA)
xg_diet_lefse$LDA.abs<-ifelse(xg_diet_lefse$Class=="n14", -std_cef_lefse$LDA, std_cef_lefse$LDA)

std_cef_lefse$family2<-std_cef_lefse$family
std_cef_lefse$family<-gsub("Bacteroidales_S24-7_group_ge","Bacteroidales_unclassified",std_cef_lefse$family)
std_cef_lefse$family<-gsub("Ruminococcaceae_UCG-008","Ruminococcaceae_unclassified",std_cef_lefse$family)
std_cef_lefse$family<-gsub("Lachnospiraceae_FCS020_group","Lachnospiraceae_unclassified",std_cef_lefse$family)
std_cef_lefse$family<-gsub("uncultured","Other",std_cef_lefse$family)
std_cef_lefse$family<-gsub("Bacteria_unclassified","Other",std_cef_lefse$family)
std_cef_lefse$family<-gsub("Lachnospiraceae_ge","Lachnospiraceae_unclassified",std_cef_lefse$family)
std_cef_lefse$family<-gsub("Lachnospiraceae_UCG-006","Lachnospiraceae_unclassified",std_cef_lefse$family)
std_cef_lefse$family<-gsub("Family_XIII_ge","Other",std_cef_lefse$family)
std_cef_lefse$family<-gsub("Family_XIII_unclassified","Other",std_cef_lefse$family)
std_cef_lefse$family<-gsub("Ruminiclostridium_6","Ruminiclostridium",std_cef_lefse$family)
std_cef_lefse$family<-gsub("Ruminococcaceae_UCG-014","Ruminococcaceae_unclassified",std_cef_lefse$family)
std_cef_lefse$family<-gsub("Ruminococcaceae_UCG-014","Ruminococcaceae_unclassified",std_cef_lefse$family)
std_cef_lefse$family<-gsub("Mollicutes_RF9_ge","Mollicutes_unclassified",std_cef_lefse$family)

xg_cef_lefse$family2<-xg_cef_lefse$family
xg_cef_lefse$family<-gsub("Bacteroidales_S24-7_group_ge","Bacteroidales_unclassified",xg_cef_lefse$family)
xg_cef_lefse$family<-gsub("Ruminococcaceae_UCG-008","Ruminococcaceae_unclassified",xg_cef_lefse$family)
xg_cef_lefse$family<-gsub("Lachnospiraceae_FCS020_group","Lachnospiraceae_unclassified",xg_cef_lefse$family)
xg_cef_lefse$family<-gsub("uncultured","Other",xg_cef_lefse$family)
xg_cef_lefse$family<-gsub("Bacteria_unclassified","Other",xg_cef_lefse$family)
xg_cef_lefse$family<-gsub("Lachnospiraceae_ge","Lachnospiraceae_unclassified",xg_cef_lefse$family)
xg_cef_lefse$family<-gsub("Lachnospiraceae_UCG-006","Lachnospiraceae_unclassified",xg_cef_lefse$family)
xg_cef_lefse$family<-gsub("Ruminiclostridium_6","Ruminiclostridium",xg_cef_lefse$family)
xg_cef_lefse$family<-gsub("Lachnospiraceae_NK4A136_group","Lachnospiraceae_unclassified",xg_cef_lefse$family)
xg_cef_lefse$family<-gsub("Tyzzerella_3","Tyzzerella",xg_cef_lefse$family)

xg_diet_lefse$family2<-xg_diet_lefse$family
xg_diet_lefse$family<-gsub("Bacteroidales_S24-7_group_ge","Bacteroidales_unclassified",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("Ruminococcaceae_UCG-008","Ruminococcaceae_unclassified",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("Lachnospiraceae_FCS020_group","Lachnospiraceae_unclassified",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("Ruminiclostridium_6","Ruminiclostridium",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("uncultured","Other",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("Ruminococcaceae_UCG-014","Ruminococcaceae_unclassified",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("Lachnospiraceae_UCG-006","Lachnospiraceae_unclassified",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("Bacteria_unclassified","Other",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("Lachnospiraceae_ge","Lachnospiraceae_unclassified",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("Lachnospiraceae_UCG-006","Lachnospiraceae_unclassified",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("Lachnospiraceae_NK4A136_group","Lachnospiraceae_unclassified",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("Tyzzerella_3","Tyzzerella",xg_diet_lefse$family)
xg_diet_lefse$family<-gsub("Mollicutes_RF9_ge","Mollicutes_unclassified",xg_diet_lefse$family)

std_cef_lefse$otu_tax<-paste0(std_cef_lefse$OTU,std_cef_lefse$family)
xg_cef_lefse$otu_tax<-paste0(xg_cef_lefse$OTU,xg_cef_lefse$family)
xg_diet_lefse$otu_tax<-paste0(xg_diet_lefse$OTU,xg_diet_lefse$family)

std_cef_lefse$phylum<-factor(std_cef_lefse$phylum, levels=c("Bacteroidetes","Firmicutes","Actinobacteria","Proteobacteria","Tenericutes","Verrucomicrobia","Bacteria_unclassified"))
xg_cef_lefse$phylum<-factor(xg_cef_lefse$phylum, levels=c("Bacteroidetes","Firmicutes","Actinobacteria","Proteobacteria","Tenericutes","Verrucomicrobia","Bacteria_unclassified"))
xg_diet_lefse$phylum<-factor(xg_diet_lefse$phylum, levels=c("Bacteroidetes","Firmicutes","Actinobacteria","Proteobacteria","Tenericutes","Verrucomicrobia","Bacteria_unclassified"))
colors.l<-c("#006d2c","#045a8d","#fd8d3c","#fde157","#bd0026","#6a51a3","gray31")
colors.l<-c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc")

# write.csv(std_cef_lefse, file="lefse/lefse_prepost_cef_20190419/std_cef_lefse_results.txt", row.names = FALSE, quote = FALSE)
# write.csv(std_diet_lefse, file="lefse/lefse_prepost_cef_20190419/std_diet_lefse_results.txt", row.names = FALSE, quote = FALSE)
# write.csv(xg_cef_lefse, file="lefse/lefse_prepost_cef_20190419/xg_cef_lefse_results.txt", row.names = FALSE, quote = FALSE)
# write.csv(xg_diet_lefse, file="lefse/lefse_prepost_cef_20190419/xg_diet_lefse_results.txt", row.names = FALSE, quote = FALSE)


limits.l<-c("Bacteroidetes","Firmicutes","Actinobacteria","Proteobacteria","Tenericutes","Verrucomicrobia","Bacteria_unclassified")
labels.l<-c("Bacteroidetes","Firmicutes","Actinobacteria","Proteobacteria","Tenericutes","Verrucomicrobia","Unclassified Bacteria")

ggplot(std_cef_lefse, aes(x=reorder(OTU,LDA.abs), y=LDA.abs,fill=phylum)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_manual(values=colors.l,breaks=limits.l,labels=labels.l)

ggplot(xg_cef_lefse, aes(x=reorder(OTU,LDA.abs), y=LDA.abs,fill=phylum)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_manual(values=colors.l,breaks=limits.l,labels=labels.l)

ggplot(xg_diet_lefse, aes(x=reorder(OTU,LDA.abs), y=LDA.abs,fill=phylum)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_manual(values=colors.l,breaks=limits.l,labels=labels.l)

unique(subset(xg_diet_lefse, Class=="n12")$family)
unique(subset(xg_diet_lefse, Class=="n14")$family)
subset(xg_diet_lefse,Class=="n12")$OTU

unique(subset(std_cef_lefse, Class=="p00")$family)
unique(subset(std_cef_lefse, Class=="n12")$family)
unique(subset(std_cef_lefse, Class=="n12")$phylum)
subset(std_cef_lefse, Class=="p00")$OTU
subset(std_cef_lefse, Class=="n12")$OTU

# 116 total different OTUs
std_cef_lefse[std_cef_lefse$Class=="n12" & std_cef_lefse$phylum=="Firmicutes",]$phylum
# 93 are Firmicutes
std_cef_lefse[std_cef_lefse$Class=="n12" & std_cef_lefse$phylum=="Bacteroidetes",]$phylum
std_cef_lefse[std_cef_lefse$Class=="p00" & std_cef_lefse$phylum=="Firmicutes",]$phylum

unique(subset(xg_cef_lefse, Class=="p00")$family)
unique(subset(xg_cef_lefse, Class=="p00")$phylum)
unique(subset(xg_cef_lefse, Class=="n12")$family)
unique(subset(xg_cef_lefse, Class=="n12")$phylum)

subset(xg_cef_lefse, Class=="n12")$OTU
subset(xg_cef_lefse, Class=="p00")$OTU
# 84 total different OTUs
xg_cef_lefse[xg_cef_lefse$Class=="n12" & xg_cef_lefse$phylum=="Firmicutes",]$phylum
# 64 are Firmicutes
xg_cef_lefse[xg_cef_lefse$Class=="n12" & xg_cef_lefse$phylum=="Bacteroidetes",]$phylum
# 11 Bacteroidetes
xg_cef_lefse[xg_cef_lefse$Class=="p00" & xg_cef_lefse$phylum=="Firmicutes",]$phylum

# comparing the taxa that are reduced in both xg and std chows
xg_cef_lefse.12<-xg_cef_lefse[xg_cef_lefse$Class %in% "n12",]
std_cef_lefse.12<-std_cef_lefse[std_cef_lefse$Class=="n12",]

# 112 OTUs negatively associated with cef treatment in std treatment
# 83 OTUs negatively associated with Cef treatment in xg treatment

xg_cef_lefse.12[xg_cef_lefse.12$OTU %in% std_cef_lefse.12$OTU,]$family
# 48 OTUs are shared between the two groups (42% or 48/83 of the OTUs in the xtd chow group, 57% or 48/112 of OTUs in the std chow group)

xg_cef_lefse.12[xg_cef_lefse.12$OTU %in% std_cef_lefse.12$OTU,]$phylum

xg_cef_lefse[xg_cef_lefse$pValue<=0.01,]$OTU
xg_cef_lefse$OTU

std_cef_lefse[std_cef_lefse$pValue<=0.01,]$OTU
std_cef_lefse$pValue
xg_diet_lefse[xg_diet_lefse$pValue<=0.01,]$OTU

#####
# OTU shared rel abundance
setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/experiment_1_20171127_youngmice/mothur_analysis/")

lefse.preprocess<-function(lefse,shared,outname){
  shared.tax<-merge(lefse,shared,by.all="OTU")
  shared.tax.otu<-shared.tax[,c(1,6:7,3,2,4,5)]
  shared.tax<-shared.tax[,c(1,3,6:7,11:159)]
  shared.tax[,5:153]<-lapply(lapply(shared.tax[,5:153],as.character),as.numeric)
  assign(paste("shared",outname,"otu",sep="."),shared.tax,envir=.GlobalEnv)
  shared.tax<-aggregate(.~phylum+family+Class,data=shared.tax,FUN=sum)
  shared.tax$Class<-as.character(shared.tax$Class)
  assign(paste("shared",outname,"tax",sep="."),shared.tax,envir=.GlobalEnv)
}
lefse.preprocess.time<-function(x,t1,t2,outname,exclude){
  # t1 indicates the timepoint OTUs are associated with, t2 is the comparison timepoint
  shared.tax.t1<-x[x$Class %in% t1,]
  shared.tax.t1.t<-data.frame(t(shared.tax.t1))
  colnames(shared.tax.t1.t)<-lapply(shared.tax.t1.t[2,],as.character)
  shared.tax.t1.t$time<-substrRight(rownames(shared.tax.t1.t),3)
  shared.tax.t1.t$group<-substrLeft(rownames(shared.tax.t1.t),4)
  shared.tax.t1.t$group2<-shared.tax.t1.t$group
  for (i in c("1738","1739")){
    shared.tax.t1.t$group<-gsub(i,"std",shared.tax.t1.t$group)
  }
  for (i in c("1737","1740","1744","1743")){
    shared.tax.t1.t$group<-gsub(i,"xg",shared.tax.t1.t$group)
  }
  shared.t1<-shared.tax.t1.t[shared.tax.t1.t$time %in% c(t1,t2),]
  shared.t1<-subset(shared.t1,group!=exclude)
  shared.t1.melt<-melt(shared.t1, id=c("time","group","group2"))
  shared.t1.melt$value<-as.numeric(as.character(shared.t1.melt$value))
  print(unique(shared.t1.melt$variable))
  assign(paste("shared",outname,t1,"melt",sep="."),shared.t1.melt,envir=.GlobalEnv)
}
insert_minor <- function(major_labs, n_minor){
  labs <- c(sapply(major_labs, function(x) c(x, rep("", n_minor))))
  labs[1:(length(labs)-n_minor)]
}
lefse.50.plot<-function(d,t1,t2,lim1,lab1,g1name,g2name,outname){
  shared.p<-ggplot(d,aes(x=variable,y=value,group=time,color=time))+
    geom_jitter(position=position_jitterdodge(dodge.width = 0.8,jitter.width = 0.2)) +
    scale_y_continuous(breaks= seq(0,50,by=2),labels = function(x) insert_minor(seq(0,50,by=10), 4),limits = c(0,50)) +
    scale_x_discrete(labels=rev(lab1),limits=rev(lim1))+
    stat_summary(aes(group=time),fun.y=mean,fun.ymin=mean,fun.ymax=mean,geom="crossbar",width=0.5,colour="black",position=position_dodge(width=0.9)) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),legend.title=element_text(size=15),legend.text=element_text(size=12)) +ylab("Relative abundance (%)") +xlab("") +
    theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.8), angle = 00) ) +
    coord_flip()+
    scale_color_manual(values=c("#08519c","#FFCB05"),name="Time", limits=c(t1,t2),labels=c(g1name,g2name)) +ggtitle(paste(outname,"LEfSe",sep=" "))
  plot(shared.p)
}

shared.rel<-read.table(file="xanthan_name.final.relabund.shared",sep="\t",header=TRUE)

xg_diet_lefse<-read.csv(file="lefse/lefse_prepost_cef_20190419/xg_diet_lefse_results.txt",header = TRUE)
xg_cef_lefse<-read.csv(file="lefse/lefse_prepost_cef_20190419/xg_cef_lefse_results.txt",header = TRUE)
std_cef_lefse<-read.csv(file="lefse/lefse_prepost_cef_20190419/std_cef_lefse_results.txt",header=TRUE)
std_diet_lefse<-read.csv(file="lefse/lefse_prepost_cef_20190419/std_diet_lefse_results.txt",header=TRUE)

# colnames(shared.rel)<-as.character(shared.rel[1,])
shared.rel.t<-data.frame(t(shared.rel))
colnames(shared.rel.t)<-lapply(shared.rel.t[2,1:149],as.character)
shared.rel.t<-shared.rel.t[-c(1:3),]
shared.rel.t$OTU<-rownames(shared.rel.t)
shared.rel.t[,1:149]<-lapply(lapply(shared.rel.t[,1:149], as.character),as.numeric)
shared.rel.t[,1:149]<-shared.rel.t[,1:149]*100

# shared.tax<-merge(tax_file,shared.rel.t,by.all="OTU")
# shared.tax<-subset(shared.tax, select=-OTU)
# shared.tax[,3:151]<-lapply(shared.tax[,3:151], as.character)
# shared.tax[,3:151]<-lapply(shared.tax[,3:151], as.numeric)
# shared.tax.sum<-aggregate(.~phylum+family,data=shared.tax,FUN=sum)
#####
# XG Diet Lefse Plot
lefse.preprocess(xg_diet_lefse,shared.rel.t,"xgdiet")
lefse.preprocess.time(shared.xgdiet.tax,"n12","n14","xgdiet","std")
lefse.preprocess.time(shared.xgdiet.tax,"n14","n12","xgdiet","std")
# Warning is from the melt function, doesn't affect output

limits.xgdiet.n12<-c("Bacteroidales_unclassified","Caproiciproducens","Fusicatenibacter","Tyzzerella","Lachnospiraceae_unclassified","Ruminococcaceae_unclassified","Verrucomicrobiaceae_unclassified","Other")
limits.xgdiet.n14<-c("Alistipes","Bacteroidales_unclassified","Coprobacillus","Lachnoclostridium","Lactobacillus","Ruminiclostridium","Syntrophococcus","Clostridiales_unclassified","Firmicutes_unclassified","Lachnospiraceae_unclassified","Ruminococcaceae_unclassified","Eggerthella","Coriobacteriaceae_unclassified","Mollicutes_unclassified","Other")

# Bacteroidetes: "Alistipes","Bacteroidales_unclassified"
# Firmicutes: "Caproiciproducens","Coprobacillus","Fusicatenibacter","Lachnoclostridium","Lactobacillus","Ruminiclostridium","Syntrophococcus","Tyzzerella","Clostridiales_unclassified","Firmicutes_unclassified",,"Lachnospiraceae_unclassified","Ruminococcaceae_unclassified",
# Actinobacteriaceae: "Eggerthella","Coriobacteriaceae_unclassified"
# Tenericutes:,"Mollicutes_unclassified"
# Verrucomicrobiaceae: "Verrucomicrobiaceae_unclassified"
# "Other"

## Lefse OTUS in each taxa group
# xgdiet.otus<-shared.tax.xgdiet.otus[,c(1,4,2,3)]
# xgdiet.otus.n14<-xgdiet.otus[xgdiet.otus$Class %in% "n14",]
# xgdiet.otus.n12<-xgdiet.otus[xgdiet.otus$Class %in% "n12",]

labels.xgdiet.n12<-c("Unclassified Bacteroidales (7)","Caproiciproducens (2)","Fusicatenibacter (2)","Tyzzerella (1)","Unclassified Lachnospiraceae (19)","Unclassified Ruminococcaceae (1)","Unclassified Verrucomicrobiaceae (1)","Unclassified Bacteria (2)")
labels.xgdiet.n14<-c("Alistipes (2)","Unclassified Bacteroidales (6)","Coprobacillus (1)","Lachnoclostridium (1)","Lactobacillus (3)","Ruminiclostridium (1)","Syntrophococcus (1)","Unclassified Clostridiales (1)","Unclassified Firmicutes (1)","Unclassified Lachnospiraceae (5)","Unclassified Ruminococcaceae (2)","Eggerthella (1)","Unclassified Coriobacteriaceae (2)","Unclassified Mollicutes(2)","Unclassified Bacteria (3)")

# Plotting
shared.xgdiet.n14.melt$time<-factor(shared.xgdiet.n14.melt$time,levels=c("n12","n14"))
shared.xgdiet.n12.melt$time<-factor(shared.xgdiet.n12.melt$time,levels=c("n12","n14"))

# lefse.50.plot(shared.xgdiet.n12.melt,"n14","n12",limits.xgdiet.n12,labels.xgdiet.n12,"Pre-Diet Change","Post-Diet Change","XG Diet")
lefse.50.plot(shared.xgdiet.n14.melt,"n14","n12",limits.xgdiet.n14,labels.xgdiet.n14,"Pre-Diet Change","Post-Diet Change","XG Diet")

shared.xgdiet.p<-ggplot(shared.xgdiet.n12.melt,aes(x=variable,y=value,group=time,color=time))+
  geom_jitter(position=position_jitterdodge(dodge.width = 0.7,jitter.width = 0.08)) +
  scale_y_continuous(breaks= seq(0,50,by=2),labels = function(x) insert_minor(seq(0,50,by=10), 4),limits = c(0,50)) +
  scale_x_discrete(labels=rev(labels.xgdiet.n12),limits=rev(limits.xgdiet.n12))+
  stat_summary(aes(group=time),fun.y=mean,fun.ymin=mean,fun.ymax=mean,geom="crossbar",width=0.3,colour="black",position=position_dodge(width=0.7)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),legend.title=element_text(size=15),legend.text=element_text(size=12)) +ylab("Relative abundance (%)") +xlab("") +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00) ) +
  coord_flip()+
  scale_color_manual(values=c("#08519c","#FFCB05"),name="Time", limits=c("n14","n12"),labels=c("Pre-Xanthan Gum","Post-Xanthan Gum")) +ggtitle("XG Diet")
plot(shared.xgdiet.p)


######
# XG Cef Lefse Plot
lefse.preprocess(xg_cef_lefse,shared.rel.t,"xgcef")
lefse.preprocess.time(shared.xgcef.tax,"n12","p00","xgcef","std")
lefse.preprocess.time(shared.xgcef.tax,"p00","n12","xgcef","std")

## Lefse OTUS in each taxa group
# Xg Cef
# xgcef.otus<-shared.tax.xgcef.otu[,c(1,4,2,3)]
# xgcef.otus.n12<-xgcef.otus[xgcef.otus$Class %in% "n12",]
# xgcef.otus.p00<-xgcef.otus[xgcef.otus$Class %in% "p00",]

limits.xgcef.n12<-c("Alistipes","Bacteroidales_unclassified","Blautia","Butyrivibrio","Fusicatenibacter","Lachnoclostridium","Lachnospiraceae_unclassified","Lactobacillus","Roseburia","Ruminiclostridium","Ruminococcaceae_unclassified","Ruminococcus","Syntrophococcus","Turicibacter","Tyzzerella","Other","Bifidobacterium","Eggerthella","Other.1")
limits.xgcef.p00<-c("Bacteroidales_unclassified","Lachnospiraceae_unclassified","Verrucomicrobiaceae_unclassified","Verrucomicrobiales_unclassified")

labels.xgcef.n12<-c("Alistipes (2)","Unclassified Bacteroidales (9)","Blautia (1)","Butyrivibrio (1)","Fusicatenibacter (2)","Lachnoclostridium (1)","Unclassified Lachnospiraceae(41)","Lactobacillus (2)","Roseburia (1)","Ruminiclostridium (1)","Unclassified Ruminococcaceae (4)","Ruminococcus (2)","Syntrophococcus (2)","Turicibacter (1)","Tyzzerella (2)","Unclassified Firmicutes (5)","Bifidobacterium (1)","Eggerthella (1)","Unclassified Bacteria (1)")
labels.xgcef.p00<-c("Unclassified Bacteroidales (1)","Unclassified Lachnospiraceae (1)","Unclassified Verrucomicrobiaceae (1)","Unclassified Verrucomicrobiales(1)")

# Bacteroides: "Alistipes","Bacteroidales_unclassified"
# Firmicutes: "Blautia","Butyrivibrio","Fusicatenibacter","Lachnoclostridium","Lachnospiraceae_unclassified","Lactobacillus","Roseburia","Ruminiclostridium",,"Ruminococcaceae_unclassified","Ruminococcus","Syntrophococcus","Turicibacter","Tyzzerella","Other"
# Actinobacteria: "Bifidobacterium","Eggerthella",
# "Other.1"

shared.xgcef.n12.melt$time<-factor(shared.xgcef.n12.melt$time,levels=c("p00","n12"))
shared.xgcef.p00.melt$time<-factor(shared.xgcef.p00.melt$time,levels=c("p00","n12"))

lefse.50.plot(shared.xgcef.n12.melt,"n12","p00",limits.xgcef.n12,labels.xgcef.n12,"Pre-Cefoperazone","Post-Cefoperazone","XG Cef")
# lefse.50.plot(shared.xgcef.p00.melt,"n12","p00",limits.xgcef.p00,labels.xgcef.p00,"Pre-Cefoperazone","Post-Cefoperazone","XG Cef")

shared.xgcef.p<-ggplot(shared.xgcef.p00.melt,aes(x=variable,y=value,group=time,color=time))+
  geom_jitter(position=position_jitterdodge(dodge.width = 0.7,jitter.width = 0.1)) +
  scale_y_continuous(breaks= seq(0,50,by=2),labels = function(x) insert_minor(seq(0,50,by=10), 4),limits = c(0,50)) +
  scale_x_discrete(labels=rev(labels.xgcef.p00),limits=rev(limits.xgcef.p00))+
  stat_summary(aes(group=time),fun.y=mean,fun.ymin=mean,fun.ymax=mean,geom="crossbar",width=0.15,colour="black",position=position_dodge(width=0.7)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),legend.title=element_text(size=15),legend.text=element_text(size=12)) +ylab("Relative abundance (%)") +xlab("") +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00) ) +
  coord_flip()+
  scale_color_manual(values=c("#08519c","#FFCB05"),name="Time", limits=c("n12","p00"),labels=c("Pre-Cefoperazone","Post-Cefoperazone")) +ggtitle("XG Cef")
plot(shared.xgcef.p)


######
# Std Cef Lefse Plot
lefse.preprocess(std_cef_lefse,shared.rel.t,"stdcef")
lefse.preprocess.time(shared.stdcef.tax,"n12","p00","stdcef","xg")
lefse.preprocess.time(shared.stdcef.tax,"p00","n12","stdcef","xg")

## Lefse OTUS in each taxa group
# Std Cef
# stdcef.otus<-shared.tax.stdcef.otu[,c(1,4,2,3)]
stdcef.otus.n12<-stdcef.otus[stdcef.otus$Class %in% "n12",]
# stdcef.otus.p00<-stdcef.otus[stdcef.otus$Class %in% "p00",]

limits.stdcef.n12<-c("Bacteroides","Bacteroidales_unclassified","Blautia","Butyricicoccus","Butyrivibrio","Caproiciproducens","Clostridiales_unclassified","Coprobacillus","Escherichia-Shigella","Fusicatenibacter","Lachnospiraceae_unclassified","Lactobacillus","Mollicutes_unclassified","Ruminiclostridium","Ruminococcaceae_unclassified","Ruminococcus","Syntrophococcus","Firmicutes_unclassified","Other","Coriobacteriaceae_unclassified","Eggerthella","Olsenella","Verrucomicrobiaceae_unclassified","Other.1")

labels.stdcef.n12<-c("Bacteroides (1)","Unclassified Bacteroidales (5)","Blautia (1)","Butyricicoccus (1)","Butyrivibrio (1)","Caproiciproducens (2)","Unclassified Clostridiales (4)","Coprobacillus (1)","Escherichia-Shigella (1)","Fusicatenibacter (2)","Unclassified Lachnospiraceae (49)","Lactobacillus (3)","Unclassified Mollicutes (1)","Ruminiclostridium (1)","Unclassified Ruminococcaceae (14)","Ruminococcus (1)","Syntrophococcus (3)","Unclassified Firmicutes (1)","Unclassified Firmicutes (8)","Unclassified Coriobacteriaceae (1)","Eggerthella (3)","Olsenella (2)","Unclassified Verrucomicrobiaceae (4)","Unclassified Bacteria (2)")

limits.stdcef.p00<-c("Lactobacillus")
labels.stdcef.p00<-c("Lactobacillus (1)")

# Bacteroides:"Bacteroidales_unclassified","Bacteroides",
# Firmicutes: "Blautia","Butyricicoccus","Butyrivibrio","Caproiciproducens","Clostridiales_unclassified","Coprobacillus","Escherichia-Shigella","Fusicatenibacter","Lachnospiraceae_unclassified","Lactobacillus","Mollicutes_unclassified","Ruminiclostridium","Ruminococcaceae_unclassified","Ruminococcus","Syntrophococcus","Firmicutes_unclassified","Other",
# Actinobacteria:,"Coriobacteriaceae_unclassified","Eggerthella","Olsenella",
# Verrucobacteria:"Verrucomicrobiaceae_unclassified"
# "Other.1",

shared.stdcef.n12.melt$time<-factor(shared.stdcef.n12.melt$time,levels=c("p00","n12"))
shared.stdcef.p00.melt$time<-factor(shared.stdcef.p00.melt$time,levels=c("p00","n12"))

lefse.50.plot(shared.stdcef.n12.melt,"n12","p00",limits.stdcef.n12,labels.stdcef.n12,"Pre-Cefoperazone","Post-Cefoperazone","Std Cef")

shared.stdcef.p00.p<-ggplot(shared.stdcef.p00.melt,aes(x=variable,y=value,group=time,color=time))+
  geom_jitter(position=position_jitterdodge(dodge.width = 0.8,jitter.width = 0.07)) +
  scale_y_continuous(breaks= seq(0,100,by=5),labels = function(x) insert_minor(seq(0,100,by=25), 4),limits = c(0,100)) +
  scale_x_discrete(labels=rev(labels.stdcef.p00),limits=rev(limits.stdcef.p00))+
  stat_summary(aes(group=time),fun.y=mean,fun.ymin=mean,fun.ymax=mean,geom="crossbar",width=0.07,colour="black",position=position_dodge(width=0.8)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),legend.title=element_text(size=15),legend.text=element_text(size=12)) +ylab("Relative abundance (%)") +xlab("") +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00) ) +
  coord_flip()+
  scale_color_manual(values=c("#08519c","#FFCB05"),name="Time", limits=c("n12","p00"),labels=c("Pre-Cefoperazone","Post-Cefoperazone")) +ggtitle("XG Cef LEfSe p00")
plot(shared.stdcef.p00.p)


write.table(std_cef_lefse[,c(1:4,8,5,6,7,9)],file="lefse/lefse_prepost_cef_20190419/table_s2_lefse_stdcef.tsv",sep="\t")
write.table(xg_cef_lefse[,c(1:4,8,5,6,7,9)],file="lefse/lefse_prepost_cef_20190419/table_s3_lefse_xgcef.tsv",sep="\t")
write.table(xg_diet_lefse[,c(1:4,8,5,6,7,9)],file="lefse/lefse_prepost_cef_20190419/table_s4_lefse_xgdiet.tsv",sep="\t")
