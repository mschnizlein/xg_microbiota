######
# Using the outputs of the "xanthan_names" run on 20180503 the following plots were made.
# Exp 1 Analysis for Publication

##### Files needed:
setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/experiment_1_20171127_youngmice/mothur_analysis/")

library(ggplot2)
library(reshape2)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n){
  substr(x, 1, n)
}

data<-read.table(file="xanthan.allmeasures.txt", header=TRUE)

data$timepoint<-gsub(-11,-12,data$timepoint)
data<-subset(data, timepoint!=-3)

data$timepoint<-factor(data$timepoint,levels=c(0,2,14,15,17,18,20,21,22,23,24))

data <- data[order(data$group2),]
unique(data$group2)

df<-data[,3:29]
# subset for a slightly smaller df to work with
# data<-droplevels(df)

######
# Diversity Metrics
data.subset.pl<-df
# data.subset.pl$group<-gsub("xanthan_cdiff","xanthan",data.subset.pl$group)
data.subset.pl<-data.subset.pl[data.subset.pl$group %in% c("standard_cdiff","xanthan_cdiff"),]

invsimp.p<-ggplot(data=data.subset.pl, aes(x=timepoint, y=invsimpson_03, color=group,group=group)) +
  scale_x_discrete(breaks=c(0,2,14,15,17,18,20,21,22,23)) +
  scale_y_continuous(limits=c(0,45)) +
  stat_summary(fun.y=mean, geom="line", size=1.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se,size=1,width=0.15)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"),limits=c("xanthan_cdiff","standard_cdiff"), labels=c("Xanthan Gum Chow","Standard Chow")) + ylab("Inverse Simpson Index") + xlab("Time (Days)") + labs(color="Legend") + ggtitle("")
plot(invsimp.p)

shannon.p<-ggplot(data=data.subset.pl, aes(x=timepoint, y=shannon_03, color=group,group=group)) +
  scale_x_discrete(breaks=c(0,2,14,15,17,18,20,21,22,23)) +
  scale_y_continuous(limits=c(0,5)) +
  stat_summary(fun.y=mean, geom="line", size=1.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se,size=1,width=0.15)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"),limits=c("xanthan_cdiff","standard_cdiff"), labels=c("Xanthan Gum Chow","Standard Chow")) + ylab("Shannon Diversity Index") + xlab("Time (Days)") + labs(color="Legend") + ggtitle("")
plot(shannon.p)

######
### stats:

# comparison of pre vs post (combined):
# donors <- data[data$human_group=="xanthan", ]
# wilcox.test(invsimpson_03~timepoint, data=donors)

t.test(subset(data.subset.pl, group=="standard_cdiff" & timepoint==0)$shannon_03,subset(data.subset.pl, group=="xanthan_cdiff" & timepoint==0)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="standard_cdiff" & timepoint==2)$shannon_03,subset(data.subset.pl, group=="xanthan_cdiff" & timepoint==2)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="standard_cdiff" & timepoint==14)$shannon_03,subset(data.subset.pl, group=="xanthan_cdiff" & timepoint==14)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="standard_cdiff" & timepoint==15)$shannon_03,subset(data.subset.pl, group=="xanthan_cdiff" & timepoint==15)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="standard_cdiff" & timepoint==17)$shannon_03,subset(data.subset.pl, group=="xanthan_cdiff" & timepoint==17)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="standard_cdiff" & timepoint==18)$shannon_03,subset(data.subset.pl, group=="xanthan_cdiff" & timepoint==18)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="standard_cdiff" & timepoint==20)$shannon_03,subset(data.subset.pl, group=="xanthan_cdiff" & timepoint==20)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="standard_cdiff" & timepoint==21)$shannon_03,subset(data.subset.pl, group=="xanthan_cdiff" & timepoint==21)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="standard_cdiff" & timepoint==22)$shannon_03,subset(data.subset.pl, group=="xanthan_cdiff" & timepoint==22)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="standard_cdiff" & timepoint==23)$shannon_03,subset(data.subset.pl, group=="xanthan_cdiff" & timepoint==23)$shannon_03,alternative="less")

# 0, 2 = p = ns
# 14,15,17,18,22 = p < 0.00001
# 21,23 = p < 0.001
# 20 standard group only has one value

### for multiple days (timepoints); NOT USED IN PAPER:
df.stat<-data[data$timepoint %in% c(-12, 0), c("timepoint", "invsimpson_03", "shannon_03")]

modelList<-list()
for(i in 2:3){
  fmla <- formula(paste(names(df.stat)[i], " ~ timepoint"))
  modelList[[i]]<-wilcox.test(fmla, data = df.stat, paired = FALSE)
}
modelList

# Inverse Simpson
# for -14 to -12:
## W = 99, p-value = 0.4892
# for -12 to -0:
## W = 271, p-value = 1.10E-6

# Shannon
# for -14 to -12:
## W = 105, p-value = 0.65
# for -12 to -0:
## W = 268, p-value = 2.45E-6
# 0,2,14,15,17,18,20,21,22,23

######

# Plotting the barplot

setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/experiment_1_20171127_youngmice/mothur_analysis/")

genus.2p<-read.table(file="xanthan_genfrac2p.all_w.meta.txt")

# The contaminated sample was still in the dataframe and needs to be removed...
genus.2p<-subset(genus.2p, sampleID!="17401p01Amis")
genus.2p<-subset(genus.2p, cecum_result!="cecum")

genus.2p$timepoint<-factor(genus.2p$timepoint,levels=c(0,2,14,15,17,18,20,21,22,23,24))

# Make the dataframe smaller so we only have the data we need to plot the families
genus.2p<-genus.2p[,c(1,4,6,235:251)]
rownames(genus.2p)<-genus.2p$sampleID

# average by group and timepoint
group_final<-aggregate(.~timepoint+group, data=genus.2p, mean)
group_final<-subset(group_final, select=-sampleID)

#####
group_final.melt<-melt(group_final, id=c("timepoint","group"))
group_final.melt$variable<-as.character(group_final.melt$variable)
group_final.melt<- group_final.melt[order(group_final.melt$variable),]

unique(group_final.melt$timepoint)

group_final.melt$variable<-as.character(group_final.melt$variable)

group_final.melt$variable2<-group_final.melt$variable
group_final.melt<-na.omit(group_final.melt)
group_final.melt$variable<-gsub("Bacteria_unclassified", "Other", group_final.melt$variable)
group_final.melt$variable<-gsub("Bacteroidales_S24.7_group", "Bacteroidales_unclassified", group_final.melt$variable)
group_final.melt$variable<-gsub("other", "Other", group_final.melt$variable)
sort(unique(group_final.melt$variable))
group_final.melt$variable<-as.factor(group_final.melt$variable)

levels.m<-c("Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Bacteroidales_unclassified","Erysipelotrichaceae","Lachnospiraceae","Lactobacillaceae","Peptostreptococcaceae","Ruminococcaceae","Clostridiales_unclassified","Actinobacteria_unclassified","Bifidobacteriaceae","Coriobacteriaceae","Enterobacteriaceae","Alphaproteobacteria_unclassified","Verrucomicrobiaceae","Verrucomicrobia_unclassified","Verrucomicrobiales_unclassified","Other")

group_final.melt$variable<-factor(group_final.melt$variable, levels=c("Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Bacteroidales_unclassified","Erysipelotrichaceae","Lachnospiraceae","Lactobacillaceae","Peptostreptococcaceae","Ruminococcaceae","Clostridiales_unclassified","Actinobacteria_unclassified","Bifidobacteriaceae","Coriobacteriaceae","Enterobacteriaceae","Alphaproteobacteria_unclassified","Verrucomicrobiaceae","Verrucomicrobia_unclassified","Verrucomicrobiales_unclassified","Other"))

labels.m<-c("Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Bacteroidales (unclassified)","Erysipelotrichaceae","Lachnospiraceae","Lactobacillaceae","Peptostreptococcaceae","Ruminococcaceae","Clostridiales (unclassified)","Actinobacteria (unclassified)","Bifidobacteriaceae","Coriobacteriaceae","Enterobacteriaceae","Alphaproteobacteria (unclassified)","Verrucomicrobiaceae","Verrucomicrobia (unclassified)","Verrucomicrobiales (unclassified)","Other")

# Colors chosen by brewer at this website: http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=9
## 4 bacteroidetes (green): "Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Bacteroidales_unclassified"
## 6 firmicutes (blue): "Erysipelotrichaceae","Lachnospiraceae","Lactobacillaceae","Peptostreptococcaceae","Ruminococcaceae","Clostridiales_unclassified"
## 3 Actinos (orange): "Actinobacteria_unclassified","Bifidobacteriaceae","Coriobacteriaceae"
## 2 proteos (yellow): "Enterobacteriaceae","Alphaproteobacteria (unclassified)"
## 3 Verrucos (purple): "Verrucomicrobiaceae","Verrucomicrobia_unclassified","Verrucomicrobiales_unclassified"
## Other (gray/black): "Other"

colors.m<-c(
  "#00441b","#006d2c","#238b45","#41ae76",
  "#023858","#045a8d","#0570b0","#3690c0","#74a9cf","#a6bddb",
  "#a63603","#d94801","#fd8d3c",
  "#C6A710","#fde157",
  "#54278f","#6a51a3","#807dba",
  "gray31")

group_final.melt$value<-as.numeric(group_final.melt$value)

xanthan.group<-subset(group_final.melt, group=="xanthan")
xanthan.group<- xanthan.group[order(xanthan.group$variable),]
xanthan_cdiff.group<-subset(group_final.melt, group=="xanthan_cdiff")
xanthan_cdiff.group<- xanthan_cdiff.group[order(xanthan_cdiff.group$variable),]
standard_cdiff.group<-subset(group_final.melt, group=="standard_cdiff")
standard_cdiff.group<- standard_cdiff.group[order(standard_cdiff.group$variable),]

insert_minor <- function(major_labs, n_minor){
  labs <- c(sapply(major_labs, function(x) c(x, rep("", n_minor))))
  labs[1:(length(labs)-n_minor)]
}

xanthan_cdiff.group.plot<-ggplot(data=xanthan_cdiff.group, aes(x =rev(timepoint), y=rev(value), group=timepoint, fill =rev(variable))) +
  geom_bar(width=0.4,stat="identity")+
  scale_y_continuous(breaks= seq(0,100,by=5),labels = function(x) paste0(insert_minor(seq(0,100,by=25), 4), "%"),limits = c(0,101)) +
  scale_x_discrete(breaks=c(0,2,4,7,10,14,15,17,18,20,21,22,23,24))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),legend.title=element_text(size=15),legend.text=element_text(size=12)) +ylab("Relative abundance (%)") +xlab("Time (Days)") +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
  scale_fill_manual(values=colors.m,name="Taxonomic Group", limits=levels.m,labels=labels.m) +ggtitle("XG Cefoperazone Model (>2% genera)") 
plot(xanthan_cdiff.group.plot)

std.cdiff.group.plot<-ggplot(data=standard_cdiff.group, aes(x =rev(timepoint), y=rev(value), group=timepoint, fill =rev(variable))) +
  geom_bar(width=0.4,stat="identity")+
  scale_y_continuous(breaks= seq(0,100,by=5),labels = function(x) paste0(insert_minor(seq(0,100,by=25), 4), "%"),limits = c(0,101)) +
  scale_x_discrete(breaks=c(0,2,4,7,10,14,15,17,18,20,21,22,23,24))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),legend.title=element_text(size=15),legend.text=element_text(size=12)) +ylab("Relative abundance (%)") +xlab("Days of Infection") +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
  scale_fill_manual(values=colors.m,name="Taxonomic Group", limits=levels.m) +ggtitle("Std+630g Cefoperazone Model (>2% genera)")
plot(std.cdiff.group.plot)

#####
# Beta Diversity
beta<-read.table(file="xanthan.betasummary.txt",header=TRUE)

beta$comparison<-as.character(beta$comparison)
beta$timepoint_comp<-substrRight(beta$comparison,3)
# beta$timepoint_comp<-gsub("c","",beta$timepoint_comp)
beta$timepoint_comp<-gsub("n","-",beta$timepoint_comp)
beta$timepoint_comp<-gsub("p","",beta$timepoint_comp)
beta<-beta[!is.na(beta$timepoint_comp),]

# exclude samples that were from cecal contents, mistakes, or duplicates
beta<-subset(beta, timepoint_comp!="-03")
beta<-subset(beta, timepoint_comp!="mis")
beta<-subset(beta, timepoint_comp!="4cr")
beta<-subset(beta, timepoint_comp!="04c")
beta<-subset(beta, timepoint_comp!="01r")
beta<-subset(beta, timepoint_comp!="08r")
beta<-subset(beta, timepoint_comp!="01B")
beta<-subset(beta, timepoint_comp!="0cr")
beta<-subset(beta, timepoint_comp!="10c")
# beta<-subset(beta, timepoint_comp!="n03")
# beta$timepoint_comp<-gsub(-11,-12,beta$timepoint_comp)
beta$timepoint_comp<-gsub("-11","-12",beta$timepoint_comp)
beta$timepoint_comp<-as.numeric(beta$timepoint_comp)

beta.n14.n12<-beta[beta$timepoint %in% c("-14","-12") & beta$timepoint_comp %in% c("-14","-12"),]
beta.n12.p00<-beta[beta$timepoint %in% c("-12","0") & beta$timepoint_comp %in% c("-12","0"),]
beta.p00.p01<-beta[beta$timepoint %in% c("0","1") & beta$timepoint_comp %in% c("0","1"),]
beta.p01.p03<-beta[beta$timepoint %in% c("1","3") & beta$timepoint_comp %in% c("1","3"),]
beta.p03.p04<-beta[beta$timepoint %in% c("3","4") & beta$timepoint_comp %in% c("3","4"),]
beta.p04.p06<-beta[beta$timepoint %in% c("4","6") & beta$timepoint_comp %in% c("4","6"),]
beta.p06.p07<-beta[beta$timepoint %in% c("6","7") & beta$timepoint_comp %in% c("6","7"),]
beta.p07.p08<-beta[beta$timepoint %in% c("7","8") & beta$timepoint_comp %in% c("7","8"),]
beta.p08.p09<-beta[beta$timepoint %in% c("8","9") & beta$timepoint_comp %in% c("8","9"),]

# beta$time_comppair<-paste(beta$timepoint,beta$timepoint_comp,sep="")

hist(beta.n14.n12[beta.n14.n12$group %in% c("xanthan_cdiff"),]$braycurtis)
hist(beta.n14.n12[beta.n14.n12$group %in% c("standard_cdiff"),]$braycurtis)

mean(beta.n14.n12[beta.n14.n12$group %in% c("xanthan_cdiff"),]$braycurtis)
mean(beta.n14.n12[beta.n14.n12$group %in% c("standard_cdiff"),]$braycurtis)

hist(beta.n14.n12[beta.n14.n12$group %in% c("xanthan_cdiff"),]$braycurtis, beta.n14.n12[beta.n14.n12$group %in% c("standard_cdiff"),]$braycurtis)
hist(beta.n12.p00[beta.n12.p00$group %in% c("xanthan_cdiff"),]$braycurtis, beta.n12.p00[beta.n12.p00$group %in% c("standard_cdiff"),]$braycurtis)

mean(beta.n12.p00[beta.n12.p00$group %in% c("xanthan_cdiff"),]$braycurtis)
mean(beta.n12.p00[beta.n12.p00$group %in% c("standard_cdiff"),]$braycurtis)

wilcox.test(beta.n14.n12[beta.n14.n12$group %in% c("xanthan_cdiff"),]$braycurtis, beta.n14.n12[beta.n14.n12$group %in% c("standard_cdiff"),]$braycurtis)
# p=0.09
wilcox.test(beta.n12.p00[beta.n12.p00$group %in% c("xanthan_cdiff"),]$braycurtis, beta.n12.p00[beta.n12.p00$group %in% c("standard_cdiff"),]$braycurtis)
# p<1e-8

mean(beta.p00.p01[beta.p00.p01$group %in% c("xanthan_cdiff"),]$braycurtis)
mean(beta.p00.p01[beta.p00.p01$group %in% c("standard_cdiff"),]$braycurtis)

wilcox.test(beta.p00.p01[beta.p00.p01$group %in% c("xanthan_cdiff"),]$braycurtis, beta.p00.p01[beta.p00.p01$group %in% c("standard_cdiff"),]$braycurtis)
# p<1e-10
mean(beta.p01.p03[beta.p01.p03$group %in% c("xanthan_cdiff"),]$braycurtis)
mean(beta.p01.p03[beta.p01.p03$group %in% c("standard_cdiff"),]$braycurtis)

wilcox.test(beta.p01.p03[beta.p01.p03$group %in% c("xanthan_cdiff"),]$braycurtis, beta.p01.p03[beta.p01.p03$group %in% c("standard_cdiff"),]$braycurtis)
# p<1e-15

mean(beta.p03.p04[beta.p03.p04$group %in% c("xanthan_cdiff"),]$braycurtis)
mean(beta.p03.p04[beta.p03.p04$group %in% c("standard_cdiff"),]$braycurtis)

wilcox.test(beta.p03.p04[beta.p03.p04$group %in% c("xanthan_cdiff"),]$braycurtis, beta.p03.p04[beta.p03.p04$group %in% c("standard_cdiff"),]$braycurtis)
# p<1e-15

mean(beta.p04.p06[beta.p04.p06$group %in% c("xanthan_cdiff"),]$braycurtis)
mean(beta.p04.p06[beta.p04.p06$group %in% c("standard_cdiff"),]$braycurtis)

wilcox.test(beta.p04.p06[beta.p04.p06$group %in% c("xanthan_cdiff"),]$braycurtis, beta.p04.p06[beta.p04.p06$group %in% c("standard_cdiff"),]$braycurtis)
# p<1e-11

mean(beta.p06.p07[beta.p06.p07$group %in% c("xanthan_cdiff"),]$braycurtis)
mean(beta.p06.p07[beta.p06.p07$group %in% c("standard_cdiff"),]$braycurtis)

wilcox.test(beta.p06.p07[beta.p06.p07$group %in% c("xanthan_cdiff"),]$braycurtis, beta.p06.p07[beta.p06.p07$group %in% c("standard_cdiff"),]$braycurtis)
# p<0.001

mean(beta.p07.p08[beta.p07.p08$group %in% c("xanthan_cdiff"),]$braycurtis)
mean(beta.p07.p08[beta.p07.p08$group %in% c("standard_cdiff"),]$braycurtis)

wilcox.test(beta.p07.p08[beta.p07.p08$group %in% c("xanthan_cdiff"),]$braycurtis, beta.p07.p08[beta.p07.p08$group %in% c("standard_cdiff"),]$braycurtis)
# p<1e-4

mean(beta.p08.p09[beta.p08.p09$group %in% c("xanthan_cdiff"),]$braycurtis)
mean(beta.p08.p09[beta.p08.p09$group %in% c("standard_cdiff"),]$braycurtis)

wilcox.test(beta.p08.p09[beta.p08.p09$group %in% c("xanthan_cdiff"),]$braycurtis, beta.p08.p09[beta.p08.p09$group %in% c("standard_cdiff"),]$braycurtis)
# p<0.001