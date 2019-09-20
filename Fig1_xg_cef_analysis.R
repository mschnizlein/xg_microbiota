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

beta$seqID<-as.character(beta$seqID)
beta$seqID_comp<-as.character(beta$seqID_comp)
beta$filter<-substrRight(beta$seqID,3)
beta$filter_comp<-substrRight(beta$seqID_comp,3)
# beta$timepoint_comp<-gsub("c","",beta$timepoint_comp)
#beta$timepoint_comp<-gsub("n","-",beta$timepoint_comp)
#beta$timepoint_comp<-gsub("p","",beta$timepoint_comp)
# beta<-beta[!is.na(beta$timepoint_comp),]

# exclude samples that were from cecal contents, mistakes, or duplicates
beta<-subset(beta, !(filter %in% c("-03","mis","4cr","04c","01r","08r","01B","0cr","10c")))
beta<-subset(beta, !(filter_comp %in% c("-03","mis","4cr","04c","01r","08r","01B","0cr","10c")))


# beta<-subset(beta, timepoint_comp!="n03")
# beta$timepoint_comp<-gsub(-11,-12,beta$timepoint_comp)
# beta$timepoint_comp<-gsub("-11","-12",beta$timepoint_comp)
beta$timepoint<-as.numeric(beta$timepoint)
beta$timepoint_comp<-as.numeric(beta$timepoint_comp)

beta.nn.nn<-function(d1,d2){
  beta.day<-beta[beta$timepoint %in% c(d1,d2) & beta$timepoint_comp %in% c(d1,d2),]
    beta.day<-beta.day[!(beta.day$timepoint == d1 & beta.day$timepoint_comp == d1) & !(beta.day$timepoint == d2 & beta.day$timepoint_comp == d2), ]
        assign(paste("beta",d1,d2,sep="."), beta.day,envir = globalenv())
}

beta.nn.nn(0,2)
beta.nn.nn(2,14)
beta.nn.nn(14,15)
beta.nn.nn(15,17)
beta.nn.nn(17,18)
beta.nn.nn(18,20)
beta.nn.nn(20,21)
beta.nn.nn(21,22)
beta.nn.nn(22,23)

# beta$time_comppair<-paste(beta$timepoint,beta$timepoint_comp,sep="")

hist(beta.0.2[beta.0.2$group %in% c("xanthan_cdiff") & beta.0.2$group_comp %in% c("xanthan_cdiff"),]$braycurtis)
hist(beta.0.2[beta.0.2$group %in% c("standard_cdiff") & beta.0.2$group_comp %in% c("standard_cdiff"),]$braycurtis)

# beta.0.2[beta.0.2$group %in% c("xanthan_cdiff") & beta.0.2$group_comp %in% c("xanthan_cdiff") & beta.0.2$braycurtis <0.6,]

beta.0.2.xg<-beta.0.2[beta.0.2$group %in% c("xanthan_cdiff") & beta.0.2$group_comp %in% c("xanthan_cdiff"),]$braycurtis
beta.0.2.std<-beta.0.2[beta.0.2$group %in% c("standard_cdiff") & beta.0.2$group_comp %in% c("standard_cdiff"),]$braycurtis

mean(beta.0.2.xg)
mean(beta.0.2.std)
wilcox.test(beta.0.2.xg,beta.0.2.std)
# p=0.87

# hist(beta.00.02[beta.00.02$group %in% c("xanthan_cdiff") & beta.00.02$group_comp %in% c("xanthan_cdiff"),]$braycurtis, beta.00.02[beta.00.02$group %in% c("standard_cdiff") & beta.00.02$group_comp %in% c("standard_cdiff"),]$braycurtis)
# hist(beta.02.14[beta.02.14$group %in% c("xanthan_cdiff"),]$braycurtis, beta.02.14[beta.02.14$group %in% c("standard_cdiff"),]$braycurtis)

beta.2.14.xg<-beta.2.14[beta.2.14$group %in% c("xanthan_cdiff") & beta.2.14$group_comp %in% c("xanthan_cdiff"),]$braycurtis
beta.2.14.std<-beta.2.14[beta.2.14$group %in% c("standard_cdiff") & beta.2.14$group_comp %in% c("standard_cdiff"),]$braycurtis

mean(beta.2.14.xg)
mean(beta.2.14.std)
wilcox.test(beta.2.14.xg,beta.2.14.std)
# p<2.2e-16

beta.14.15.xg<-beta.14.15[beta.14.15$group %in% c("xanthan_cdiff") & beta.14.15$group_comp %in% c("xanthan_cdiff"),]$braycurtis
beta.14.15.std<-beta.14.15[beta.14.15$group %in% c("standard_cdiff") & beta.14.15$group_comp %in% c("standard_cdiff"),]$braycurtis

mean(beta.14.15.xg)
mean(beta.14.15.std)
wilcox.test(beta.14.15.xg,beta.14.15.std)
# p=0.97
beta.15.17.xg<-beta.15.17[beta.15.17$group %in% c("xanthan_cdiff") & beta.15.17$group_comp %in% c("xanthan_cdiff"),]$braycurtis
beta.15.17.std<-beta.15.17[beta.15.17$group %in% c("standard_cdiff") & beta.15.17$group_comp %in% c("standard_cdiff"),]$braycurtis

mean(beta.15.17.xg)
mean(beta.15.17.std)
wilcox.test(beta.15.17.xg,beta.15.17.std)
# p=0.0006

beta.17.18.xg<-beta.17.18[beta.17.18$group %in% c("xanthan_cdiff") & beta.17.18$group_comp %in% c("xanthan_cdiff"),]$braycurtis
beta.17.18.std<-beta.17.18[beta.17.18$group %in% c("standard_cdiff") & beta.17.18$group_comp %in% c("standard_cdiff"),]$braycurtis

mean(beta.17.18.xg)
mean(beta.17.18.std)
wilcox.test(beta.17.18.xg,beta.17.18.std)
# p=0.6

beta.18.20.xg<-beta.18.20[beta.18.20$group %in% c("xanthan_cdiff") & beta.18.20$group_comp %in% c("xanthan_cdiff"),]$braycurtis
beta.18.20.std<-beta.18.20[beta.18.20$group %in% c("standard_cdiff") & beta.18.20$group_comp %in% c("standard_cdiff"),]$braycurtis

# beta.18.20[!(beta.18.20$timepoint == "18" & beta.18.20$timepoint_comp == "18"), ]

mean(beta.18.20.xg)
mean(beta.18.20.std)
wilcox.test(beta.18.20.xg,beta.18.20.std)
# p=0.1429

beta.20.21.xg<-beta.20.21[beta.20.21$group %in% c("xanthan_cdiff") & beta.20.21$group_comp %in% c("xanthan_cdiff"),]$braycurtis
beta.20.21.std<-beta.20.21[beta.20.21$group %in% c("standard_cdiff") & beta.20.21$group_comp %in% c("standard_cdiff"),]$braycurtis

mean(beta.20.21.xg)
mean(beta.20.21.std)
wilcox.test(beta.20.21.xg,beta.20.21.std)
# p=0.1667

beta.21.22.xg<-beta.21.22[beta.21.22$group %in% c("xanthan_cdiff") & beta.21.22$group_comp %in% c("xanthan_cdiff"),]$braycurtis
beta.21.22.std<-beta.21.22[beta.21.22$group %in% c("standard_cdiff") & beta.21.22$group_comp %in% c("standard_cdiff"),]$braycurtis

mean(beta.21.22.xg)
mean(beta.21.22.std)
wilcox.test(beta.21.22.xg,beta.21.22.std)
# p=0.607

beta.22.23.xg<-beta.22.23[beta.22.23$group %in% c("xanthan_cdiff") & beta.22.23$group_comp %in% c("xanthan_cdiff"),]$braycurtis
beta.22.23.std<-beta.22.23[beta.22.23$group %in% c("standard_cdiff") & beta.22.23$group_comp %in% c("standard_cdiff"),]$braycurtis

mean(beta.22.23.xg)
mean(beta.22.23.std)
wilcox.test(beta.22.23.xg,beta.22.23.std)
# p=0.1483
