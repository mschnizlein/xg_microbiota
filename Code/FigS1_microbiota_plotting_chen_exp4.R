# XG Manuscript
# Chen Mouse Model Microbiota Analysis
# Fig S1

library(ggplot2)
library(reshape2)
library(dplyr)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n){
  substr(x, 1, n)
}

######
# Using the outputs of the "xg_bioreac" sequencing run in 201809 the following plots were made

##### Files needed:
setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/experiment_4_20180405_youngmice")

data<-read.table(file="../experiment_5_recurrence_20180501_youngmice/mothur/xanthan.allmeasures.txt", header=TRUE)

data <- data[order(data$group),]
unique(data$group2)
unique(data$Time)
data$Time<-as.character(data$Time)
data$Time<-as.numeric(data$Time)

# Set time relative to the experiment start (rather than C. difficile inoculation)

unique(data$Time)
data <- data[order(data$Time),]
unique(data$Time)



# subset for a slightly smaller df to work with
df_pre<-data[,3:28]

df_pre<-droplevels(df_pre)

## code to loop something over time, split by group:
# names(data)[7]<-"group"

colnames(df_pre)<-gsub("\\.x","", colnames(df_pre))

# Subset your data into the appropriate groups

chen<-subset(df_pre, df_pre$experiment=="chen_model")
chen$Time<-gsub(28,14, chen$Time)
chen$Time<-as.numeric(chen$Time)

df<-chen

unique(df$Time)

# genus.2p$timepoint<-gsub("-8","a",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("-6","b",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("-1","c",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("0","d",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("1","e",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("4","f",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("7","g",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("14","h",genus.2p$timepoint)

# genus.2p$timepoint<-gsub("a","0",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("b","2",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("c","7",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("d","8",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("e","9",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("f","12",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("g","15",genus.2p$timepoint)
# genus.2p$timepoint<-gsub("h","22",genus.2p$timepoint)

# Plot Each of the three diversity measures
data.subset.pl<-df
data.subset.pl<-data.subset.pl[data.subset.pl$group %in% c("std","xg_cdiff"),]

invsimp.p<-ggplot(data=data.subset.pl, aes(x=Time, y=invsimpson_03, color=group,group=group)) +
  scale_x_discrete(limits=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14), labels=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14)) +
  scale_y_continuous(limits=c(0,45)) +
  stat_summary(fun.y=mean, geom="line", size=1.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se,size=1,width=0.15)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"),limits=c("xg_cdiff","std"), labels=c("Xanthan Gum Chow","Standard Chow")) + ylab("Inverse Simpson Index") + xlab("Time (Days)") + labs(color="Legend") + ggtitle("")
plot(invsimp.p)

shannon.p<-ggplot(data=data.subset.pl, aes(x=Time, y=shannon_03, color=group,group=group)) +
  scale_x_discrete(limits=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14), labels=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14)) +
  scale_y_continuous(limits=c(0,5)) +
  stat_summary(fun.y=mean, geom="line", size=1.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se,size=1,width=0.15)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"),limits=c("xg_cdiff","std"), labels=c("Xanthan Gum Chow","Standard Chow")) + ylab("Shannon Diversity Index") + xlab("Time (Days)") + labs(color="Legend") + ggtitle("")
plot(shannon.p)

# Plotting all three groups (not used in publication)
sobs<-ggplot(data=df, aes(x=Time, y=sobs_03, color=group)) + geom_jitter() +
  scale_x_discrete(limits=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14), labels=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14)) +
  stat_summary(fun.y=mean, geom="line") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  scale_color_hue(limits=c("std", "xg_cdiff", "xg"), labels=c("Std Chow +630g","5% XG Chow +630g","5% XG Chow")) + ylab("SOBS") + xlab("Time (Days)") + labs(color="Legend") + ggtitle("")
plot(sobs)

shannon<-ggplot(data=df, aes(x=Time, y=shannon_03, color=group)) + geom_jitter() +
  scale_x_discrete(limits=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14), labels=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14)) +
  stat_summary(fun.y=mean, geom="line") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  scale_color_hue(limits=c("std", "xg_cdiff", "xg"), labels=c("Std Chow +630g","5% XG Chow +630g","5% XG Chow")) + ylab("Shannon") + xlab("Time (Days)") + labs(color="Legend") + ggtitle("")
plot(shannon)

inv.simp<-ggplot(data=df, aes(x=Time, y=invsimpson_03, color=group)) + 
  scale_x_discrete(limits=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14), labels=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14)) +
  stat_summary(fun.y=mean, geom="line",size=1.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se,size=1,width=0.5)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  scale_color_hue(limits=c("std", "xg_cdiff", "xg"), labels=c("Standard Chow","5% Xanthan Gum Chow","Negative Control")) + ylab("Inverse Simpson") + xlab("Time (Days)") + labs(color="Legend") + ggtitle("")
plot(inv.simp)

#######
### stats:
t.test(subset(data.subset.pl, group=="std" & Time==-8)$shannon_03,subset(data.subset.pl, group=="xg_cdiff" & Time==-8)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="std" & Time==-6)$shannon_03,subset(data.subset.pl, group=="xg_cdiff" & Time==-6)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="std" & Time==-1)$shannon_03,subset(data.subset.pl, group=="xg_cdiff" & Time==-1)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="std" & Time==0)$shannon_03,subset(data.subset.pl, group=="xg_cdiff" & Time==0)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="std" & Time==1)$shannon_03,subset(data.subset.pl, group=="xg_cdiff" & Time==1)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="std" & Time==4)$shannon_03,subset(data.subset.pl, group=="xg_cdiff" & Time==4)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="std" & Time==7)$shannon_03,subset(data.subset.pl, group=="xg_cdiff" & Time==7)$shannon_03,alternative="less")
t.test(subset(data.subset.pl, group=="std" & Time==14)$shannon_03,subset(data.subset.pl, group=="xg_cdiff" & Time==14)$shannon_03,alternative="less")

# Day -8,-6,-1,0 all not significant
# Day 1 p = 0.0037
# Day 4 p = 0.0001515
# Day 7 p = 1e-7
# Day 14 p = 0.0001751

######
# Plotting the barplot

setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/experiment_4_20180405_youngmice/")

genus.2p<-read.table(file="../experiment_5_recurrence_20180501_youngmice/mothur/xanthan_genfrac2p.all_w.meta.txt")

genus.2p<-genus.2p[,c(1:29, 410:433)]
genus.2p<-genus.2p[,c(1:6, 30:53)]

chen<-subset(genus.2p, genus.2p$experiment=="chen_model")
chen$Time<-gsub(28,14, chen$Time)
chen$Time<-as.numeric(chen$Time)

df<-chen

df<-df[,5:30]

# average by group2 and timepoint
group_final<-aggregate(.~Time+group,data=df,mean)

# for (i in 3:99){
#  colnames(group_final)[i]<-ifelse(colSums(group_final[i])<=2,"Other", colnames(group_final[i]))
#}
#colnames(group_final)<-make.unique(colnames(group_final))

######
group_final.melt<-melt(group_final, id=c("Time","group"))
group_final.melt$variable<-gsub("\\..*","", group_final.melt$variable)
group_final.melt<-aggregate(value~Time+group+variable,data=group_final.melt,sum)
# group_final.melt$mouseID<-as.factor(group_final.melt$mouseID)
group_final.melt$variable<-as.character(group_final.melt$variable)
group_final.melt<- group_final.melt[order(group_final.melt$variable),]
unique(group_final.melt$Time)
# unique(genus.2p.group.melt$mouseID)
group_final.melt$variable<-as.character(group_final.melt$variable)

group_final.melt$variable2<-group_final.melt$variable
# group_final.melt$variable[group_final.melt$value<=1.99]<-"Other"
# group_final.melt$variable[group_final.melt$value<=0.99]<-NA
group_final.melt<-na.omit(group_final.melt)
group_final.melt$variable<-gsub("Bacteria_unclassified", "Other", group_final.melt$variable)
group_final.melt$variable<-gsub('Bacteroidales_S24', "Bacteroidales_unclassified",group_final.melt$variable)
group_final.melt$variable<-gsub("Family_XIII","Other",group_final.melt$variable)

# group_final.melt$variable2[group_final.melt$value<=0.001]<-"Other"
sort(unique(group_final.melt$variable))
levels.m<-c("Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Bacteroidales_unclassified","Enterococcaceae","Erysipelotrichaceae","Lachnospiraceae","Lactobacillaceae","Paenibacillaceae","Peptostreptococcaceae","Ruminococcaceae","Staphylococcaceae","Clostridiales_unclassified","Bifidobacteriaceae","Coriobacteriaceae","Actinobacteria_unclassified","Enterobacteriaceae","Pseudomonadaceae","Verrucomicrobiaceae","Verrucomicrobiales_unclassified","Other")

group_final.melt$variable<-factor(group_final.melt$variable, levels=c("Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Bacteroidales_unclassified","Enterococcaceae","Erysipelotrichaceae","Lachnospiraceae","Lactobacillaceae","Paenibacillaceae","Peptostreptococcaceae","Ruminococcaceae","Staphylococcaceae","Clostridiales_unclassified","Bifidobacteriaceae","Coriobacteriaceae","Actinobacteria_unclassified","Enterobacteriaceae","Pseudomonadaceae","Verrucomicrobiaceae","Verrucomicrobiales_unclassified","Other"))

labels.m<-c("Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Bacteroidales (unclassified)","Enterococcaceae","Erysipelotrichaceae","Lachnospiraceae","Lactobacillaceae","Paenibacillaceae","Peptostreptococcaceae","Ruminococcaceae","Staphylococcaceae","Clostridiales (unclassified)","Bifidobacteriaceae","Coriobacteriaceae","Actinobacteria (unclassified)","Enterobacteriaceae","Pseudomonadaceae","Verrucomicrobiaceae","Verrucomicrobiales_unclassified","Other")

# Colors chosen by brewer at this website: http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=9
## 4 bacteroidetes (green): "Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Bacteroidales_unclassified"
## 6 firmicutes (blue): "Enterococcaceae","Erysipelotrichaceae","Lachnospiraceae","Lactobacillaceae","Paenibacillaceae","Peptostreptococcaceae","Ruminococcaceae","Staphylococcaceae","Clostridiales_unclassified"
## 3 Actinos (orange): "Bifidobacteriaceae","Coriobacteriaceae","Actinobacteria_unclassified",
## 2 proteos (yellow): "Enterobacteriaceae","Pseudomonadaceae",
## 3 Verrucos (purple): "Verrucomicrobiaceae","Verrucomicrobiales_unclassified"
## Other (gray/black): "Other"






# "#023858","#045a8d","#0570b0","#3690c0","#74a9cf","#a6bddb",





"#6baed6","#9ecae1",


colors.m<-c(
  "#00441b","#006d2c","#238b45","#41ae76",
  "#023858","#045a8d","#0570b0","#3690c0","#74a9cf","#a6bddb","#08519c","#2171b5","#4292c6",
  "#a63603","#d94801","#fd8d3c",
  "#C6A710","#fde157",
  "#54278f","#807dba",
  "gray31")

# relapse<-relapse[!is.na(relapse$Time), ]

group_final.melt$value<-as.numeric(group_final.melt$value)

group_final.melt$Time<-factor(group_final.melt$Time,levels=c(-8,-6,-1,0,1,4,7,14))

xg.group<-subset(group_final.melt, group=="xg")
xg.group<- xg.group[order(xg.group$variable),]
xg_cdiff.group<-subset(group_final.melt, group=="xg_cdiff")
xg_cdiff.group<- xg_cdiff.group[order(xg_cdiff.group$variable),]
std_cdiff.group<-subset(group_final.melt, group=="std")
std_cdiff.group<- std_cdiff.group[order(std_cdiff.group$variable),]

insert_minor <- function(major_labs, n_minor){
  labs <- c(sapply(major_labs, function(x) c(x, rep("", n_minor))))
  labs[1:(length(labs)-n_minor)]
}

xg.group$variable<-factor(xg.group$variable, levels=rev(c("Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Bacteroidales_unclassified","Enterococcaceae","Erysipelotrichaceae","Lachnospiraceae","Lactobacillaceae","Paenibacillaceae","Peptostreptococcaceae","Ruminococcaceae","Staphylococcaceae","Clostridiales_unclassified","Bifidobacteriaceae","Coriobacteriaceae","Actinobacteria_unclassified","Enterobacteriaceae","Pseudomonadaceae","Verrucomicrobiaceae","Verrucomicrobiales_unclassified","Other")))

xg.group.plot<-ggplot(data=xg.group, aes(x =rev(Time), y=rev(value), group=Time, fill =rev(variable))) +
  geom_bar(width=0.4,stat="identity")+
  scale_y_continuous(breaks= seq(0,100,by=5),labels = function(x) paste0(insert_minor(seq(0,100,by=25), 4), "%"),limits = c(0,101)) +
  scale_x_discrete(breaks=c(-8,-6,-1,0,1,4,7,14), labels=c(-8,-6,-1,0,1,4,7,14))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),legend.title=element_text(size=15),legend.text=element_text(size=12)) +ylab("Relative abundance (%)") +xlab("Time (Days)") +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
  scale_fill_manual(values=colors.m,name="Taxonomic Group", limits=levels.m,labels=labels.m) +ggtitle("XG Antibiotic Cocktail Model (>2% genera)")
plot(xg.group.plot)

xg_cdiff.group.plot<-ggplot(data=xg_cdiff.group, aes(x =rev(Time), y=rev(value), group=Time, fill =rev(variable))) +
  geom_bar(width=0.4,stat="identity")+
  scale_y_continuous(breaks= seq(0,100,by=5),labels = function(x) paste0(insert_minor(seq(0,100,by=25), 4), "%"),limits = c(0,101)) +
  scale_x_discrete(breaks=c(-8,-6,-1,0,1,4,7,14), labels=c(-8,-6,-1,0,1,4,7,14))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),legend.title=element_text(size=15),legend.text=element_text(size=12)) +ylab("Relative abundance (%)") +xlab("Time (Days)") +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
  scale_fill_manual(values=colors.m,name="Taxonomic Group", limits=levels.m,labels=labels.m) +ggtitle("XG+Cd Antibiotic Cocktail Model (>2% genera)") 
plot(xg_cdiff.group.plot)

std_cdiff.group.plot<-ggplot(data=std_cdiff.group, aes(x =rev(Time), y=rev(value), group=Time, fill =rev(variable))) +
  geom_bar(width=0.4,stat="identity")+
  scale_y_continuous(breaks= seq(0,100,by=5),labels = function(x) paste0(insert_minor(seq(0,100,by=25), 4), "%"),limits = c(0,101)) +
  scale_x_discrete(breaks=c(-8,-6,-1,0,1,4,7,14), labels=c(-8,-6,-1,0,1,4,7,14))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),legend.title=element_text(size=15),legend.text=element_text(size=12)) +ylab("Relative abundance (%)") +xlab("Time (Days)") +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
  scale_fill_manual(values=colors.m,name="Taxonomic Group", limits=levels.m,labels=labels.m) +ggtitle("Std Antibiotic Cocktail Model (>2% genera)") 
plot(std_cdiff.group.plot)


#####
# Beta Diversity
beta<-read.table(file="../experiment_5_recurrence_20180501_youngmice/mothur/xanthan.betasummary.txt",header=TRUE)

beta$seqID<-as.character(beta$seqID)
beta$seqID_comp<-as.character(beta$seqID_comp)
beta$filter<-substrRight(beta$seqID,3)
beta$filter_comp<-substrRight(beta$seqID_comp,3)
beta<-beta[beta$experiment %in% "chen_model" & beta$experiment_comp %in% "chen_model",]

# beta$Time_comp<-gsub("c","",beta$Time_comp)
# beta$Time_comp<-gsub("n","-",beta$Time_comp)
# beta$Time_comp<-gsub("p","",beta$Time_comp)
beta$filter<-as.character(beta$filter)
beta$filter_comp<-as.character(beta$filter_comp)
# beta$Time<-subset(beta,  filter != "cec" & filter_comp != "cec" & Time != "28" & Time_comp != "28")


beta$Time<-gsub("28","14",beta$Time)
beta$Time_comp<-gsub("28","14",beta$Time_comp)
beta<-beta[!is.na(beta$Time_comp),]

beta$Time<-as.numeric(beta$Time)
beta$Time_comp<-as.numeric(beta$Time_comp)

beta.nn.nn<-function(d1,d2){
  beta.day<-beta[beta$Time %in% c(d1,d2) & beta$Time_comp %in% c(d1,d2),]
  beta.day<-beta.day[!(beta.day$Time == d1 & beta.day$Time_comp == d1) & !(beta.day$Time == d2 & beta.day$Time_comp == d2), ]
  assign(paste("beta",gsub("-","n0",d1),gsub("-","n0",d2),sep="."), beta.day,envir = globalenv())
}

beta.nn.nn(-8,-6)
beta.nn.nn(-6,-1)
beta.nn.nn(-1,0)
beta.nn.nn(0,1)
beta.nn.nn(1,4)
beta.nn.nn(4,7)
beta.nn.nn(7,14)

# beta$time_comppair<-paste(beta$Time,beta$Time_comp,sep="")

hist(beta.n08.n06[beta.n08.n06$group %in% c("xg_cdiff"),]$braycurtis)
hist(beta.n08.n06[beta.n08.n06$group %in% c("std"),]$braycurtis)

beta.n08.n06.xg<-beta.n08.n06[beta.n08.n06$group %in% c("xg_cdiff") & beta.n08.n06$group_comp %in% c("xg_cdiff"),]$braycurtis
beta.n08.n06.std<-beta.n08.n06[beta.n08.n06$group %in% c("std") & beta.n08.n06$group_comp %in% c("std"),]$braycurtis

mean(beta.n08.n06.xg)
mean(beta.n08.n06.std)
wilcox.test(beta.n08.n06.xg,beta.n08.n06.std)
# p = 0.002997

beta.n06.n01.xg<-beta.n06.n01[beta.n06.n01$group %in% c("xg_cdiff") & beta.n06.n01$group_comp %in% c("xg_cdiff"),]$braycurtis
beta.n06.n01.std<-beta.n06.n01[beta.n06.n01$group %in% c("std") & beta.n06.n01$group_comp %in% c("std"),]$braycurtis

mean(beta.n06.n01.xg)
mean(beta.n06.n01.std)
wilcox.test(beta.n06.n01.xg,beta.n06.n01.std)
# p =5.407e-8

beta.n01.0.xg<-beta.n01.0[beta.n01.0$group %in% c("xg_cdiff") & beta.n01.0$group_comp %in% c("xg_cdiff"),]$braycurtis
beta.n01.0.std<-beta.n01.0[beta.n01.0$group %in% c("std") & beta.n01.0$group_comp %in% c("std"),]$braycurtis

mean(beta.n01.0.xg)
mean(beta.n01.0.std)
wilcox.test(beta.n01.0.xg,beta.n01.0.std)
# p = 0.01288

beta.0.1.xg<-beta.0.1[beta.0.1$group %in% c("xg_cdiff") & beta.0.1$group_comp %in% c("xg_cdiff"),]$braycurtis
beta.0.1.std<-beta.0.1[beta.0.1$group %in% c("std") & beta.0.1$group_comp %in% c("std"),]$braycurtis

mean(beta.0.1.xg)
mean(beta.0.1.std)
wilcox.test(beta.0.1.xg,beta.0.1.std)
# p = 0.003999

beta.1.4.xg<-beta.1.4[beta.1.4$group %in% c("xg_cdiff") & beta.1.4$group_comp %in% c("xg_cdiff"),]$braycurtis
beta.1.4.std<-beta.1.4[beta.1.4$group %in% c("std") & beta.1.4$group_comp %in% c("std"),]$braycurtis

mean(beta.1.4.xg)
mean(beta.1.4.std)
wilcox.test(beta.1.4.xg,beta.1.4.std)
# p = 0.04592

beta.4.7.xg<-beta.4.7[beta.4.7$group %in% c("xg_cdiff") & beta.4.7$group_comp %in% c("xg_cdiff"),]$braycurtis
beta.4.7.std<-beta.4.7[beta.4.7$group %in% c("std") & beta.4.7$group_comp %in% c("std"),]$braycurtis

mean(beta.4.7.xg)
mean(beta.4.7.std)
wilcox.test(beta.4.7.xg,beta.4.7.std)
# p = 0.1787

# Cecal Contents
beta.7.14.xg<-beta.7.14[beta.7.14$group %in% c("xg_cdiff") & beta.7.14$group_comp %in% c("xg_cdiff"),]$braycurtis
beta.7.14.std<-beta.7.14[beta.7.14$group %in% c("std") & beta.7.14$group_comp %in% c("std"),]$braycurtis

mean(beta.7.14.xg)
mean(beta.7.14.std)
wilcox.test(beta.7.14.xg,beta.7.14.std)

#######
# PCoA Plots (not used in publication)

#Color palette by treatment group over time:

pcoa1<-ggplot(data=df, aes(y=pcoa03_axis1, x=Time, color=group)) + geom_jitter() +
  scale_x_discrete(limits=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14), labels=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14)) +
  stat_summary(fun.y=mean, geom="line") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  scale_color_hue(limits=c("std", "xg_cdiff", "xg"), labels=c("Std Chow +630g","5% XG Chow +630g","5% XG Chow")) + ylab("PCoA Axis 1 (8.0%)") + xlab("Time (Days)") + labs(color="Legend") + ggtitle("")
plot(pcoa1)

pcoa2<-ggplot(data=df, aes(y=pcoa03_axis2, x=Time, color=group)) + geom_jitter() +
  scale_x_discrete(limits=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14), labels=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14)) +
  stat_summary(fun.y=mean, geom="line") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  scale_color_hue(limits=c("std", "xg_cdiff", "xg"), labels=c("Std Chow +630g","5% XG Chow +630g","5% XG Chow")) + ylab("PCoA Axis 2 (5.0%)") + xlab("Time (Days)") + labs(color="Legend") + ggtitle("")
plot(pcoa2)

pcoa3<-ggplot(data=df, aes(y=pcoa03_axis3, x=Time, color=group)) + geom_jitter() +
  scale_x_discrete(limits=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14), labels=c(-8,-6,-4,-2,0,2,4,6,8,10,12,14)) +
  stat_summary(fun.y=mean, geom="line") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  scale_color_hue(limits=c("std", "xg_cdiff", "xg"), labels=c("Std Chow +630g","5% XG Chow +630g","5% XG Chow")) + ylab("PCoA Axis 3 (3.2%)") + xlab("Time (Days)") + labs(color="Legend") + ggtitle("")
plot(pcoa3)

# pcoa1.2<-ggplot(data=df, aes(y=pcoa03_axis1, x=pcoa03_axis2, color=group)) + geom_jitter() +
#  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
#  scale_color_hue(limits=c("rel", "rel_cdiff", "rel_cdiff_fmt","rel_cdiff_xg"), labels=c("No intervention, No 630g","No intervention, 630g", "FMT, 630g","XG Chow, 630g")) + ylab("PCoA Axis 1 (8.0%)") + xlab("PCoA Axis 2 (5.0 %)") + labs(color="Legend") + ggtitle("")
# plot(pcoa1.2)

# note: you can add the % variance explained per axis by looking at the .loadings file produced alongside your .axes file in mothur