# Xanthan Mouse Experiment Data Analysis (ALL COMBINED)
# 03/07/2019
# Matt Schnizlein
# This code brings together all of the data from the xanthan gum experiments 

library(ggplot2)
library(reshape2)
library(scales)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n){
  substr(x, 1, n)
}

setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/data_summary/")

#######
# CFU DATA

cfu_summary <- read.csv("mouse_cfu_all_experiments.csv", header = TRUE)

cfu<-cfu_summary
#####
cfu_avg<-cfu
# cfu_avg<-cfu_1[!is.na(cfu_1$d1_avg_cfu), ]
melt.cfu<-melt(cfu_avg, id=c("mouse_id","cage","model","experiment","mouse_source","group"))
melt.cfu$value<-as.character(melt.cfu$value)
melt.cfu$value<-as.numeric(melt.cfu$value)
melt.cfu$variable<-as.character(melt.cfu$variable)

melt.cfu$variable<-substrLeft(melt.cfu$variable,3)
melt.cfu<-melt.cfu[!is.na(melt.cfu$value),]


# melt$group<-as.character(melt$group)

melt.cfu[melt.cfu ==0]<-50

cef.cfu<-subset(melt.cfu, model=="cefoperazone")
abx.cfu<-subset(melt.cfu, model=="abx_cocktail")

cef.cfu$variable<-gsub("p00","14",cef.cfu$variable)
cef.cfu$variable<-gsub("p01","15",cef.cfu$variable)
cef.cfu$variable<-gsub("p02","16",cef.cfu$variable)
cef.cfu$variable<-gsub("p03","17",cef.cfu$variable)
cef.cfu$variable<-gsub("p04","18",cef.cfu$variable)
cef.cfu$variable<-gsub("p06","20",cef.cfu$variable)
cef.cfu$variable<-gsub("p07","21",cef.cfu$variable)
cef.cfu$variable<-gsub("p08","22",cef.cfu$variable)
cef.cfu$variable<-gsub("p09","23",cef.cfu$variable)
cef.cfu$variable<-gsub("p10","24",cef.cfu$variable)
cef.cfu$variable<-gsub("p11","25",cef.cfu$variable)
cef.cfu$variable<-gsub("p14","28",cef.cfu$variable)
cef.cfu$variable<-as.numeric(cef.cfu$variable)

abx.cfu$variable<-gsub("p00","8",abx.cfu$variable)
abx.cfu$variable<-gsub("p01","9",abx.cfu$variable)
abx.cfu$variable<-gsub("p03","11",abx.cfu$variable)
abx.cfu$variable<-gsub("p04","12",abx.cfu$variable)
abx.cfu$variable<-gsub("p07","15",abx.cfu$variable)
abx.cfu$variable<-gsub("p11","19",abx.cfu$variable)
abx.cfu$variable<-gsub("p14","22",abx.cfu$variable)
abx.cfu$variable<-gsub("p15","23",abx.cfu$variable)
abx.cfu$variable<-as.numeric(abx.cfu$variable)

diet.trans.cfu.cef<-subset(cef.cfu, group=="std_630_d7xan")
diet.trans.cfu.abx<-subset(abx.cfu, group=="std_630_d7xan")
diet.trans.cfu<-rbind(diet.trans.cfu.cef,diet.trans.cfu.abx)

cef.cfu<-subset(cef.cfu, group!="std_630_d7xan")
cef.cfu$group<-factor(cef.cfu$group,levels = c("std_630","xan_630","xan_neg"))
abx.cfu<-subset(abx.cfu, group!="std_630_d7xan")

# Exculuded timepoints because we only had one experimental replicate for it
cef.cfu.sub<-subset(cef.cfu, variable!="16" & variable!="20" & variable!="22" & variable!="23" & variable!="28")
# 15 16 17 18 20 21 22 23 24 25 28
cef.cfu$variable<-as.numeric(cef.cfu$variable)
# cef.cfu$group<-as.factor(cef.cfu$variable)

cfu.cef.p<-ggplot(cef.cfu.sub, aes(x=variable, y=value, color=group)) +
  stat_summary(fun.y=mean, geom="line", size=1.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se,size=1,width=0.3)+
  scale_x_discrete(limits=c(15:28), labels=c(15:29)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits=c(1e1, 1e9)) +
  annotation_logticks(base=10,sides = "l",scaled=TRUE)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) + xlab("Days Post-Infection") + ylab("CFU/g Content") +
  scale_color_manual(values=c("#08519C","#FFCB05","#1C140E"),name="Condition",breaks=c("std_630","xan_630","xan_neg"), labels=c("Standard Chow + C. difficile (N=9)","Xanthan Gum Chow + C. difficile (N=12)","Xanthan Gum Chow (N=9)")) +ggtitle("Cefoperazone (All Exps)")
plot(cfu.cef.p)

#  geom_boxplot(data=subset(cef.cfu, group=="std_630"), aes(x = variable, y=value, group=variable, fill = group), stat ="boxplot", position = "dodge") +
#  geom_boxplot(data=subset(cef.cfu, group=="xan_630"), aes(x = variable, y=value, group=variable, fill = group), stat ="boxplot", position = "dodge") +
#  geom_boxplot(data=subset(cef.cfu, group=="xan_neg"), aes(x = variable, y=value, group=variable, fill = group), stat ="boxplot", position = "dodge")

abx.cfu<-subset(abx.cfu, variable!="29")
  
cfu.abx.p<-ggplot(data=abx.cfu, aes(x=variable, y=value, group=group,color=group)) +
  stat_summary(fun.y=mean, geom="line", size=1.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se,size=1,width=0.3)+
  scale_x_discrete(limits=c(9:23), labels=c(9:23)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits=c(1e1, 1e9)) +
  annotation_logticks(base=10,sides = "l",scaled=TRUE)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) + xlab("Days Post-Infection") + ylab("CFU/g Content") +
  scale_color_manual(values=c("#08519C","#FFCB05","#1C140E"),name="Condition",breaks=c("std_630","xan_630","xan_neg"), labels=c("Standard Chow + C. difficile (N=9)","Xanthan Gum Chow + C. difficile (N=12)","Xanthan Gum Chow (N=9)")) +ggtitle("Abx Cocktail (All Exps)")
plot(cfu.abx.p)



# Stats:
unique(cef.cfu.sub$variable)
t.test(subset(cef.cfu.sub, group=="xan_630" & variable==15)$value,subset(cef.cfu.sub, group=="std_630" & variable==15)$value,alternative="less")
t.test(subset(cef.cfu.sub, group=="xan_630" & variable==17)$value,subset(cef.cfu.sub, group=="std_630" & variable==17)$value,alternative="less")
t.test(subset(cef.cfu.sub, group=="xan_630" & variable==18)$value,subset(cef.cfu.sub, group=="std_630" & variable==18)$value,alternative="less")
t.test(subset(cef.cfu.sub, group=="xan_630" & variable==21)$value,subset(cef.cfu.sub, group=="std_630" & variable==21)$value,alternative="less")
t.test(subset(cef.cfu.sub, group=="xan_630" & variable==24)$value,subset(cef.cfu.sub, group=="std_630" & variable==24)$value,alternative="less")
t.test(subset(cef.cfu.sub, group=="xan_630" & variable==25)$value,subset(cef.cfu.sub, group=="std_630" & variable==25)$value,alternative="less")

# day 15: not significant
# day 17: p = 0.01322
# day 18: p = 0.03569
# day 21: p = 0.006474
# day 24: p = 0.004417
# day 25: p = 0.0153

unique(abx.cfu$variable)
t.test(subset(abx.cfu, group=="xan_630" & variable==9)$value,subset(abx.cfu, group=="std_630" & variable==9)$value,alternative="less")
t.test(subset(abx.cfu, group=="xan_630" & variable==11)$value,subset(abx.cfu, group=="std_630" & variable==11)$value,alternative="less")
t.test(subset(abx.cfu, group=="xan_630" & variable==12)$value,subset(abx.cfu, group=="std_630" & variable==12)$value,alternative="less")
t.test(subset(abx.cfu, group=="xan_630" & variable==15)$value,subset(abx.cfu, group=="std_630" & variable==15)$value,alternative="less")
t.test(subset(abx.cfu, group=="xan_630" & variable==19)$value,subset(abx.cfu, group=="std_630" & variable==19)$value,alternative="less")
t.test(subset(abx.cfu, group=="xan_630" & variable==22)$value,subset(abx.cfu, group=="std_630" & variable==22)$value,alternative="less")
t.test(subset(abx.cfu, group=="xan_630" & variable==23)$value,subset(abx.cfu, group=="std_630" & variable==23)$value,alternative="less")

# day 9,11,12: not significant
# day 15: p = 0.01674
# day 19: p = 0.03569
# day 22: p = 0.01553
# day 23: p = 0.0153

# Diet Transition at Day 7
cfu.cefabx.p<-ggplot(data=diet.trans.cfu, aes(x=variable, y=value, group=model,color=model)) +
  stat_summary(fun.y=mean, geom="line", size=1) +
  stat_summary(geom = "errorbar", fun.data = mean_se,size=1,width=0.5,position = "dodge")+
  scale_x_discrete(limits=c(1:15), labels=c(1:15)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)),limits=c(1e1, 1e9)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) + xlab("Days Post-Infection") + ylab("CFU/g Content") +
  scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"),name="Mouse Model",breaks=c("abx_cocktail","cefoperazone"), labels=c("Antibiotic Cocktail +Clindamycin (N=6)","Cefoperazone (N=6)")) +ggtitle("Std -> XG Chow @Day 7")
plot(cfu.cefabx.p)
