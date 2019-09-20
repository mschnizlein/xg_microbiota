# Xanthan gum mouse exp 1
# SCFA analysis 05/23/2019

# 
setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/experiment_1_20171127_youngmice/scfa_analysis/")

library(ggplot2)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n){
  substr(x, 1, n)
}

key<-read.csv("scfa_batch_key.csv")
data<-read.csv("xg_scfa_results_detA_20190517_updated.tsv.csv")

data<-merge(key, data, by.x="scfa_batch_names",by.y="Sample.ID",all.y=TRUE)
data$time<-substrLeft(data$sample_name,4)
data$time<-gsub("dp","+",data$time)
data$time<-gsub("dn","-",data$time)
data$time<-as.numeric(data$time)
data$type<-substrLeft(data$scfa_batch_names,3)
data$Area<-as.numeric(as.character(data$Area))
# Standard curves

lm.results<-NULL
lm.results<-matrix(nrow=10,ncol=3)
colnames(lm.results)<-c("compound","intercept","slope")
compound.names<-c("sialic_acid","fructose_manitol_sorbitol","lactate","formate","acetate","propionate","isobutyrate","butyrate","isovalerate","valerate")
lm.results[,1]<-compound.names

# No succinate std so excluded from analysis. Isobutyrate peak was consistently called at the wrong retention time so was excluded. All other compounds except for Acetate, propionate, and butyrate were not consistently present across all of the samples and were excluded.

stds<-subset(data, type=="STD")
stds$Sample.Name.2<-substrLeft(stds$Sample.Name,3)
stds$Sample.Name.2<-gsub("1 m","1",stds$Sample.Name.2)
stds$Sample.Name.2<-gsub("5 m","5",stds$Sample.Name.2)
stds$Sample.Name.2<-as.numeric(stds$Sample.Name.2)
stds[3,]
# excluded std samples based on linearity of their fit (i.e. if they don't fit well and are far outside of the expected value, we will exclude them since it's likely they are erroneous)
stds<-stds[-c(3),]

std.curve<-for (i in 1:10){
  std.scfa<-subset(stds, compound==compound.names[i])
  std.scfa$Sample.Name.2<-as.numeric(std.scfa$Sample.Name.2)
  std.lm<-lm(Sample.Name.2~Area,data=std.scfa)
  lm.results[i,2:3]<-std.lm$coefficients
}

lm.results<-data.frame(lm.results)
lm.results[,2:3]<-lapply(lapply(lm.results[,2:3],as.character),as.numeric)

data.scfa_conc<-for (i in 1:10){
  data.scfa<-subset(data, compound==compound.names[i])
  data.scfa$scfa_conc<-data.scfa$Area*lm.results[i,3]+lm.results[i,2]
  assign(paste("data",compound.names[i],sep="."), data.scfa)
}

data.scfa.f<-rbind(data.acetate,data.butyrate,data.formate,data.fructose_manitol_sorbitol,data.isobutyrate,data.isovalerate,data.lactate,data.propionate,data.valerate,data.sialic_acid)

data.scfa.f$group<-substrRight(substrLeft(data.scfa.f$sample_name,6),2)

for (i in c("05","07","13","38","39")){
  data.scfa.f$group<-gsub(i,"std",data.scfa.f$group)
}
for (i in c("12","11","06","37","40","44","43")){
  data.scfa.f$group<-gsub(i,"xg",data.scfa.f$group)
}

data.scfa.f$experiment<-substrRight(substrLeft(data.scfa.f$sample_name,6),2)
for (i in c("05","07","13","12","06","11")){
  data.scfa.f$experiment<-gsub(i,"2",data.scfa.f$experiment)
}
for (i in c("37","40","44","38","39","43")){
  data.scfa.f$experiment<-gsub(i,"1",data.scfa.f$experiment)
}

data.scfa.f$scfa_conc
exp1<-subset(data.scfa.f, experiment=="1")
exp2<-subset(data.scfa.f, experiment=="2")
exp1$scfa_conc_norm<-exp1$scfa_conc*10
exp2$scfa_conc_norm<-exp2$scfa_conc*10

# Limit of detection = 0.1 mM
exp1[is.na(exp1$scfa_conc_norm),]$scfa_conc_norm<-1
exp2[is.na(exp2$scfa_conc_norm),]$scfa_conc_norm<-1

# Plotting
scfa.plot<-function(dat,expnum,ident,ylimit){
  subset.data.exp2<-subset(dat, compound==ident)
  subset.data.exp2$time<-as.factor(subset.data.exp2$time)
  exp2.time.p<-ggplot(subset.data.exp2, aes(x=time, y=scfa_conc_norm,fill=group,color=group)) +
    geom_dotplot(binaxis="y",stackdir="center",position=position_dodge(width=0.7),dotsize=0.7) +
    scale_y_continuous(limits=c(0,ylimit)) +
    stat_summary(aes(group=group),fun.y=mean,fun.ymin=mean,fun.ymax=mean,geom="crossbar",width=0.3,colour="black",position=position_dodge(width=0.7)) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) + ylab("SCFA Concentration (mM)") +xlab("") +
    ggtitle(paste("Exp",expnum,ident, sep = " ")) +
    scale_color_manual(values=c("#08519c","#FFCB05"))+
  scale_fill_manual(values=c("#08519c","#FFCB05"),limits=c("std","xg"),labels=c("Standard Chow","Xanthan Gum Chow"))
  plot(exp2.time.p)
}

scfa.plot(exp2,"2","butyrate",40)
scfa.plot(exp2,"2","propionate",40)
scfa.plot(exp2,"2","acetate",400)

scfa.plot(exp1,"1","butyrate",4)
scfa.plot(exp1,"1","propionate",4)
scfa.plot(exp1,"1","acetate",20)

scfa.t.test<-function(x,t){
  t.test(subset(exp2, compound==x & time==t & group=="std")$scfa_conc,subset(exp2, compound==x & time==t & group=="xg")$scfa_conc,alternative = "two.sided")
}

# Stats
scfa.t.test("butyrate","-14")
scfa.t.test("propionate","-14")
scfa.t.test("acetate","-14")
# Day -14 not significant

scfa.t.test("butyrate","-12")
# p=0.003817
scfa.t.test("propionate","-12")
# p=0.006257
scfa.t.test("acetate","-12")
# p=0.92

scfa.t.test("butyrate","0")
scfa.t.test("propionate","0")
scfa.t.test("acetate","0")
# Day 0 not significant

#####
exp1.n01<-subset(exp1, time=="-1")
exp.n01<-ggplot(exp1.n01, aes(x=group, y=scfa_conc, fill=group)) + geom_boxplot(size=1.5) + facet_grid(.~compound) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("SCFA Concentration (mM)") +xlab("") +
  ggtitle("XG Day -1")
plot(exp.n01)


exp1.n10<-subset(exp1, time=="-10")
exp.n10<-ggplot(exp1.n10, aes(x=group, y=scfa_conc, fill=group)) + geom_boxplot(size=1.5) + facet_grid(.~compound) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("SCFA Concentration (mM)") +xlab("") +
  ggtitle("XG Day -10")
plot(exp.n10)

# exp$group2<-substrLeft(exp1$group,4)

exp2.n14<-subset(exp2, time=="-14")
exp.n14<-ggplot(exp2.n14, aes(x=group, y=scfa_conc, fill=group)) + geom_boxplot(size=1.5) + facet_grid(.~compound) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("SCFA Concentration (mM)") +xlab("") +
  ggtitle("XG Exp 2 Day -14")
plot(exp.n14)

exp2.n12<-subset(exp2, time=="-12")
exp.n12<-ggplot(exp2.n12, aes(x=group, y=scfa_conc, fill=group)) + geom_boxplot(size=1.5) + facet_grid(.~compound) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("SCFA Concentration (mM)") +xlab("") +
  ggtitle("XG Exp 2 Day -12")
plot(exp.n12)

exp2.n01<-subset(exp2, time=="-1")
exp.n01<-ggplot(exp2.n01, aes(x=group, y=scfa_conc, fill=group)) + geom_boxplot(size=1.5) + facet_grid(.~compound) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("SCFA Concentration (mM)") +xlab("") +
  ggtitle("XG Exp 2 Day -1")
plot(exp.n01)

exp2.p00<-subset(exp2, time=="0")
exp.p00<-ggplot(exp2.p00, aes(x=group, y=scfa_conc, fill=group)) + geom_boxplot(size=1.5) + facet_grid(.~compound) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("SCFA Concentration (mM)") +xlab("") +
  ggtitle("XG Exp 2 Day 0")
plot(exp.p00)

######

exp1.p<-ggplot(exp1, aes(x=group, y=scfa_conc, fill=group)) + geom_boxplot(size=1.5) + facet_grid(.~compound) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("SCFA Concentration (mM)") +xlab("") +
  ggtitle("Xanthan Gum Relapse Model")
plot(exp1.p)

exp2.subset<-exp2[exp2$compound %in% c("butyrate","propionate","lactate"),]
exp2.subset
exp2.p<-ggplot(exp2.subset, aes(x=group, y=scfa_conc, fill=group)) + geom_boxplot(size=1.5) + facet_grid(.~compound) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("SCFA Concentration (mM)") +xlab("") +
  ggtitle("Xanthan Gum Relapse Model")
plot(exp2.p)
