
setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/qPCR/16s_qPCR_exp2/")

library(ggplot2)
library(scales)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n){
  substr(x, 1, n)
}

data.16s.p1<-read.csv(file="16s_qPCR_xg_exp2_datafile_plate1.csv")
data.16s.p2<-read.csv(file="16s_qPCR_xg_exp2_datafile_plate2.csv")
data.16s.p3<-read.csv(file="16s_qPCR_xg_exp1-2_datafile_plate3.csv")
data.16s.p4<-read.csv(file="16s_qPCR_xg_exp1_datafile_plate4.csv")
data.16s.p5<-read.csv(file="16s_qPCR_xg_exp1_datafile_plate5.csv")

data.16s.p1$plate<-1
colnames(data.16s.p1)[5]<-"Sample.Type"
data.16s.p2$plate<-2
data.16s.p3$plate<-3
data.16s.p4$plate<-4
data.16s.p5$plate<-5

data.16s<-rbind(data.16s.p1,data.16s.p2,data.16s.p3,data.16s.p4,data.16s.p5)

data.16s$Cq<-as.numeric(data.16s$Cq)
data.16s<-data.16s[!is.na(data.16s$Cq), ]

data.16s$Sample.Type<-gsub("Unknown","sample",data.16s$Sample.Type)
data.16s$Sample.Type<-gsub("Positive control","standard",data.16s$Sample.Type)
data.16s$Sample.Type<-gsub("Negative control","neg",data.16s$Sample.Type)
data.16s$Sample.Type<-gsub(" ","",data.16s$Sample.Type)
data.16s$Sample.Type<-gsub("Standard","standard",data.16s$Sample.Type)
#####
# renaming samples from plate 1
data.16s$Sample.Name<-as.factor(data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d-14 11-1","n141101",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d-14 11-2","n141102",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d-14 12-1","n141201",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d-14 12-2","n141202",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d-14 13-1","n141301",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d-14 13-2","n141302",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d-14 5-1","n140501",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d-14 5-2","n140502",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d-14 6-2","n140602",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d-14 7-1","n140701",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d-14 7-2","n140702",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 11-1","p001101",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 11-2","p001102",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 11-3","p001103",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 12-1","p001201",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 12-2","p001202",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 13-1","p001301",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 13-3","p001303",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 5-1","p000501",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 5-2","p000502",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 6-1","p000601",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 6-2","p000602",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 6-3","p000603",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d0 7-1","p000701",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d4 7-1","p040701",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d4 7-2","p040702",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d4 5-1","p040501",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d4 5-2","p040502",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d4 6-2","p040602",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d4 6-3","p040603",data.16s$Sample.Name)
data.16s$Sample.Name<-gsub("d4 7-1","p040701",data.16s$Sample.Name)

data.16s$Sample.Name<-gsub(" ","",data.16s$Sample.Name)

#####
unique(data.16s$Sample.Name)
data.16s$Cq.Mean<-as.numeric(data.16s$Cq.Mean)
data.16s$Cq.Error<-as.numeric(data.16s$Cq.Error)

data.16s.f<-subset(data.16s, select=-(fecal_mass_g))
data.16s.f<-aggregate(.~Sample.Name+Sample.Type+plate, data=data.16s.f, FUN=mean)
data.16s.f$time<-substrLeft(data.16s.f$Sample.Name,3)
# data.16s.f$time<-substrRight(data.16s.f$time,3)
# data.16s.f$time<-gsub("n","-",data.16s.f$time)
# data.16s.f$time<-gsub("p","",data.16s.f$time)
# data.16s.f$time<-as.numeric(data.16s.f$time)

data.16s.f<-subset(data.16s.f, select=c(Sample.Name,Cq,Sample.Type,time,plate))

data.16s.f$time<-gsub("n14","0",data.16s.f$time)
data.16s.f$time<-gsub("n12","2",data.16s.f$time)
data.16s.f$time<-gsub("n10","4",data.16s.f$time)
data.16s.f$time<-gsub("n07","7",data.16s.f$time)
data.16s.f$time<-gsub("n04","10",data.16s.f$time)
data.16s.f$time<-gsub("n01","13",data.16s.f$time)
data.16s.f$time<-gsub("p00","14",data.16s.f$time)
data.16s.f$time<-gsub("p04","18",data.16s.f$time)
data.16s.f$time<-gsub("p07","21",data.16s.f$time)

data.16s.f$time<-as.numeric(data.16s.f$time)

data.16s.f$time<-factor(data.16s.f$time,levels=c(0,2,4,7,10,13,14,18,21))

# Standard curves

lm.results<-NULL
lm.results<-matrix(nrow=5,ncol=2)
colnames(lm.results)<-c("intercept","slope")

std.curve<-for (i in 1:5){
  std.16s<-subset(data.16s.f, Sample.Type=="standard" & plate==i)
  std.16s<-std.16s[-5,] # remove 10^-8 because not linear
  std.16s$dna_conc<-c(9.22E+06,9.22E+04,9.22E+02,9.22E+00,9.22E+07) # assign values for std samples
  std.16s$log10_dna_conc<-log10(std.16s$dna_conc)
  std.lm<-lm(log10_dna_conc~Cq, data=std.16s) # get coefficients for standard curve
  lm.results[i,]<-std.lm$coefficients
  dna_conc<-subset(data.16s.f, plate==i & Sample.Type=="sample") # subset by plate
  dna_conc$log.dna.conc<-(lm.results[i,2]*dna_conc$Cq+lm.results[i,1]) # perform std dna calcs by plate
  dna_conc$dna.conc<-10^dna_conc$log.dna.conc
  assign(paste("data.16s_conc.p",i,sep=""), dna_conc) # output based on plate
}
# ggplot(data=std.16s.p1, aes(x=Cq, y=log10_dna_conc))+geom_point()+geom_smooth(method=lm)
data.16s.sam<-rbind(data.16s_conc.p1,data.16s_conc.p2,data.16s_conc.p3,data.16s_conc.p4,data.16s_conc.p5)

data.16s.sam<-data.16s.sam[!is.na(data.16s.sam$time), ]

data.16s.sam$group<-substrRight(substrLeft(data.16s.sam$Sample.Name,5),2)
for (i in c("05","07","13","38","39")){
  data.16s.sam$group<-gsub(i,"std",data.16s.sam$group)
}
for (i in c("12","11","06","37","40","44","43")){
  data.16s.sam$group<-gsub(i,"xg",data.16s.sam$group)
}

data.16s.sam$experiment<-substrRight(substrLeft(data.16s.sam$Sample.Name,5),2)
for (i in c("05","07","13","12","06","11")){
  data.16s.sam$experiment<-gsub(i,"2",data.16s.sam$experiment)
}
for (i in c("37","40","44","38","39","43")){
  data.16s.sam$experiment<-gsub(i,"1",data.16s.sam$experiment)
}

#add my data from fecal weights
data.fecalweights <- read.csv("xg_exp1-2_fecal_weights.csv")
#View(data.fecalweights)
data.fecalweights$Mouse.ID <- as.character(data.fecalweights$Mouse.ID)
data.fecalweights <- subset(data.fecalweights, select = c(Mouse.ID,feces.only))
data.fecalweights <- data.fecalweights[-30,]
names(data.fecalweights) <- c("Sample.Name" , "Fecal.Weight")
#merge files
data.16s.sam <- merge(data.16s.sam, data.fecalweights, by.all="Sample.Name")

#normalize the dna concentration to fecal weight
data.16s.sam$dna.conc.norm <- data.16s.sam$dna.conc/data.16s.sam$Fecal.Weight
data.16s.sam <- data.16s.sam[!is.na(data.16s.sam$dna.conc.norm), ]

# Plot the data

data.16s.sam.exp1<-subset(data.16s.sam, experiment=="1")
data.16s.sam.exp2<-subset(data.16s.sam, experiment=="2")

qpcr.plot<-ggplot(data=data.16s.sam.exp1,aes(x=time, y=dna.conc.norm, color=group, group=group)) +
  stat_summary(fun.y=mean, geom="line", size=1.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se,size=1,width=0.5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits=c(1e3,1e11)) +
  scale_x_discrete(limits=c(-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8), labels=c(-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) + xlab("Days Post-Infection") + ylab("16S Gene Copies/g feces (log10)") +
  scale_color_manual(values=c("#d95f02","#1b9e77"),name="Condition",breaks=c("std","xg"), labels=c("Standard Chow","Xanthan Gum Chow"))
plot(qpcr.plot)

qpcr.plot<-ggplot(data=data.16s.sam.exp2,aes(x=time, y=dna.conc.norm, color=group, group=group)) +
  stat_summary(fun.y=mean, geom="line", size=1.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se,size=1,width=0.15)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits=c(1e5,1e11)) +
  scale_x_discrete(breaks=c(0,2,4,7,10,13,14,18,21)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) + xlab("Days Post-Infection") + ylab("16S Gene Copies/g feces (log10)") +
  scale_color_manual(values=c("#d95f02","#1b9e77"),name="Condition",breaks=c("std","xg"), labels=c("Standard Chow","Xanthan Gum Chow"))
plot(qpcr.plot)

# 0,2,14,15,17,18,20,21,22,23,24

# Possibly use a t.test correction with p-adjust???

wilcox.test(subset(data.16s.sam.exp2, group=="std" & time==0)$dna.conc.norm,subset(data.16s.sam.exp2, group=="xg" & time==0)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam.exp2, group=="std" & time==2)$dna.conc.norm,subset(data.16s.sam.exp2, group=="xg" & time==2)$dna.conc.norm,alternative="less")
# wilcox.test(subset(data.16s.sam.exp2, group=="std" & time==4)$dna.conc.norm,subset(data.16s.sam.exp2, group=="xg" & time==4)$dna.conc.norm,alternative="less")
# wilcox.test(subset(data.16s.sam.exp2, group=="std" & time==7)$dna.conc.norm,subset(data.16s.sam.exp2, group=="xg" & time==7)$dna.conc.norm,alternative="less")
# wilcox.test(subset(data.16s.sam.exp2, group=="std" & time==10)$dna.conc.norm,subset(data.16s.sam.exp2, group=="xg" & time==10)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam.exp2, group=="std" & time==13)$dna.conc.norm,subset(data.16s.sam.exp2, group=="xg" & time==13)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam.exp2, group=="std" & time==14)$dna.conc.norm,subset(data.16s.sam.exp2, group=="xg" & time==14)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam.exp2, group=="std" & time==18)$dna.conc.norm,subset(data.16s.sam.exp2, group=="xg" & time==18)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam.exp2, group=="std" & time==21)$dna.conc.norm,subset(data.16s.sam.exp2, group=="xg" & time==21)$dna.conc.norm,alternative="less")

# 0,21 p=ns
# 2,18 p<0.05
# 13 p<0.01
# 14, p<0.001

# TO be used to test by experiment
wilcox.test(subset(data.16s.sam, group=="std" & time==0)$dna.conc.norm,subset(data.16s.sam, group=="xg" & time==0)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam, group=="std" & time==2)$dna.conc.norm,subset(data.16s.sam, group=="xg" & time==2)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam, group=="std" & time==4)$dna.conc.norm,subset(data.16s.sam, group=="xg" & time==4)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam, group=="std" & time==7)$dna.conc.norm,subset(data.16s.sam, group=="xg" & time==7)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam, group=="std" & time==10)$dna.conc.norm,subset(data.16s.sam, group=="xg" & time==10)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam, group=="std" & time==13)$dna.conc.norm,subset(data.16s.sam, group=="xg" & time==13)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam, group=="std" & time==14)$dna.conc.norm,subset(data.16s.sam, group=="xg" & time==14)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam, group=="std" & time==18)$dna.conc.norm,subset(data.16s.sam, group=="xg" & time==18)$dna.conc.norm,alternative="less")
wilcox.test(subset(data.16s.sam, group=="std" & time==21)$dna.conc.norm,subset(data.16s.sam, group=="xg" & time==21)$dna.conc.norm,alternative="less")

data.16s.sam$time.fact<-as.factor(data.16s.sam$time)
qpcr.plot<-ggplot(data=data.16s.sam,aes(x=time.fact, y=dna.conc.norm, color=group, group=group)) +
  stat_summary(fun.y=mean, geom="line", size=1.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se,size=1,width=0.2)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits=c(1e5,1e11)) +
  scale_x_discrete(breaks=c(-14,-12,-10,-7,-4,-1,0,4,7), labels=c(-14,-12,-10,-7,-4,-1,0,4,7)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) + xlab("Days Post-Infection") + ylab("16S Gene Copies/g feces (log10)") +
  scale_color_manual(values=c("#d95f02","#1b9e77"),name="Condition",breaks=c("std","xg"), labels=c("Standard Chow","Xanthan Gum Chow"))
plot(qpcr.plot)

# P-values
# -14 = 0.7977
# -12 = 0.009699
# -10 = 0.03404
# -7 = 0.08902
# -4 = 0.0834
# -1 = 0.3506
# 0 = 0.1482
# 4 = 0.04908
# 7 = 0.1485

# P-values
# -14 = 0.0504
# -12 = 0.0003912
# -10 = 0.02397
# -7 = 0.09608
# -4 = 0.06877
# -1 = 
