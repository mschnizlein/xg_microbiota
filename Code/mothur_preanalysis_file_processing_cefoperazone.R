######
##### Mothur output files and resulting modified files:
#	- alpha diversity measures --> xanthan_name_summary.txt
#		- input file is a summary file of alpha diversity measures created in mothur
#		- output file combines these with metadata
#	- xanthan_name.final.0.03.cons.taxonomy --> xanthan_name.taxonomy.names.txt
#		- input file lists the taxonomic classifications for each of the OTUs at 0.03 cutoff
#		- output file was modified to reflect new OTU names that includes the taxonomy, as well as cleaner classifications
#	- xanthan_name.final.0.03.pick.0.03.filter.summary --> xanthan_name.betasummary.txt
#		- input file represents pairwise distances (beta diversity) generated in mothur, filtered to include only relevant samples
#		- output file adds sample meta data to sampleIDs
#	- xanthan_name.final.0.03.pick.0.03.filter.shared, xanthan_name_metadata.txt --> xanthan_name_otus.w.meta.txt
#		- input .shared file was previously filtered using the specified measures
#		- output file combines metadata with OTU counts
#	- combining all data (metadata, alpha metrics, relative abundance of OTUs, and metabolites) --> xanthan_name_allmeasures.txt
#		- all data input files
#		- output file has all variables, combined in one file
#	- --> xanthan_name_genfrac2p.all_w.meta.txt
#		- input file is a file produced in mothur classifying sequences directly to the RDP database
#		- output file is a file with phylotype information (relative abundance of genus-level sequence assignments) 

###### alpha summaries --> xanthan_name_summary.txt
#	- xanthan_name.final.0.03.pick.0.03.pick.groups.summary
#	- xanthan_name.final.0.03.pick.0.03.pick.thetayc.0.03.lt.pcoa.axes
#	- xanthan_name.final.0.03.pick.0.03.pick.thetayc.0.03.lt.nmds.axes
#	- xanthan_name_metadata.txt
	
#####

setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/experiment_1_20171127_youngmice/mothur_analysis/")

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n){
  substr(x, 1, n)
}

#####
## read in files and merge together:
meta<-read.table(file="r_analysis/xanthan_metadata.tsv", header=TRUE)

meta$timepoint<-as.character(meta$sampleID)
meta$timepoint<-substring(meta$timepoint,6)
meta$cecum_result<-substrRight(meta$timepoint,2)
meta$cecum_result<-gsub("cr","cecum",meta$cecum_result)
meta$cecum_result<-gsub("0c","cecum",meta$cecum_result)
meta$cecum_result<-gsub("4c","cecum",meta$cecum_result)
meta$cecum_result<-gsub("[^cecum]+", "feces", meta$cecum_result)

meta$timepoint<-substrLeft(meta$timepoint, 3)

# One sample was mislabelled as d-11 when it was actually d-12
meta$timepoint<-gsub("n11","n12",meta$timepoint)

meta$timepoint<-gsub("p", "", meta$timepoint)
meta$timepoint<-gsub("d", "", meta$timepoint)
meta$timepoint<-gsub("n","-", meta$timepoint)
meta$timepoint<-as.numeric(meta$timepoint)

meta$timepoint<-gsub("-14","a",meta$timepoint)
meta$timepoint<-gsub("-12","b",meta$timepoint)
meta$timepoint<-gsub("-3","l",meta$timepoint)
meta$timepoint<-gsub("3","e",meta$timepoint)
meta$timepoint<-gsub("4","f",meta$timepoint)
meta$timepoint<-gsub("6","g",meta$timepoint)
meta$timepoint<-gsub("7","h",meta$timepoint)
meta$timepoint<-gsub("8","i",meta$timepoint)
meta$timepoint<-gsub("9","j",meta$timepoint)
meta$timepoint<-gsub("10","k",meta$timepoint)
meta$timepoint<-gsub("1","d",meta$timepoint)
meta$timepoint<-gsub("0","c",meta$timepoint)

meta$timepoint<-gsub("a","0",meta$timepoint)
meta$timepoint<-gsub("b","2",meta$timepoint)
meta$timepoint<-gsub("c","14",meta$timepoint)
meta$timepoint<-gsub("d","15",meta$timepoint)
meta$timepoint<-gsub("e","17",meta$timepoint)
meta$timepoint<-gsub("f","18",meta$timepoint)
meta$timepoint<-gsub("g","20",meta$timepoint)
meta$timepoint<-gsub("h","21",meta$timepoint)
meta$timepoint<-gsub("i","22",meta$timepoint)
meta$timepoint<-gsub("j","23",meta$timepoint)
meta$timepoint<-gsub("k","24",meta$timepoint)
meta$timepoint<-gsub("l","-3",meta$timepoint)

meta$group2<-meta$group
meta$group2<-as.character(meta$group2)
meta$group2<-as.factor(meta$group2)
pcoa<-read.table(file="xanthan_name.final.0.03.pick.0.03.filter.thetayc.0.03.lt.pcoa.axes", header=TRUE)
	pcoa<-pcoa[,1:4]
	colnames(pcoa)[2:4] <- paste("pcoa03", colnames(pcoa)[2:4], sep = "_")
	colnames(pcoa)[1]<-"sampleID"
nmds<-read.table(file="xanthan_name.final.0.03.pick.0.03.filter.thetayc.0.03.lt.nmds.axes", header=TRUE)
	nmds<-nmds[1:4]
	colnames(nmds)[2:4] <- paste("nmds03", colnames(nmds)[2:4], sep = "_")
	colnames(nmds)[1]<-"sampleID"
sum<-read.table(file="xanthan_name.final.0.03.pick.groups.summary", header=TRUE)
	sum<-subset(sum, select=-c(label))
	colnames(sum)[2:16] <- paste(colnames(sum)[2:16], "03", sep = "_")
	colnames(sum)[1]<-"sampleID"
combined.pcoa<-merge(meta, pcoa, by.x=c("seqID"), by.y=c("sampleID"))
combined.nmds<-merge(combined.pcoa, nmds, by.x=c("seqID"), by.y=c("sampleID"))
combined.sum<-merge(combined.nmds, sum, by.x=c("seqID"), by.y=c("sampleID"))
write.table(combined.sum, 'xanthan_summary.txt',quote=FALSE,sep="\t", col.names=TRUE, row.names=FALSE)

#####
## xanthan_name.final.0.03.cons.taxonomy --> xanthan_name.taxonomy.names.txt

taxonomy_file<-read.table(file="xanthan_name.final.0.03.cons.taxonomy", header=TRUE)
# taxname to OTU:
tax <- taxonomy_file$Taxonomy
tax <- gsub("\\(\\d*\\)", "", tax)
tax <- gsub(";unclassified", "", tax)
tax <- gsub("_1", "", tax)
tax <- gsub(";$", "", tax)
tax <- gsub("/.*", "", tax)
tax <- gsub(".*;", "", tax)
tax.names <-paste(taxonomy_file$OTU, tax, sep="_")
tax.names <-gsub("000", "", tax.names)
taxonomy_file$taxname<-tax.names
# phylum variable:
phylum <- taxonomy_file$Taxonomy
phylum <- gsub("\\(\\d*\\)", "", phylum)
phylum <- gsub("Bacteria;", "", phylum)
phylum <- gsub(";$", "", phylum)
phylum <- gsub(";.*", "", phylum)
taxonomy_file$phylum<-phylum
# family variable:
fam <- taxonomy_file$Taxonomy
fam <- gsub("\\(\\d*\\)", "", tax)
fam <- gsub(";unclassified", "", tax)
fam <- gsub("_1", "", tax)
fam <- gsub(";$", "", tax)
fam <- gsub("/.*", "", tax)
fam <- gsub(".*les;", "", tax)
fam <- gsub(";.*", "", tax)
taxonomy_file$family<-fam
write.table(taxonomy_file, file="xanthan.taxonomy.names.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# note: some family level marks may still be messed up due to naming scheme
# these can be addressed when necessary

#####
## Subsetting pairwise distances specific to human files
#	-  xanthan_name.final.0.03.pick.0.03.filter.0.03.pick.summary,  xanthan_metadata.txt --> xanthan.betasummary.txt
	
mdist<-read.table(file="xanthan_name.final.0.03.pick.0.03.filter.summary", header=TRUE, sep="\t", row.names = NULL)
mdist<-mdist[,-c(1,4)]

# now add metadata for each sample comparison:
var<-read.table(file='xanthan_summary.txt', header=TRUE)
var.sm<-var[,c(1,4,6)]
var.comp<-var.sm
colnames(var.comp)<-paste0(colnames(var.comp),"_comp")	
# merge files:
m1<-merge(var.sm, mdist, by.x=c("seqID"), by.y=c("label")) #this merges the data based on the sampleID/group1 match
m2<-merge(var.comp, m1, by.x=c("seqID_comp"), by.y=c("comparison")) #this merges the data based on the sampleID/group1 match

write.table(m2, file="xanthan.betasummary.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#####
# Creating an OTU count file:
## heatmap with metadata and taxonomy, and clustering:
	# read in files:
meta<-read.table(file="xanthan_summary.txt", header=TRUE)
shared<-read.table(file="xanthan_name.final.0.03.pick.0.03.filter.shared", header=TRUE, row.names=2)
dim(shared)
shared$seqID<-rownames(shared)

# merge with meta:
sum.shared<-merge(meta, shared, by.x="seqID", by.y="seqID")
sum.shared<-subset(sum.shared, select =-c(label, numOtus) )
sum.shared<-droplevels(sum.shared)

# Output gives raw OTU sequence counts
write.table(sum.shared, file="xanthan_otus.w.meta.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#####
# Creating OTU rel abund file:
meta<-read.table(file="xanthan_summary.txt", header=TRUE)
shared<-read.table(file="xanthan_name.final.0.03.pick.0.03.filter.shared", header=TRUE, row.names=2)
	otu<-subset(shared, select =-c(label, numOtus) )
	otu.rel<-otu/rowSums(otu)
	otu.rel$sampleID<-rownames(otu.rel)

combo2<-merge(meta, otu.rel, by="sampleID", all.x=TRUE)

write.table(combo2, file="xanthan.allmeasures.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#####
###### xanthan_name_genfrac2p.all_w.meta.txt
#	- xanthan_name.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary
#	- xanthan_name_all.genera.txt	#created previously with mouse data)
#	- xanthan_name_summary.txt (as metadata)

# step 1: create a 'phylotype' file with phylum levels
# read in mothur file; get genus-level assignments and assign phyla
tax<-read.table(file="xanthan_name.trim.contigs.good.unique.good.filter.unique.precluster.pick.seed_v128.wang.tax.summary", header=TRUE)
# get phylum designations for level 6 (genera) rows, and curate levels for graphing (later):
tax3<-tax[which(tax$taxlevel==3), ]
tax3[, c("rankID", "taxon")]
tax6<-tax[which(tax$taxlevel==6), ]

tax6$rankID<-gsub("^0.1.1.*", "20_Euryarchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.2.1\\..*", "04_Actinobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.2\\..*", "11_Bacteria_Unclassified", tax6$rankID)
tax6$rankID<-gsub("^0.2.3\\..*", "01_Bacteroidetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.4\\..*", "20_Cyanobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.5\\..*", "20_Firmicutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.6\\..*", "20_Fusobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.7\\..*", "20_Proteobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.8\\..*", "20_Saccharibacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.9\\..*", "20_Tenericutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.10..*", "20_Verrucomicrobia", tax6$rankID)

colnames(tax6)[2]<-"phylum"
	# remove samples w/ <5000:
subtax6<-subset(tax6, select=-c(taxlevel, daughterlevels))
subtax6<-subtax6[order(subtax6$phylum, -subtax6$total), ]
taxmatrix<-subtax6[, c(4:ncol(subtax6))]

duplicated(subtax6$taxon)			#identify any duplicated taxon names
subtax6$taxon<-as.character(subtax6$taxon)
subtax6$taxon[15]<-"Actinobacteria_unclassified2"
subtax6$taxon[30]<-"Cyanobacteria_unclassified2"

subtax6$taxon<-as.factor(subtax6$taxon)


rownames(taxmatrix)<-make.names(subtax6$taxon)
genera<- taxmatrix[, colSums(taxmatrix)>5000,]
genera<-genera[,1:143]
	# get rel. abund fraction:
genmatrix<-as.data.frame(t(genera))
genera.fr<-genmatrix/rowSums(genmatrix)*100
genus.fr<-t(genera.fr)
all.genera<-cbind(subtax6[1:3], genus.fr)
write.table(all.genera, file="allxanthan_all.genera.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

### step 2: combine with metadata and filter out data relevant to human stuff:
# read in files
combined<-read.table(file="xanthan.allmeasures.txt", header = TRUE)
meta<-combined[, 1:9]

all.genera<-read.table(file="allxanthan_all.genera.txt", header = TRUE)
genbar<-all.genera

rownames(genbar)<-make.names(genbar$taxon, unique=TRUE)

mice<-genbar
phyla<-subset(genbar, select=c(phylum,taxon,total))
mice<-subset(mice, select=-c(phylum,taxon,total))
# now filter to 1 or 2%:
mice[] <- lapply(mice[,1:143],as.numeric)
genus1<- mice[rowSums(mice>=1)>=1,]

# mice[,1:3] <- lapply(mice[,1:3],as.factor)
# namelist1p<-as.character(rownames(genus1))
# phyla1p<-phyla[phyla$taxon %in% namelist1p, ]
# genera1<-cbind(phyla, genus1)

# get top 2%
genus2<- mice[rowSums(mice>=2)>=2,]

namelist2p<-as.character(rownames(genus2))
phyla2p<-phyla[phyla$taxon %in% namelist2p, ]

genera2<-genus2

#write.table(genera2, file="xanthan_name_genfrac2p.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# read in file and combine with meta:
barg<-NULL
genbar<-genera2
#	rm_g<-subset(genbar, select =-c(phylum, total) )
#	barg<-as.data.frame(t(rm_g))
barg<-as.data.frame(t(genbar))
# 	barg[4:148,] <- sapply(barg[4:148,],as.numeric)
# 	barg$other<-100-rowSums(barg[4:148,])
#	barg<-barg[-1,]
	barg[]<-lapply(barg[], as.character)
	barg[]<-lapply(barg[], as.numeric)
	barg$other<-100-rowSums(barg)
	barg$sampleID<-rownames(barg)
	# col.gen<-c(as.character(genbar$color), "grey47")
	barg$sampleID<-gsub("X", "", barg$sampleID)

bar<-merge(combined, barg, by.x=c("sampleID"), by.y=c("sampleID"), all.y = TRUE, all.x = TRUE)
write.table(bar, 'xanthan_genfrac2p.all_w.meta.txt',quote=FALSE,sep="\t", col.names=NA)

######
# if you want all genera (including the rarer guys), do this:
meta<-meta
genbar<-read.table(file="allxanthan_all.genera.txt", header=TRUE, row.names=NULL)
	#genbar5<- genbar[rowSums(genbar[ ,3:ncol(genbar)]>=5)>=5,]
	rm_g<-subset(genbar, select =-c(phylum, total) )
	
	barg<-as.data.frame(t(rm_g))
	taxon.name<-barg[1,]
	taxon.name<-as.data.frame(t(taxon.name))
	colnames(barg)<-taxon.name$taxon
	barg$sampleID<-rownames(barg)
bar2<-merge(meta, barg, by.x="sampleID", by.y="sampleID", all.y = TRUE)
write.table(bar2, 'allxanthan_allgenera_w.meta.txt',quote=FALSE,sep="\t", col.names=NA)