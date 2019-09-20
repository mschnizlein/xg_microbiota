######
# This code takes the sequencing output from both antibiotic cocktail and relapse experiments and summarizes the data into several outputs
## Only the abx cocktail experiment was published

##### Mothur output files and resulting modified files:
#	- alpha diversity measures --> xg_bioreac_summary.txt
#		- input file is a summary file of alpha diversity measures created in mothur
#		- output file combines these with metadata
#	- xg_bioreac.final.0.03.cons.taxonomy --> xg_bioreac.taxonomy.names.txt
#		- input file lists the taxonomic classifications for each of the OTUs at 0.03 cutoff
#		- output file was modified to reflect new OTU names that includes the taxonomy, as well as cleaner classifications
#	- xg_bioreac.final.0.03.pick.0.03.filter.summary --> xg_bioreac.betasummary.txt
#		- input file represents pairwise distances (beta diversity) generated in mothur, filtered to include only relevant samples
#		- output file adds sample meta data to sampleIDs
#	- xg_bioreac.final.0.03.pick.0.03.filter.shared, xg_bioreac_metadata.txt --> xg_bioreac_otus.w.meta.txt
#		- input .shared file was previously filtered using the specified measures
#		- output file combines metadata with OTU counts
#	- combining all data (metadata, alpha metrics, relative abundance of OTUs, and metabolites) --> xg_bioreac_allmeasures.txt
#		- all data input files
#		- output file has all variables, combined in one file
#	- --> xg_bioreac_genfrac2p.all_w.meta.txt
#		- input file is a file produced in mothur classifying sequences directly to the RDP database
#		- output file is a file with phylotype information (relative abundance of genus-level sequence assignments) 

###### alpha summaries --> xg_bioreac_summary.txt
#	- xg_bioreac.final.0.03.pick.0.03.pick.groups.summary
#	- xg_bioreac.final.0.03.pick.0.03.pick.thetayc.0.03.lt.pcoa.axes
#	- xg_bioreac.final.0.03.pick.0.03.pick.thetayc.0.03.lt.nmds.axes
#	- xg_bioreac_metadata.txt
	
#####

## To generate the metadata file, the following commands were used to grab the appropriate information from the sequence ID (DONE IN UNIX)

# sed 's/.....//' sample_names_up2.txt > timepoints.txt
# sed 's/n/-/' timepoints.txt > timepoints_neg.txt
# sed 's/p//' timepoints_neg.txt > timepoints_all.txt
# sed 's/r//' timepoints_all.txt | sed 's/c//' | sed 's/B//' > timepoints_final.txt
# sed -r 's/(.{5}).*/\1/' sample_names_up2.txt > names.txt

setwd("C:/Users/mksch/Box Sync/working_folder_schnizlein/collaborations/xanthan_gum_community_martens/experiment_5_recurrence_20180501_youngmice/mothur/")

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# With a metadata file containing the Sample and Seq IDs
meta.init<-read.table(file="metadata_xg4-5_bioreac1.txt", header=TRUE)

# Read in design file for the treatment conditions
groups.init<-read.table(file="../sequencing_exp4-5/xg_bioreac.design",header=TRUE,fill=TRUE)




# Change cec to day of experiment (this is the day for the relapse experiment; changed to correct day in the abx cocktail model while plotting)
meta.init$Time<-gsub("cec", "p28", meta.init$sampleID)
meta.init$Time<-substrRight(meta.init$Time, 3)

meta.init$Time<-gsub("p", "", meta.init$Time)
meta.init$Time<-gsub("d", "", meta.init$Time)
meta.init$Time<-gsub("n","-", meta.init$Time)

# To induce NAs into non-numeric strings in the Time column
meta.init$Time<-as.numeric(meta.init$Time)
meta.init$group<-groups.init$treatment

write.csv(meta.init, file="metadata_xg4-5_bioreac1_final.csv", quote=FALSE, row.names=FALSE)

#####
## read in files and merge together:
meta<-read.csv(file="metadata_xg4-5_bioreac1_final.csv", header=TRUE)
meta$group2<-meta$group
meta$group2<-as.character(meta$group2)
meta$group2<-as.factor(meta$group2)
pcoa<-read.table(file="xg_bioreac.final.0.03.pick.0.03.filter.thetayc.0.03.lt.pcoa.axes", header=TRUE)
	pcoa<-pcoa[,1:4]
	colnames(pcoa)[2:4] <- paste("pcoa03", colnames(pcoa)[2:4], sep = "_")
	colnames(pcoa)[1]<-"sampleID"
nmds<-read.table(file="xg_bioreac.final.0.03.pick.0.03.filter.thetayc.0.03.lt.nmds.axes", header=TRUE)
	nmds<-nmds[1:4]
	colnames(nmds)[2:4] <- paste("nmds03", colnames(nmds)[2:4], sep = "_")
	colnames(nmds)[1]<-"sampleID"
sum<-read.table(file="xg_bioreac.final.0.03.pick.groups.summary", header=TRUE)
	sum<-subset(sum, select=-c(label))
	colnames(sum)[2:16] <- paste(colnames(sum)[2:16], "03", sep = "_")
	colnames(sum)[1]<-"sampleID"
combined.pcoa<-merge(meta, pcoa, by.x=c("seqID"), by.y=c("sampleID"))
combined.nmds<-merge(combined.pcoa, nmds, by.x=c("seqID"), by.y=c("sampleID"))
combined.sum<-merge(combined.nmds, sum, by.x=c("seqID"), by.y=c("sampleID"))
write.table(combined.sum, 'xg4-5_bioreac1_summary.txt',quote=FALSE,sep="\t", col.names=TRUE, row.names=FALSE)

######
###### xg_bioreac.final.0.03.cons.taxonomy --> xg_bioreac.taxonomy.names.txt

taxonomy_file<-read.table(file="xg_bioreac.final.0.03.cons.taxonomy", header=TRUE)
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
write.table(tax, file="xanthan.taxonomy.names.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# note: some family level marks may still be messed up due to naming scheme
# these can be addressed when necessary

#####
###### Subsetting pairwise distances
#	-  xg_bioreac.final.0.03.pick.0.03.filter.0.03.pick.summary,  xanthan_metadata.txt --> xanthan.betasummary.txt
	
# read in file:
mdist<-read.table(file="xg_bioreac.final.0.03.pick.0.03.filter.summary", header=TRUE, sep="\t", row.names = NULL)
mdist<-mdist[,-c(1,4)]

# now add metadata for each sample comparison:
var<-read.table(file='xg4-5_bioreac1_summary.txt', header=TRUE)
var.sm<-var[,c(1,4,5,6)]
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
meta<-read.table(file="xg4-5_bioreac1_summary.txt", header=TRUE)
shared<-read.table(file="xg_bioreac.final.0.03.pick.0.03.filter.shared", header=TRUE, row.names=2)
dim(shared)
shared$seqID<-rownames(shared)
	# merge with meta:
sum.shared<-merge(meta, shared, by.x="seqID", by.y="seqID")
sum.shared<-subset(sum.shared, select =-c(label, numOtus) )
sum.shared<-droplevels(sum.shared)
write.table(sum.shared, file="xanthan_otus.w.meta.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#####
###### combining all data --> xg_bioreac_allmeasures.txt
#	- xg_bioreac_summary.txt
#	- xg_bioreac_metadata.txt
#	- mothurfiles/xg_bioreac.final.0.03.pick.0.03.pick.0.03.filter.shared
#####
combined.sum<-read.table(file="xanthan_otus.w.meta.txt", header=TRUE)
combined.raw<-combined.sum
	combined<-combined.raw[, c(2, 8:ncol(combined.raw))]
meta.raw<-read.table(file="xg4-5_bioreac1_summary.txt", header=TRUE)
	meta<-meta.raw[!is.na(meta.raw$patientID), ]

shared<-read.table(file="xg_bioreac.final.0.03.pick.0.03.filter.shared", header=TRUE, row.names=2)
	otu<-subset(shared, select =-c(label, numOtus) )
	otu.rel<-otu/rowSums(otu)
	colnames(otu.rel)<-paste0("rel_",colnames(otu.rel))
	otu.rel$sampleID<-rownames(otu.rel)
	
# combine all together
meta<-meta[,1:7]
combo<-merge(meta, combined, by="sampleID", all.x=TRUE)
combo2<-merge(combo, otu.rel, by="sampleID", all.x=TRUE)
write.table(combo2, file="xanthan.allmeasures.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#####
###### xg_bioreac_genfrac2p.all_w.meta.txt
#	- xg_bioreac.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary
#	- xg_bioreac_all.genera.txt	#created previously with mouse data)
#	- xg_bioreac_summary.txt (as metadata)

# step 1: create a 'phylotype' file with phylum levels
# read in mothur file; get genus-level assignments and assign phyla
tax<-read.table(file="xg_bioreac.trim.contigs.good.unique.good.filter.unique.precluster.pick.seed_v128.wang.tax.summary", header=TRUE)
# get phylum designations for level 6 (genera) rows, and curate levels for graphing (later):
tax2<-tax[which(tax$taxlevel==3), ]
tax2[, c("rankID", "taxon")]
tax6<-tax[which(tax$taxlevel==6), ]
# tax6$rankID<-gsub("^0.1.1.*", "20_Archaea_unclassified", tax6$rankID)
tax6$rankID<-gsub("^0.1.1.*", "Euryarchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.2.1\\..*", "Actinobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.2\\..*", "Bacteria_Unclassified", tax6$rankID)
tax6$rankID<-gsub("^0.2.3\\..*", "Bacteroidetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.4..*", "Cyanobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.5..*", "Firmicutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.6..*", "Fusobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.7..*", "Proteobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.8..*", "Saccharibacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.9..*", "Synergistetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.10..*", "Tenericutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.11..*", "Verrucomicrobia", tax6$rankID)
tax6$rankID<-gsub("^0.3.1..*", "unknown_unclassified", tax6$rankID)
colnames(tax6)[2]<-"phylum"
	# remove samples w/ <5000:
subtax6<-subset(tax6, select=-c(taxlevel, daughterlevels))
subtax6<-subtax6[order(subtax6$phylum, -subtax6$total), ]
taxmatrix<-subtax6[, c(4:ncol(subtax6))]
duplicated(subtax6$taxon)			#identify any duplicated taxon names
subtax6$taxon<-as.character(subtax6$taxon)

subtax6$taxon<-as.factor(subtax6$taxon)

rownames(taxmatrix)<-make.names(subtax6$taxon, unique = TRUE)
genera<- taxmatrix[, colSums(taxmatrix)>5000,]
	# get rel. abund fraction:
genmatrix<-as.data.frame(t(genera))
genera.fr<-genmatrix/rowSums(genmatrix)*100
genus.fr<-t(genera.fr)
all.genera<-cbind(subtax6[1:3], genus.fr)
write.table(all.genera, file="allxanthan_all.genera.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# note: this file included a larger dataset, and will be filtered to show the samples relevant to this project
### step 2: combine with metadata
# read in files
combined<-read.table(file="xanthan.allmeasures.txt", header = TRUE)
all.genera<-read.table(file="allxanthan_all.genera.txt", header = TRUE)
# xg_bioreac_all.genera.txt = all.genera
genbar<-all.genera
meta<-combined[, 1:7]
# filter out human samples only:
rownames(genbar)<-make.names(genbar$taxon, unique=TRUE)

mice<-genbar
phyla<-subset(genbar, select=c(phylum,taxon,total))
mice<-subset(mice, select=-c(phylum,taxon,total))
# now filter to 1 or 2%:
phyla<-subtax6[1:3]
mice[] <- lapply(mice[,1:350],as.character)
# mice.test<- mice[,4:353]
genus1<- mice[rowSums(mice>=1)>=1,]

# mice[,1:3] <- lapply(mice[,1:3],as.factor)
namelist<-as.character(rownames(genus1))
phyla1p<-phyla[phyla$taxon %in% namelist, ]
genera1<-cbind(phyla1p, genus1)
rownames(genera1)<-make.names(genera1$taxon,unique = TRUE)

# get top 2%
genus2<- mice[rowSums(mice>=2)>=2,]

namelist2p<-as.character(rownames(genus2))
phyla2p<-phyla[phyla$taxon %in% namelist2p, ]

genera2<-cbind(phyla2p, genus2)
rownames(genera2)<-make.names(genera2$taxon,unique = TRUE)

phylanames<-summary(as.factor(genera2$taxon))

#write.table(genera2, file="xg_bioreac_genfrac2p.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# read in file and combine with meta:
barg<-NULL
genbar<-genera2
	rm_g<-subset(genbar, select =-c(phylum,taxon, total) )
	barg<-as.data.frame(t(rm_g))
# 	barg[4:148,] <- sapply(barg[4:148,],as.numeric)
# 	barg$other<-100-rowSums(barg[4:148,])
#	barg<-barg[-1,]
	barg[]<-lapply(barg[], as.character)
	barg[]<-lapply(barg[], as.numeric)
	barg$Other<-100-rowSums(barg)
	barg$sampleID<-rownames(barg)
	# col.gen<-c(as.character(genbar$color), "grey47")
	barg$sampleID<-gsub("X", "", barg$sampleID)

	#meta.raw<-read.table(file="xg4-5_bioreac1_summary.txt", header=TRUE)
	#combined_wmeta<-merge(meta.raw, combined, by.x=c("sampleID"), by.y=c("sampleID"), all.y = TRUE, all.x = FALSE)
	
	# combined_wmeta<-combined_wmeta[!is.na(combined_wmeta$patientID.x), ]
# 	combined_wmeta$sampleID <- sub("^", "X", combined_wmeta$sampleID)
bar<-merge(combined, barg, by.x=c("sampleID"), by.y=c("sampleID"), all.y = TRUE, all.x = TRUE)
write.table(bar, 'xanthan_genfrac2p.all_w.meta.txt',quote=FALSE,sep="\t", col.names=NA)
# see continued analysis below

######
# if you want all genera (including the rarer guys), do this:
meta<-meta
genbar<-read.table(file="allxanthan_all.genera.txt", header=TRUE, row.names=NULL)
	#genbar5<- genbar[rowSums(genbar[ ,3:ncol(genbar)]>=5)>=5,]
	barg.all<-as.data.frame(t(mice))
	# taxon.name<-barg[1,]
	# taxon.name<-as.data.frame(t(taxon.name))
	# colnames(barg)<-taxon.name$taxon
	barg.all$sampleID<-rownames(barg.all)
	barg.all$sampleID<-gsub("X", "", barg.all$sampleID)
bar2<-merge(meta, barg.all, by.x=c("sampleID"), by.y=c("sampleID"), all.y = TRUE)
write.table(bar2, 'allxanthan_allgenera_w.meta.txt',quote=FALSE,sep="\t", col.names=NA)



