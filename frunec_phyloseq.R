# Example phyloseq analysis

# R Vannette 
# 5/21/15

#########################
# Bacterial phyloseq ##
#########################

#install.packages("phyloseq")
require(phyloseq)
require(ape)
#adding a line of code here

#change this to the folder that contains the bacterial data
setwd('~/Documents/2010 STANFORD/2015_3 summer/frunec')

# OTU table
biom<-import_biom("RMBL_otus_table_wTaxv2_fn.biom", header=T)

# change column names to match those in the sample data and taxonomy file: these need to be UNIQUE and CONSISTENT
colnames(biom)<-substring(colnames(biom),6,11)

#colnames(biom)
OTU = otu_table(biom, taxa_are_rows = TRUE)


###
# input taxonomy 
# may need to modify this for other formats, but this is formatted for RDP output.
rdptaxonomy<-read.table("~/Bioinformatics/MiSeq/2012_comb/prep/fixrank_otus_repsetOUT.fasta_classified.txt",sep="\t",
                fill=FALSE, strip.white=T)

###########
#code to split up the rdp taxonomy string #
# modify as necessary for your particular taxonomy file #
########

tax<-as.character(rdptaxonomy[-(1:6),])
#head(tax)
OTU.n<-NA
#root<-NA
Kingdom<-NA
Phylum<-NA
#Subphylum<-NA
Class<-NA
Order<-NA
Family<-NA
Genus<-NA
#Species<-NA
#Gen_ID<-NA


linaean<-data.frame(OTU.n,Kingdom,Phylum,Class,Order,Family,Genus)

for (i in 1:length(tax)){
  
  split.tax<-strsplit(tax[i],split=c(";"))
  linaean[i,1]<-substr(split.tax[[1]][1], 1,30)
  linaean[i,2]<-substr(split.tax[[1]][3], 1,30)
  linaean[i,3]<-substr(split.tax[[1]][5], 1,30)
  linaean[i,4]<-substr(split.tax[[1]][7], 1,30)
  linaean[i,5]<-substr(split.tax[[1]][9], 1,30)
  linaean[i,6]<-substr(split.tax[[1]][11], 1,30)
 # linaean[i,7]<-substr(split.tax[[1]][19], 1,30)
 # linaean[i,8]<-substr(split.tax[[1]][21], 1,30)
  #linaean[i,9]<-substr(split.tax[[1]][23], 1,30)

}
####

rownames(linaean)<- linaean$OTU.n
#head(linaean)
TAX<-tax_table(as.matrix(linaean))



# Tree file
setwd('~/Bioinformatics/MiSeq/2012_comb/')

tre<-read.tree("./prep/alignment/rep_set_tree.tre")
# may need to relabel tips to match OTU names above. 
#tre$tip.label<-paste( tre$tip.label, sep="")

# mapping file with sample info
setwd("~/Desktop/R files/2012/Jasper Ridge/")
dat<-read.csv("merged_2012_JRBP_all.csv", na.strings=c("NA", ".", "#DIV/0!"))

# files are merged by rowname, so make this consistent with above. 
rownames(dat)<-dat$Flower.ID.Match
samp<-sample_data(dat)


########################
###############
# make phyloseq object 
################

myphy<-phyloseq(OTU, samp, TAX, tre)


#head(tax_table(myphy))
noeuks<-subset_taxa(myphy, Phylum!="Cyanobacteria")
noeuks<-subset_taxa(noeuks, Kingdom!="Archaea")
noeuks<-subset_taxa(noeuks, Kingdom!="")
noeuks<-subset_taxa(noeuks, Phylum!="")
noeuks<-subset_taxa(noeuks, Family!="mitochondria")
noeuks<-subset_taxa(noeuks, Phylum!="unclassified_Root")
noeuks<-subset_taxa(noeuks, Phylum!="unclassified_Bacteria")
noeuks<-subset_taxa(noeuks, Phylum!="unclassified_Bacteria")
noeuks<-subset_taxa(noeuks, Class!="Chloroplast")
noeuks<-subset_taxa(noeuks, Class!="unclassified_Cyanobacteria/Chloroplast")
noeuks<-subset_taxa(noeuks, Class!="unclassified_Cyanobacteria/Chloroplast")

# can deal with negatives as you wish: either subtract out OTUs, subtract # of sequences, just look at them...

####################################
#further data cleaning: check to see what you're getting rid of...
####################################
# remove all samples with less than 50 total reads
noeuks = prune_samples(sample_sums(noeuks)>=50, noeuks)

zero_samples=which(taxa_sums(noeuks)==0)
length(zero_samples)
noeuks=prune_taxa(taxa_sums(noeuks)>0,noeuks)
# remove taxa not seen more than 3 times in at least 1% of samples
noeuks_f = filter_taxa(noeuks, function(x) sum(x > 3) > (0.01*length(x)), TRUE)
dim(otu_table(noeuks))
dim(otu_table(noeuks_f))


#####################################
# Examine diversity patterns
#####################################
rare<-rarefy_even_depth(noeuks_f, rngseed=1, verbose=T)
sample_data(rare)$month<-reorder.factor(sample_data(rare)$month, new.order=c("May", "June", "LJune"))

plot_richness(rare, x="Treatment", "month", measures=c("Observed", "Shannon", "InvSimpson"))+geom_boxplot()

#####################################
# Patterns of community composition
#####################################

library(vegan)
proportions<-transform_sample_counts(noeuks_f, function(x) {x/sum(x)})
proportions = prune_samples(sample_sums(proportions)>=1, proportions)
GPdist = phyloseq::distance(proportions, "bray")
GP.ord <- ordinate(proportions, "NMDS", "bray", weighted=T)

p1 = plot_ordination(proportions, GP.ord, color="LogBactCFU", shape="Treatment", 
                     title="Bacterial OTUs NMDS (Bray-Curtis)")+theme_bw()+facet_wrap(~month)
quartz()
p1

adonis(GPdist~month*Treatment,  as(sample_data(proportions), "data.frame"))


