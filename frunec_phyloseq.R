# Analysis of Costa Rican nectars and wild fruit: chemical and microbial survey!


###############################
###############################
### PART 1: Chemistry data! ###
###############################
###############################

load('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/frunec_chem.Rdata', verbose=T)

rownames(sample_data(rare))
# chem$SampleName <- paste('fn', chem$SampleName, sep='')


# mapping file
map0 <- read.csv('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/mapping_fn.csv', header=T)

## NB: chem$SampleName is equivalent to map$chempartner

intersect(chem$SampleName, map0$chempartner)
length(intersect(chem$SampleName, map0$chempartner)) # 42 samples exist in the mapping file and the chemistry data.

setdiff(chem$SampleName, map0$chempartner)  # These are in chem but not the mapping file (i.e. they don't have partners that were sequenced.)
setdiff(map0$chempartner, chem$SampleName) # These are in the mapping file but not chem (i.e. they don't have chemistry partners.)


# Merge the chemistry data with the mapping file.
map <- merge(map0, chem, by.x='chempartner', by.y='SampleName', all.x=T)
names(map)[2] <- 'seqpartner'
rownames(map) <- map$seqpartner



############################
############################
## PART 2: Bacterial OTUs ##
############################
############################

# load packages
require(phyloseq)
require(ape)
library(ggplot2)


#change this to the folder that contains the data
setwd('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec')

# 16s OTU table
samps <- read.table('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/16s_otutable.txt',row.names=1, header=T)
mat <- data.matrix(samps, rownames.force = T)
samps <- otu_table(mat, taxa_are_rows=T)

# 16s taxonomy
tax <- read.csv('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/16s_taxonomy_RDP.csv', row.names=1, header=T)
taxmat <- as.matrix(tax,rownames.force=T)
taxtab <- tax_table(taxmat)

# mapping file
map <- sample_data(map)

# 16s tree file
tree <- read.tree('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/16s_rep_set.tre')

# glue them all together
bacteria <- merge_phyloseq(samps,taxtab,map,tree)



# Remove undesireable OTUs.
head(tax_table(bacteria))
table(tax$kingdom) # There are 18 archaeans.
table(tax$phylum) # There are 342 OTUs with no assigned phylum.  I am removing them for now.
table(tax$class) # There are 62 chloroplast OTUs
levels(tax$order)
levels(tax$family)

noeuks<-subset_taxa(bacteria, kingdom!="Unclassified")
noeuks<-subset_taxa(noeuks, kingdom!="Archaea")
noeuks<-subset_taxa(noeuks, phylum!="")
noeuks<-subset_taxa(noeuks, class!="Chloroplast")
noeuks<-subset_taxa(noeuks, family!="mitochondria")


# can deal with negatives as you wish: either subtract out OTUs, subtract # of sequences, just look at them...


IMPORTANT TO DO: look at the negative controls!  Subtract out likely contaminants.

Negative control samples are:
	pb: pcr blanks
	eb: extraction blanks
	water: water blanks
	ctab: ctab blanks
	sterile water blank: like it says.
	
	(cant remember the difference between water and sterile water blank; its written down in my field notes.)


#colnames(samps)
#otu <- otu_table(samps, taxa_are_rows = TRUE)


sum(sample_sums(noeuks)) # total number of reads = 331544
hist(sample_sums(noeuks), breaks=20)
range(sample_sums(noeuks)) 


####################################
#further data cleaning: check to see what you're getting rid of...
####################################
zero_samples=which(taxa_sums(noeuks)==0)
length(zero_samples)
noeuks=prune_taxa(taxa_sums(noeuks)>0,noeuks)

# remove taxa not seen more than 3 times in at least 1% of samples
noeuks_f = filter_taxa(noeuks, function(x) sum(x > 3) > (0.01*length(x)), TRUE)
dim(otu_table(noeuks))
dim(otu_table(noeuks_f))



# Rank abundance curve
barplot(sort(taxa_sums(rare),TRUE,)[1:100],las=2,cex.axis=.7,ylab="No. Sequences",xlab="OTU",cex=0.5,main="Rank Abundance Plot for Top 100 OTUs")


#####################################
# Examine diversity patterns
#####################################
rare <- rarefy_even_depth(noeuks_f, rngseed=1, verbose=T)
sample_data(rare)$foodtype <- reorder.factor(sample_data(rare)$foodtype, new.order=c("fruit", "nectar", "control", "eb", "pb"))

plot_richness(rare, x="species", measures=c("Observed", "Shannon", "InvSimpson"))+geom_boxplot()


#####################################
# Patterns of community composition
#####################################

library(vegan)
proportions<-transform_sample_counts(noeuks_f, function(x) {x/sum(x)})
proportions = prune_samples(sample_sums(proportions)>=1, proportions)
GPdist = phyloseq::distance(proportions, "bray")
GP.ord <- ordinate(proportions, "NMDS", "bray", weighted=T)


# ugh, I'm bad at ggplot.  can't figure out how to assign colors and shapes.

p1 = plot_ordination(proportions, GP.ord, color="species", shape="foodtype", 
                     title="Bacterial OTUs NMDS (Bray-Curtis)")+theme_bw()+facet_wrap(~month)
quartz()
p1

adonis(GPdist~month*Treatment,  as(sample_data(proportions), "data.frame"))





###########################
###########################
### PART 3: Fungal OTUs ###
###########################
###########################


# ITS OTU table
samps <- read.table('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/ITS_otutable.txt',row.names=1, header=T)
mat <- data.matrix(samps, rownames.force = T)
samps <- otu_table(mat, taxa_are_rows=T)

# ITS taxonomy
tax <- read.csv('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/ITS_taxonomy_RDPwarcup_w.confidence.csv', row.names=1, header=T)
taxmat <- as.matrix(tax,rownames.force=T)
taxtab <- tax_table(taxmat)

# The mapping file has already been loaded. The dataframe is 'map'.

# glue them all together
fungi <- merge_phyloseq(samps,taxtab,map)


# Remove undesireable OTUs.
head(tax_table(fungi))
table(tax$kingdom) # All are fungi.
table(tax$phylum) # One is incertae sedis. I'm leving it in for now.
table(tax$class) # More incertae sedis's here.
levels(tax$order)
levels(tax$family)


# can deal with negatives as you wish: either subtract out OTUs, subtract # of sequences, just look at them...


IMPORTANT TO DO: look at the negative controls!  Subtract out likely contaminants.


sum(sample_sums(fungi)) # total number of reads = 405205
hist(sample_sums(fungi), breaks=20)
range(sample_sums(fungi)) 


####################################
#further data cleaning: check to see what you're getting rid of...
####################################
zero_samples=which(taxa_sums(fungi)==0)
length(zero_samples)
fungi =prune_taxa(taxa_sums(fungi)>0, fungi)

# remove taxa not seen more than 3 times in at least 1% of samples
fungi_f = filter_taxa(fungi, function(x) sum(x > 3) > (0.01*length(x)), TRUE)
dim(otu_table(fungi))
dim(otu_table(fungi_f))



#####################################
# Examine diversity patterns
#####################################


# This is weird.  Take a closer look at this.

rare <- rarefy_even_depth(fungi_f, rngseed=1, verbose=T)
sample_data(rare)$foodtype <- reorder.factor(sample_data(rare)$foodtype, new.order=c("fruit", "nectar", "control", "eb", "pb"))

plot_richness(rare, x="species", measures=c("Observed", "Shannon", "InvSimpson"))+geom_boxplot()


# Rank abundance curve
barplot(sort(taxa_sums(rare),TRUE,)[1:100],las=2,cex.axis=.7,ylab="No. Sequences",xlab="OTU",cex=0.5,main="Rank Abundance Plot for Top 100 OTUs")


#####################################
# Patterns of community composition
#####################################

to do.









