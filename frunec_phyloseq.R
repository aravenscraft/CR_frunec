# Analysis of Costa Rican nectars and wild fruit: chemical and microbial survey!

############################
############################
## PART 1: Bacterial OTUs ##
############################
############################

# load packages
require(phyloseq)
require(ape)
library(ggplot2)


#change this to the folder that contains the bacterial data
setwd('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec')

# 16s OTU table
samps <- read.table('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/frunec_16s_otutable.txt',row.names=1, header=T)
mat <- data.matrix(samps, rownames.force = T)
samps <- otu_table(mat, taxa_are_rows=T)

# 16s taxonomy
tax <- read.csv('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/repset_tax_assignments.csv', row.names=1, header=T)
taxmat <- as.matrix(tax,rownames.force=T)
taxtab <- tax_table(taxmat)

# 16s mapping file
map <- read.csv('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/mapping_fn.csv', row.names=1, header=T)
map <- sample_data(map)

# 16s tree file
tree <- read.tree('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/rep_set.tre')

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


sum(sample_sums(noeuks)) # total number of reads
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
### PART 2: Fungal OTUs ###
###########################
###########################

# TO DO:
Load the fungal OTU table here!







###############################
###############################
### PART 3: Chemistry data! ###
###############################
###############################

load('~/Documents/2010 STANFORD/2015_3 summer/frunec/CR_frunec/frunec_chem.Rdata', verbose=T)

rownames(sample_data(rare))
chem$SampleName <- paste('fn', chem$SampleName, sep='')

## NB: these don't match because each of the chem samples has a corresponding "partner" sample that was sequenced.  Partners were usually from different flowers on the same plant, collected at the same time.

# TO DO:
Load the list of correpsonding partners, then merge the chemistry data into the phyloseq mapping file.







