## This script examines community variation from Lichen Canopy Traits project





########################################################################
### Read in Data 

working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/'
sql_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Data/SQLite Tables/'

setwd(working_dir)
options(stringsAsFactors=F)

source('./Analysis/CFT_functions.R')

comm = read.csv(paste(sql_dir, 'community.csv', sep=''))
thalli = read.csv(paste(sql_dir, 'derived_traits.csv', sep=''))
samples = read.csv(paste(sql_dir, 'samples.csv', sep=''))
lichen_taxa = read.csv(paste(sql_dir, 'taxa.csv', sep=''))
trees = read.csv(paste(sql_dir, 'trees.csv', sep=''))

rownames(lichen_taxa) = lichen_taxa$TaxonID

########################################################################
### Taxonomic community composition

# Remove thalli from samples on T10-T12
use_samps = subset(samples, TreeID <10)$SampID
comm = subset(comm, SampID %in% use_samps)
thalli = subset(thalli, SampID %in% use_samps)

# Make a site X species matrix only from thalli sampled for traits
sampXsp = t(xtabs(~TaxonID+SampID, data=thalli))

rowSums(sampXsp) # richness across samples
colSums(sampXsp)

# Make a site X genus matrix
thalli$Genus = lichen_taxa[thalli$TaxonID,'Genus']
sampXgen = t(xtabs(~Genus+SampID, data=thalli))

# Make site X species matrix that includes both functional trait thalli and community samples
# NOTE: NEED TO GO THROUGH AND ID A COUPLE MISSING THALLI
taxa = unique(c(comm$TaxonID, thalli$TaxonID))
taxa = subset(lichen_taxa, 
taxa = taxa[order(taxa)]
use_samps = use_samps[order(use_samps)]

sampXsp = sapply(use_samps, function(i){
	present_taxa = unique(c(subset(comm, SampID==i)$TaxonID, subset(thalli, SampID==i)$TaxonID))
	
	taxa %in% present_taxa
})
sampXsp = t(sampXsp)*1
rownames(sampXsp) = as.character(use_samps)
colnames(sampXsp) = taxa



### Vertical abundance profiles for species and genera



# Calculate distance from top of tree
samples$Height_top = sapply(1:nrow(samples), function(i){
	trees[trees$TreeID==samples[i,'TreeID'],'Height']-samples[i,'Height']
}) 

height_cat = c(0,10,15,20,25)

# Calculate relative abundance of each species in each height class
height_bin = cut(samples[rownames(sampXsp),'Height'], height_cat)

binXsp = aggregate(1:nrow(sampXsp), list(Height=height_bin), function(i){
	these_samps = sampXsp[i,]
	colSums(these_samps)/(nrow(these_samps)*24)*100
})$x
rownames(binXsp) = levels(height_bin)

height_bin = cut(samples[rownames(sampXgen),'Height'], height_cat)
binXgen = aggregate(1:nrow(sampXgen), list(Height=height_bin), function(i){
	these_samps = sampXgen[i,]
	colSums(these_samps)/(nrow(these_samps)*24)*100
})$x
rownames(binXgen) = levels(height_bin) 

# Actually, rather than summing all samples makes more sense to treat each sample in a height bin as an observation in order to get error bars
plot(




















