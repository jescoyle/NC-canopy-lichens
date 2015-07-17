## This script examines community variation from Lichen Canopy Traits project

########################################################################
### Read in Data 

working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/'
sql_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Data/SQLite Tables/'

setwd(working_dir)
options(stringsAsFactors=F)

source('./Analysis/GitHub/NC-canopy-lichens/CFT_functions.R')

comm = read.csv(paste(sql_dir, 'community.csv', sep=''))
thalli = read.csv(paste(sql_dir, 'derived_traits.csv', sep=''))
samples = read.csv(paste(sql_dir, 'samples.csv', sep=''))
lichen_taxa = read.csv(paste(sql_dir, 'taxa.csv', sep=''))
trees = read.csv(paste(sql_dir, 'trees.csv', sep=''))

rownames(lichen_taxa) = lichen_taxa$TaxonID

########################################################################
### Taxonomic community composition

# Remove thalli from samples on T10-T12
use_samps = subset(samples, TreeID <=10)$SampID
comm = subset(comm, SampID %in% use_samps)
thalli = subset(thalli, SampID %in% use_samps)

# Make a site X species matrix only from thalli sampled for traits
# Note- we will make a better one below that combines taxa known only to genus
sampXsp = t(xtabs(~TaxonID+SampID, data=thalli))

rowSums(sampXsp) # richness across samples
colSums(sampXsp)

# Make a site X genus matrix
thalli$Genus = lichen_taxa[thalli$TaxonID,'Genus']
sampXgen = t(xtabs(~Genus+SampID, data=thalli))

# Make site X species matrix that includes both functional trait thalli and community samples
taxa = unique(c(comm$TaxonID, thalli$TaxonID)) # Find all unique taxon IDs
taxa = taxa[!(taxa %in% c('','Fol'))] # Remove taxa that could not be uniquely identified
taxa = taxa[lichen_taxa[taxa,'Lichenized']==1] # Remove all non-lichenized taxa
taxa = taxa[order(taxa)] # Order alphabetically

# Make taxa subsets for: only macrolichens, only foliose lichens, only crustose lichens
taxa_macro = taxa[lichen_taxa[taxa,'Form'] %in% c('foliose','fruitcose')]
taxa_fol = taxa[lichen_taxa[taxa,'Form'] == 'foliose']
taxa_crust = taxa[lichen_taxa[taxa,'Form'] == 'crustose']

use_samps = use_samps[order(use_samps)]

# Make a list of species present in whole samples (comm) and just thalli used for traits (thalli)
# This resolves (whenever possible) thalli IDed only to genus
sp_list = lapply(use_samps, function(i){
	all_present = unique(c(subset(comm, SampID==i)$TaxonID, subset(thalli, SampID==i)$TaxonID))
	trait_present = unique(subset(thalli, SampID==i)$TaxonID)
	
	# Use the full community data to resolve thalli only identified to genus (as long as only one other congeneric is present)
	trait_present = resolve_genera(trait_present, all_present, lichen_taxa)
	all_present = resolve_genera(all_present, trait_present, lichen_taxa)

	# Remove unidentified genera when they are represented by 
	all_present = calc_unique_taxa(all_present, lichen_taxa)
	trait_present = calc_unique_taxa(trait_present, lichen_taxa)

	list(comm = all_present, thalli = trait_present)
})

## Calculate sample X species matrices

# For all lichens in community (comm + thalli)
sampXsp = sapply(sp_list, function(x){
	sp = x$comm
	taxa %in% sp
})
rownames(sampXsp) = taxa
colnames(sampXsp) = use_samps
sampXsp = t(sampXsp)*1
sampXsp = sampXsp[as.numeric(rownames(sampXsp))>18,]

# For all macrolichens in community (comm + thalli)
sampXsp_macro = sapply(sp_list, function(x){
	sp = x$comm
	taxa_macro %in% sp
})
rownames(sampXsp_macro) = taxa_macro
colnames(sampXsp_macro) = use_samps
sampXsp_macro = t(sampXsp_macro)*1
sampXsp_macro = sampXsp_macro[as.numeric(rownames(sampXsp_macro))>18,]

# For all crustose in community (comm + thalli)
sampXsp_crust = sapply(sp_list, function(x){
	sp = x$comm
	taxa_crust %in% sp
})
rownames(sampXsp_crust) = taxa_crust
colnames(sampXsp_crust) = use_samps
sampXsp_crust = t(sampXsp_crust)*1
sampXsp_crust = sampXsp_crust[as.numeric(rownames(sampXsp_crust))>18,]

# For macrolichens in thalli used for trait analyses
sampXsp_thalli = sapply(sp_list, function(x){
	sp = x$thalli
	taxa_macro %in% sp
})
rownames(sampXsp_thalli) = taxa_macro
colnames(sampXsp_thalli) = use_samps
sampXsp_thalli = t(sampXsp_thalli)*1

# For foliose species in thalli used for trait analyses
sampXsp_fol = sapply(sp_list, function(x){
	sp = x$thalli
	taxa_fol %in% sp
})
rownames(sampXsp_fol) = taxa_fol
colnames(sampXsp_fol) = use_samps
sampXsp_fol = t(sampXsp_fol)*1

# Save matrices
save(sp_list, sampXsp, sampXsp_macro, sampXsp_crust, sampXsp_thalli, sampXsp_fol, 
	file='./Data/Derived Tables/sampXsp_matrices.RData')
load('./Data/Derived Tables/sampXsp_matrices.RData')


#########################################################################################
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




















