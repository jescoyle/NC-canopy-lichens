## This script examines community variation from Lichen Canopy Traits project

########################################################################
### Read in Data 

working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/'
sql_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Data/SQLite Tables/'
data_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Data/Derived Tables/'

setwd(working_dir)
options(stringsAsFactors=F)

source('./Analysis/GitHub/NC-canopy-lichens/CFT_functions.R')

comm = read.csv(paste(sql_dir, 'community.csv', sep=''))
thalli = read.csv(paste(sql_dir, 'derived_traits.csv', sep=''))
samples = read.csv(paste(sql_dir, 'samples.csv', sep=''))
lichen_taxa = read.csv(paste(sql_dir, 'taxa.csv', sep=''))
trees = read.csv(paste(sql_dir, 'trees.csv', sep=''))
env = read.csv(paste(data_dir, 'loggerdata.csv', sep=''), row.names=1)

rownames(lichen_taxa) = lichen_taxa$TaxonID
env$SampID = rownames(env)

# Load sample mean traits
load('./Analysis/sample mean traits.RData')

# Set working colors
mycolor = read.csv('../blue2red_10colramp.txt')[9:2,]
mycolor = rgb(mycolor, maxColorValue=255)# Set working colors

# Calculate distance from top of tree
samples$Height_top = sapply(1:nrow(samples), function(i){
	trees[trees$TreeID==samples[i,'TreeID'],'Height']-samples[i,'Height']
}) 


########################################################################
### Taxonomic community composition

# Remove thalli from samples on T11-T12
use_samps = subset(samples, TreeID <=10)$SampID
comm = subset(comm, SampID %in% use_samps)
thalli = subset(thalli, SampID %in% use_samps)

# Substitute NA for bad trait values
thalli[thalli$ThallusID %in% c(419,501,521,541), 'STA'] = NA # Values abnormally high- area measured prob doesn't correspond to mass
thalli[thalli$ThallusID==8, 'Water_capacity'] = NA # Can't have negative water capacity
thalli[(thalli$Chla2b<0) & (!is.na(thalli$Chla2b)),'Chla2b'] = NA # Can't have negative ratios. Measurements probably erroneous
thalli[thalli$Genus=='Usnea','Rhizine_length'] = NA # Usnea don't have rhizines

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
taxa_macro = taxa[lichen_taxa[taxa,'Form'] %in% c('foliose','fruticose')]
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


# Define predictors
env_vars = c('Light_mean','Light_high','Temp_max','Vpd_mean','Vpd_daysatfreq')
xvarnames = data.frame(var = env_vars)
xvarnames = expression('Mean Light Intensity (Lux)', 'Freq. Full Sun', 
	'Max. Temperature (C*degree)','Mean VPD','Freq. Daytime Saturated Air')
names(xvarnames) = env_vars

xvarnames_short = c('Mean Light','Full Sun Freq.','Max Temp.','Mean VPD','Sat. Air Freq.')
names(xvarnames_short) = env_vars


#########################################################################################
### Ordination

library(vegan)

# Define positional predictors
pos_vars = c('Height_top','Angle')

# Define variables
Y = sampXsp_thalli
X = env[as.character(use_samps),env_vars]
#X = samples[as.character(use_samps), pos_vars]
W = factor(samples[as.character(use_samps),'Year'])
names(W) = use_samps

# Scale X data so that variances are similar
X[,'Light_mean'] = X[,'Light_mean']/1000
X[,'Light_high'] = X[,'Light_high']*100
X[,'Vpd_daysatfreq'] = X['Vpd_daysatfreq']*100

# Remove species with no observed occurence in this data set
Y = Y[,colSums(Y)>0]

# Remove ambiguous taxa: NOT DONE
#Y = Y[,lichen_taxa[colnames(Y),'TaxonConcept']!='genus']

# Remove samples with no observed species
Y = Y[rowSums(Y)>0,]

# Make sure rows are in same order (i.e. from the same samples)
X = X[rownames(Y),]
W = W[rownames(Y)]

# Make a matrix that defines which species belong to which genera
taxon_hier = sapply(colnames(Y), function(x) colnames(Y) == substr(x, 1, 3))
rownames(taxon_hier) = colnames(Y)

# Define Hellinger distance function that finds 0 distance between genus sp. and species in the same genus.
# i and j are vectors of species abundances
# taxon_dependency is a binary matrix whose columns indicate the taxonomic hierarchy
#	rownames are the taxa that are ambiguous (i.e. the genera)
calc_Hellinger = function(i, j, taxon_dependency){
	# Determine which taxa are ambiguous
	ambig_taxa = which(colnames(taxon_dependency) %in% rownames(taxon_dependency))
	
	# Determine whether there is a lichen only identified to genus
	if(sum(i[ambig_taxa])>0){
	for(n in names(which(i[ambig_taxa]>0))){
		# If the other site does have a species in that genus
		# Distribute individuals among those species according to their abundance distribution at that site
		dependents = names(which(j[taxon_dependency[n,]]>0))

		if(length(dependents)>0){
			Ni = i[n]
			Nj = sum(j[dependents])
			i[dependents] = Ni*j[dependents]/Nj
			if(n!=dependents) i[n] = 0 # Checks whether dependent taxon is actually the same as the ambiguous taxon
		}
	}}

	if(sum(j[ambig_taxa])>0){
	for(n in names(which(j[ambig_taxa]>0))){
		# If the other site does have a species in that genus
		# Distribute individuals among those species according to their abundance distribution at that site
		dependents = names(which(i[taxon_dependency[n,]]>0))

		if(length(dependents)>0){
			Nj = j[n]
			Ni = sum(i[dependents])
			j[dependents] = Nj*i[dependents]/Ni
			if(n!=dependents) j[n] = 0
		}
	
	}}
	
	# Calculate Hellinger distance
	sqrt(sum((sqrt(i/sum(i)) - sqrt(j/sum(j)))^2))
}

Dmat = matrix(NA, nrow=nrow(Y), ncol=nrow(Y))
rownames(Dmat) = rownames(Y)
colnames(Dmat) = rownames(Y)
for(i in rownames(Dmat)){
for(j in colnames(Dmat)){
	Dmat[i,j] = calc_Hellinger(Y[i,], Y[j,], taxon_hier[rowSums(taxon_hier)>0,])
}}

# PCoA on Dmat to get community matrix
pcoa = cmdscale(Dmat, k = min(ncol(Y), nrow(Y)-1), eig=T)
Y.p = pcoa$points

## Constrained ordination on presence absence data

# db-RDA on Hellinger distance matrix (PCoA followed by RDA)
ord1 = rda(Y.p~., data=X)

# CCA - constrained ordination on pres-abs data preserving chi-sq distance
ord1 = cca(Y~., data=X)

anova(ord1) # Test significance of whole model
anova(ord1, by='axis')
anova(ord1, by='term')

spscores = scores(ord1)$species
spscores[order(spscores[,'CCA1']),]
spscores[order(spscores[,'CCA2']),]

plot(ord1, display='sites', type='n', las=1)
points(ord1, display='sites', pch=1, col='black', bg='lightblue')
text(ord1, display='bp',  lwd=2)#labels=xvarnames_short,
points(ord1, display='species', pch=3, lwd=2)

plot(ord1, display='sites', type='n', las=1)
points(ord1, display='sites', pch=21, col='black', bg='lightblue')
text(ord1, 'species', select=apply(spscores, 1, function(x) sqrt(x[1]^2 + x[2]^2))>1.5)
text(ord1, display='bp', labels=c('Height','Angle'),  lwd=2, col=2)
text(ord1, 'species', select=apply(spscores, 1, function(x) sqrt(x[1]^2 + x[2]^2))>1.5)

## Partial out effects of Year
# db-RDA
ord2 = rda(Y.p~.+Condition(W), data=X)

# CCA
ord2 = cca(Y~.+Condition(W), data=X)

ord2
anova(ord2)
anova(ord2, by='axis')
anova(ord2, by='term')

spscores = scores(ord2)$species

plot(ord2, display='sites', type='n', las=1)
points(ord2, display='sites', pch=1, col='black', bg='lightblue')
text(ord2, display='bp', labels=xvarnames_short, lwd=2)
points(ord2, display='species', pch=3, lwd=2)

plot(ord2, display='sites', type='n', las=1)
points(ord2, display='sites', pch=1, col='black', bg='lightblue')
text(ord2, 'species', select=apply(spscores, 1, function(x) sqrt(x[1]^2 + x[2]^2))>1.5)

# Partition variance
varpart(Y, X, as.numeric(W))
varpart(Y.p, X, as.numeric(W))

## Each variable separately

# db-RDA
ord1_byvar = lapply(env_vars, function(x) rda(Y.p~X[,x]))
# CCA
ord1_byvar = lapply(env_vars, function(x) cca(Y~X[,x]))

names(ord1_byvar) = env_vars
lapply(ord1_byvar, anova)

# RDA: only working for sampXsp_thalli
ord2_byvar = lapply(env_vars, function(x) rda(Y.p~X[,x]+Condition(W)))
names(ord2_byvar) = env_vars

# CCA: only working for sampXsp_thalli
ord2_byvar = lapply(env_vars, function(x) cca(Y~X[,x]+Condition(W)))

names(ord2_byvar) = env_vars
lapply(ord2_byvar, anova)

## Save table of variance explained 
spord_byvar = sapply(env_vars, function(x){
	vp = varpart(Y.p, ~X[,x], ~W)

	parts = vp$part$indfract[,'Adj.R.squared']
	names(parts) = c('Env','Shared','Site','Unexplained')

	mod1 = ord1_byvar[[x]]
	mod2 = ord2_byvar[[x]]
	P = anova(mod1)[1,'Pr(>F)']
	P_Site = anova(mod2)[1,'Pr(>F)']

	c(parts, P=P, P_Site=P_Site)
})

write.csv(t(spord_byvar), './Analysis/Figures/db-RDA species of trait thalli varpart site env.csv', row.names=T)




## Plot unconstrained ordination colored by tree
ord0 = rda(Y.p)
use_year = samples[rownames(Y),'Year']
use_tree = samples[rownames(Y), 'TreeID']

# Colored by 
#col1 = rgb(0.941, 0.98, .804)
#col2 = rgb(0.341, 0.427, .039)
use_col = colorRampPalette(mycolor)(10)
color_fact = factor(use_tree)

#use_col = c('white','grey50')



par(mar=c(4,4,1,1))

use_ord = ord0
ord_sum = summary(use_ord)$cont$importance

ev = envfit(use_ord, X, choices=1:2)
pcts = paste(format(ord_sum[2,1:2]*100, digits=2), '%', sep='')

op = ordiplot(use_ord, c(1,2), type='n', xlab=paste('CA1 (',pcts[1], ')', sep=''), ylab=paste('CA2 (',pcts[2], ')', sep=''), las=1)
points(op, 'sites', pch=21, bg=use_col[factor(use_tree)], cex=0.8)
for(g in 1:10){
	ordihull(op, groups=color_fact, kind='sd', draw='polygon', col=use_col[g], show.groups=g, alpha=50)
	ordihull(op, groups=color_fact, kind='sd', draw='line',lwd=2, col=use_col[g], show.groups=g)
}
par(xpd=T)
plot(ev, labels=xvarnames_short[colnames(X)], col='grey20', add=T, cex=0.9)#arrow.mul=3,
par(xpd=F)




#####################################################################################
### Models of the distributions of individual species

library(MASS)
library(lm.beta)
library(MuMIn)
library(reshape2)

## Determine which species to model based on frequency in data set
spfreq = colSums(sampXsp_macro)
spfreq[order(spfreq, decreasing=T)]

# Focus on taxa occuring in at least 10% of samples, but exclude taxa not IDed to species
focal_taxa = names(spfreq[spfreq>0.1*nrow(sampXsp_macro)])
focal_taxa = focal_taxa[lichen_taxa[focal_taxa,'TaxonConcept']=='species']

# Loop through species
comm = sampXsp_macro
Xdata = env[rownames(comm),env_vars]

# Re-scale frequency data to be % observation
Xdata$Light_high = Xdata$Light_high*100
Xdata$Vpd_daysatfreq = Xdata$Vpd_daysatfreq*100

# Re-scale light data to kilo-Lux
Xdata$Light_mean = Xdata$Light_mean/1000

# MAY WANT TO SCALE DATA, OR, CALCULATE STD COEFS AFTER?

modlist = lapply(focal_taxa, function(i){
	yvar  = cbind(comm[,i], 1-comm[,i])
	mod = glm(yvar ~ ., data=Xdata, family=binomial(logit)) 	
})
names(modlist) = focal_taxa

modlist_single = lapply(focal_taxa, function(i){
	yvar  = cbind(comm[,i], 1-comm[,i])
	this_modlist = lapply(env_vars, function(j){
		xvar = Xdata[,j]
		mod = glm(yvar ~ xvar, family=binomial(logit))
	})
	names(this_modlist) = env_vars
	this_modlist	
})
names(modlist_single) = focal_taxa

## Calculate array of coefficients, se, ci's
modarray_single = array(NA, dim=list(length(focal_taxa), length(env_vars), 9),
	dimnames=list(focal_taxa, env_vars, c('b0','b1','se','ci.up95','ci.low95','P','b1.std','AICc','Deviance'))
)

for(i in focal_taxa){
for(j in env_vars){
	this_mod = modlist_single[[i]][[j]]
	this_coef = summary(this_mod)$coef

	modarray_single[i,j,c('b0','b1')] = this_coef[,'Estimate']
	modarray_single[i,j,'se'] = this_coef['xvar','Std. Error']
	modarray_single[i,j,c('ci.low95','ci.up95')] = confint(this_mod, 'xvar', level=0.95)
	modarray_single[i,j,'P'] = with(this_mod, pchisq(null.deviance-deviance, df=df.null-df.residual, lower.tail=F))
	modarray_single[i,j,'b1.std'] = coef(lm.beta(this_mod))['xvar']
	modarray_single[i,j,'AICc'] = AICc(this_mod)
	modarray_single[i,j,'Deviance'] = deviance(this_mod)
}}

# Make a 2-D Table of models
modmelt_single = melt(modarray_single, varnames=c('Taxon','Var','Parm'))
modtable_single = dcast(modmelt_single, Taxon+Var~...)
write.table(modtable_single, './Analysis/Figures/single species occurance prob models.txt',sep='\t', quote=F, row.names=F)

# Order table from strongest to weakest effects
modtable_single[order(modtable_single$b1.std, decreasing=T),]

# Find significant models
sigmods = subset(modtable_single, P <0.1)
sigmods_df = sigmods[,c('b1','se','ci.low95','ci.up95')]
sigmods_df = exp(sigmods_df)
sigmods_df = cbind(sigmods_df, sigmods[,c('Taxon','P','Deviance')])
sigmods_df$Predictor = xvarnames_short[sigmods$Var]
sigmods_df$n.obs = spfreq[as.character(sigmods$Taxon)]
sigmods_df = sigmods_df[,c('Taxon','Predictor','b1','se','ci.low95','ci.up95','P','Deviance','n.obs')]
write.table(format(sigmods_df, digits=3, trim=T), './Analysis/Figures/single species occurance prob models significant.txt',sep='\t', quote=F, row.names=F)


# Calculate odds ratios (odds of species occurance for every one unit increase in env variable)
mod_odds = data.frame(exp(modarray_single[,,'b1']))
format(mod_odds, scientific=F)


## Show predicted distribution changes for env variables for all taxa
mylty=c(1,2)
mycol=c('black','grey50', 'grey80')

for(j in env_vars){

pdf(paste('./Analysis/Figures/species occurance prob vs ',j,'.pdf', sep=''), height=3.5, width=3.5)
par(mar=c(4,4,1,1))
Xrange = range(Xdata[,j])
plot.new()
plot.window(xlim=Xrange, ylim=c(0,1))

for(i in focal_taxa){
	modfunc = function(x) predict(modlist_single[[i]][[j]], data.frame(xvar=x), type='response')
	
	sig = cut(modarray_single[i,j,'P'], c(0,0.05,0.1,1))
	slope = modarray_single[i,j,'b1']>0
	use_col = ifelse(i=='Usn_str', 'red', mycol[sig])
	use_lty = mylty[ifelse(as.numeric(sig)<3, 1, 2)]
	curve(modfunc, from=Xrange[1], to=Xrange[2], add=T, lwd=2, col=use_col, lty=use_lty)
	
	if(as.numeric(sig) < 3 ) text(Xrange[2], modfunc(Xrange[2]), i, adj=c(1,ifelse(slope, -.25, 1.25)))
}
axis(1)
axis(2, las=1)
mtext(xvarnames[j], 1, 3)
mtext('Occurance Probability', 2, 3)
box()
dev.off()

}

#########################################################################################
### Vertical abundance profiles for species and genera
library(MASS)

height_cat = c(6,9,12,15,36)

# Calculate relative abundance of each species in each height class
height_bin = cut(samples[use_samps,'Height_top'], height_cat)

## Models of occurance using positional variables
modlist_pos = lapply(focal_taxa, function(i){
	yvar  = cbind(comm[,i], 1-comm[,i])
	mod = glm(yvar ~ ., data=samples[rownames(comm),c('Height_top','Angle')], family=binomial(logit)) 	
})
names(modlist_pos) = focal_taxa

lapply(modlist_pos, summary) # Really only height ever has an effect

heightvar = -1*samples[rownames(comm),'Height_top']
modlist_height = lapply(focal_taxa, function(i){
	yvar  = cbind(comm[,i], 1-comm[,i])
	mod = glm(yvar ~ heightvar, family=binomial(logit)) 	
})
names(modlist_height) = focal_taxa


# Plot vertical profiles
mylty=c(1,2)
mycol=c('black','grey50', 'grey80')

Xrange = range(heightvar)
xpoints = seq(Xrange[1],Xrange[2], length.out=100)

svg('./Analysis/Figures/species height profiles num labels.svg', height=5, width=5)
par(mar=c(4,4,1,1))
plot.new()
plot.window(xlim=c(0,1), ylim=Xrange)

for(i in focal_taxa){
	this_mod = modlist_height[[i]]

	modfunc = function(x) predict(this_mod, data.frame(heightvar=x), type='response')

	P = with(this_mod, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail=F))
	sig = cut(P, c(0,0.05, 0.1, 1))
	
	#slope = modarray_single[i,j,'b1']>0
	use_col = mycol[sig]
	use_lty = mylty[ifelse(as.numeric(sig)<3, 1, 2)]

	lines(modfunc(xpoints), xpoints, lwd=2, col=use_col, lty=use_lty)

	add_on = which(focal_taxa==i)*.1
	
	#if(as.numeric(sig) < 3 ) text(modfunc(Xrange[2]-add_on), Xrange[2]-add_on, i, adj=c(-.25,0))
	if(as.numeric(sig) < 3 ) text(modfunc(Xrange[2]), Xrange[2], which(focal_taxa==i), adj=c(0.5,-0.25))


}
axis(1)
axis(2, las=1, at=seq(-7,-15,-1), labels=7:15)
mtext('Occurrence Probability', 1, 3)
mtext('Distance from Canopy Top (m)', 2, 3)
#box()
dev.off() 

# Of the species that change with height, both Pyx_sub and Pun_rud respond to an env variable (Vpd_mean)
# Check whether height distribution change still exists when controlling for env
comm = sampXsp_macro
Xdata = env[rownames(comm),env_vars]

yvar  = cbind(comm[,'Pyx_sub'], 1-comm[,'Pyx_sub'])
pyx_mod = glm(yvar ~ Vpd_mean + heightvar, data=Xdata, family=binomial(logit)) 
summary(pyx_mod)
# height still has a significant effect even with vpd - probably because VPD is independent of height
anova(pyx_mod, update(pyx_mod, .~.-heightvar), test='Chisq')

yvar  = cbind(comm[,'Pun_rud'], 1-comm[,'Pun_rud'])
pun_mod = glm(yvar ~ Vpd_daysatfreq + heightvar, data=Xdata, family=binomial(logit)) 
summary(pun_mod)
# height still has a significant effect even with vpd - probably because VPD is independent of height
anova(pun_mod, update(pun_mod, .~.-heightvar), test='Chisq')



## Plot env vars as a function of height
par(mfrow=c(1,length(env_vars)))
par(mar=c(4,4,0,0))
for(i in env_vars){
	plot(Xdata[,i], heightvar, ylab='', xlab=xvarnames[i], axes=F)
	if(i==env_vars[1]) mtext('Distance from Canopy Top (m)', 2, 3)
	axis(1)
	axis(2,las=1, at=seq(-7,-15,-1), labels=7:15)
}

cor(heightvar, Xdata)

##########################################################
### Diversity and Meta community structure

## NOT CURRENTLY USING


library(metacom)

use_comm = sampXsp_macro
use_comm = use_comm > 0

use_comm = use_comm[,colSums(use_comm)>0]

comm_ord = OrderMatrix(use_comm)

par(mar=c(1,5,6,1))
image(t(comm_ord), col=c(0,1), axes=F)
grid(ncol(comm_ord), nrow(comm_ord))
axis(3, at=(0:(ncol(comm_ord)-1))/(ncol(comm_ord)-1), labels=colnames(comm_ord), las=2)
axi

#meta = Metacommunity(use_comm)

save(meta, file='./Analysis/macro_metacomm.RData')




## Plot taxonomic richness vs environment

calc_rich = function(x, taxa){
	these_taxa = taxa[names(x[x>0]),]
	sp_gen = subset(these_taxa, TaxonConcept=='species')$Genus
	gen = subset(these_taxa, TaxonConcept=='genus')$Genus
	gen_nosp = setdiff(gen, sp_gen)

	length(sp_gen)+length(gen_nosp)
}

richness = apply(use_comm, 1, function(x) calc_rich(x, lichen_taxa))
rich_df = data.frame(richness); rich_df$SampID=rownames(rich_df)
rich_df = merge(rich_df, env[,c('SampID',env_vars)])


## Model richness as a function of environment

richmods = as.list(rep(NA, length(env_vars)))
names(richmods) = env_vars

for(i in env_vars){
	richmods[[i]] = glm(richness~rich_df[,i], data=rich_df, family=poisson(link='log'))
}

# Plot models
use_lty=c(1,2)
use_col=c('black','grey50')

pdf('./Analysis/Figures/species richness vs env all taxa.pdf', height=4, width=12)
par(mfrow=c(1,5))
for(i in env_vars){
	xvar = rich_df[,i]
	this_mod = richmods[[i]]
	
	plot(rich_df$richness ~ xvar, ylab='Species richness', xlab=xvarnames[i], las=1, ylim=c(0,18))
	
	sig = coef(summary(this_mod))[2,4] < 0.05

	modfunc = function(x) exp(coef(this_mod)%*%rbind(1,x))

	curve(modfunc(x), from=min(xvar), to=max(xvar), lwd=2, 
		lty=use_lty[ifelse(sig, 1, 2)], col=use_col[ifelse(sig, 1, 2)], add=T)

}
dev.off()

par(mfrow=c())
plot(richness~Light_mean_sum, data=env)
plot(richness~Temp_mean_sum, data=env)

mod_rich = glm(richness~Light_mean_sum, data=env, family=poisson)
summary(mod_rich)

pdf('./Analysis/Figures/species richness histogram.pdf', height=4, width=4)
par(mar=c(4,4,1,1))
plot(table(env$richness), las=1, ylab='Num. Samples', xlab='Num. Species', ylim=c(0,20), lwd=4, lend=1)
dev.off()


