## This script models canopy lichen functional traits
options(stringsAsFactors=F)

working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/'
sql_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Data/SQLite Tables/'
data_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Data/Derived Tables/'
script_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Analysis/GitHub/NC-canopy-lichens/'

setwd(working_dir)
options(stringsAsFactors=F)

# Read in functions
source(paste(script_dir, 'CFT_functions.R', sep=''))

# Set working colors
mycolor = read.csv('../blue2red_10colramp.txt')[9:2,]
mycolor = rgb(mycolor, maxColorValue=255)

# Define predictors
env_vars = c('Light_mean','Light_high','Temp_max','Vpd_mean','Vpd_daysatfreq')
xvarnames = data.frame(var = env_vars)
xvarnames$displayName = c('Mean Light Intensity (Lux)', 'Freq. Full Sun', 
	'Max. Temperature (C)','Mean VPD','Freq. Daytime Saturated Air')
rownames(xvarnames) = xvarnames$var

######################################################
### Read in Data

traits = read.csv(paste(sql_dir, 'derived_traits.csv', sep=''))
env = read.csv(paste(data_dir, 'loggerdata.csv', sep=''), row.names=1)
samples = read.csv(paste(sql_dir, 'samples.csv', sep=''))

taxa = read.csv(paste(sql_dir, 'taxa.csv', sep=''))
rownames(taxa) = taxa$TaxonID

# Merge sample and env data
env = merge(env, samples, by='SampID')
rownames(env) = paste('S', env$SampID, sep='')

# Add column for Genus to trait data
traits$Genus = sapply(traits$TaxonID, function(x) taxa[x, 'Genus'])

# Read in table of measured traits and their display names
traitnames = read.csv('trait_types.csv', row.names=1)

# Subset data to remove unidentified thalli
traits = subset(traits, is.na(TaxonID)|!(TaxonID %in% c('Fol','')))

# Merge traits and env data
use_data = merge(traits, env)

#####################################################
### Clean up trait data 

## Examine outliers
table(traits$Sexual_abun)
table(traits$Asexual_abun)
table(traits$Rhizine_abun)

numeric_traits = rownames(subset(traitnames, mode=='N'))

pdf('./Analysis/Figures/numeric trait outlier.pdf', height=8, width=8)
par(mfrow=c(3,3))
par(mar=c(3,3,1,1))
for(i in numeric_traits){
	plot(traits[,i]~ I(rank(traits[,i])/max(rank(traits[,i]))), xlab='', ylab='', main=i)
	abline(v=c(0.05, 0.95), col=2)
}
dev.off()

#  Omit thallus 22 b/c very small and lost, which made an outlier
use_data = subset(use_data, ThallusID != 22)

## Boxplots of traits across genera

use_traits = c('Water_capacity','STA','Tot_chl_DW','Chla2b','Rhizine_length')

pdf('./Analysis/Figures/trait distributions by genus.pdf', height=5, width=6)
par(mar=c(7,3,3,1))
for(i in use_traits){
	trait_means = tapply(use_data[,i], use_data$Genus, FUN=mean, na.rm=T)
	boxplot(use_data[,i]~use_data$Genus, las=2, main=traitnames[i,'displayName'])
}
dev.off()

## Assign NA to trait values that appear to be outlier / errors based on previous plot
## Decided that high chlorophyll not outliers b/c few Physcia samples
use_data[use_data$ThallusID %in% c(419,501,521,541), 'STA'] = NA # Values abnormally high- area measured prob doesn't correspond to mass
use_data[use_data$ThallusID==8, 'Water_capacity'] = NA # Can't have negative water capacity
use_data[(use_data$Chla2b<0) & (!is.na(use_data$Chla2b)),'Chla2b'] = NA # Can't have negative ratios. Measurements probably erroneous

#####################################################
### Distribution of Genera/Species across samples / environments

taxon_freq = table(use_data$TaxonID)
taxon_freq[order(taxon_freq)]

genus_freq = table(use_data$Genus)
genus_freq[order(genus_freq)]

## Summary table of genera
genus_tab = data.frame(genus_freq)
names(genus_tab) = c('Genus','Thalli')

genus_samp = data.frame(colSums(table(use_data$SampID, use_data$Genus)>0))
names(genus_samp) = 'Samples'; genus_samp$Genus=rownames(genus_samp)

genus_tree = data.frame(colSums(table(use_data$TreeID, use_data$Genus)>0))
names(genus_tree) = 'Trees'; genus_tree$Genus=rownames(genus_tree)

genus_tab = merge(genus_tab, genus_samp)
genus_tab = merge(genus_tab, genus_tree)

write.csv(genus_tab, './Analysis/genus_frequency.csv', row.names=F)

# Number of species in each genus
taxon_tab = data.frame(taxon_freq); names(taxon_tab) = c('TaxonID','Thalli')
merge(taxon_tab, unique(use_data[,c('Genus','TaxonID')]))


## Plot taxonomic richness vs light and temp

sampXsp = xtabs(~SampID+TaxonID, data=use_data)

calc_rich = function(x, taxa){
	these_taxa = taxa[names(x[x>0]),]
	sp_gen = subset(these_taxa, TaxonConcept=='species')$Genus
	gen = subset(these_taxa, TaxonConcept=='genus')$Genus
	gen_nosp = setdiff(gen, sp_gen)

	length(sp_gen)+length(gen_nosp)
}

richness = apply(sampXsp, 1, function(x) calc_rich(x, taxa))
richness = data.frame(richness); richness$SampID=rownames(richness)

env = merge(env, richness)

plot(richness~Light_mean_sum, data=env)
plot(richness~Temp_mean_sum, data=env)

mod_rich = glm(richness~Light_mean_sum, data=env, family=poisson)
summary(mod_rich)

pdf('./Analysis/Figures/species richness histogram.pdf', height=4, width=4)
par(mar=c(4,4,1,1))
plot(table(env$richness), las=1, ylab='Num. Samples', xlab='Num. Species', ylim=c(0,20), lwd=4, lend=1)
dev.off()

## Correlations among traits
cor(model_data[,use_traits], use='pairwise.complete.obs')

## Ordination of individuals
library(vegan)

noNA_data = na.omit(model_data)
traitord = rda(noNA_data[,use_traits], center=T, scale=T)

mycol=rainbow(11)

plot(traitord, type = "n")
points(traitord, display = "sites", cex = 0.8, pch=1, col=mycol[factor(noNA_data$Genus)])
text(traitord, display = "spec", cex=0.7, col="blue")



######################################################
## PCA of env

env$BranchPos = factor(env$BranchPos, levels=c('low','mid','high'))

# Plot correlations among env variables
plot(env$Light_mean~env$Light_high)
plot(env$Light_mean~env$Light_p90)
plot(env$Vpd_mean~env$Vpd_satfreq)
plot(env$Vpd_mean~env$Vpd_daysatfreq)

# PCA of chosen variables
use_vars = c('Light_mean','Light_high','Temp_max','Vpd_mean','Vpd_daysatfreq')

envpca = prcomp(env[,use_vars], scale=T)
summary(envpca);envpca


## Summarize env by position

boxplot(Temp_mean_sum~BranchPos, data=env, las=1)
boxplot(Light_mean_sum~BranchPos, data=env, las=1)
boxplot(Light_mean_sum~Pair, data=env, las=1)
boxplot(Temp_mean_sum~Pair, data=env, las=1)

mycol = c('blue','red')

colorby = factor(env$Pair)

pdf('./Analysis/Figures/sample summer light temp by height.pdf', height=6, width=4)
par(mfrow=c(2,1))
par(mar=c(1,6,5,1))
plot(Temp_mean_sum~Height, data=env, col=mycol[colorby], pch=16, axes=F, ylab='')
axis(2, las=1); box()
legend('top', levels(colorby), col=mycol, pch=16, inset=-.25, ncol=2, xpd=T, bty='n') 
mtext('Mean Temp. (May-Sept)', 2, 4)
par(mar=c(5,6,1,1))
plot(Light_mean_sum~Height, data=env, col=mycol[colorby], pch=16, las=1, ylab='')
mtext('Mean Light (May-Sept)', 2, 4)
dev.off()


######################################################
### Linear Models of Individual Traits

library(lme4)
library(LMERConvenienceFunctions)

# Make SampID, TreeID, and Genus a factor
use_data$SampID = factor(use_data$SampID)
use_data$TreeID = factor(use_data$TreeID)
use_data$Genus = factor(use_data$Genus)

# Define functional traits
functraits = c('STA','Cortex_thickness','Thallus_thickness','Tot_chl_DW','Chla2b',
	'Sexual_abun','Asexual_abun','Water_capacity','Rhizine_length','Rhizine_abun')
numtraits = rownames(subset(traitnames, mode=='N'))
cattraits = rownames(subset(traitnames, mode=='O'))


### Two-way ANOVA of each trait with sample and Genus as predictors


# TreeID, SampID and Genus are random effects with SampID nested in TreeID
trait_aov = lapply(functraits[functraits%in%numtraits], function(i){

	anovamat = matrix(NA, nrow=3, ncol=5, dimnames=list(c('Sample','Tree','Genus'), c('Variance','ChiFromNull','PFromNull','ChiFromFull','PFromFull')))
	
	modfull = lmer(use_data[,i] ~ 1+(1|Genus)+(1|TreeID)+(1|TreeID:SampID), data=use_data, REML=F)
	modtree = lmer(use_data[,i] ~ 1+(1|TreeID), data=use_data, REML=F)
	modsamp = lmer(use_data[,i] ~ 1+(1|SampID), data=use_data, REML=F)
	modgen = lmer(use_data[,i] ~ 1+(1|Genus), data=use_data, REML=F)
	modnotree = lmer(use_data[,i] ~ 1+(1|Genus)+(1|TreeID:SampID), data=use_data, REML=F)
	modnosamp = lmer(use_data[,i] ~ 1+(1|Genus)+(1|TreeID), data=use_data, REML=F)
	modnogen = lmer(use_data[,i] ~ 1+(1|TreeID)+(1|TreeID:SampID), data=use_data, REML=F)
	modnull = lm(use_data[,i] ~ 1, data=use_data)

	anovamat['Sample', 'ChiFromNull'] = as.numeric(2*(logLik(modsamp) - logLik(modnull)))
	anovamat['Sample','PFromNull'] = pchisq(anovamat['Sample', 'ChiFromNull'], 1, lower=F)
	anovamat['Tree', 'ChiFromNull'] = as.numeric(2*(logLik(modtree) - logLik(modnull)))
	anovamat['Tree','PFromNull'] = pchisq(anovamat['Tree', 'ChiFromNull'], 1, lower=F)
	anovamat['Genus', 'ChiFromNull'] = as.numeric(2*(logLik(modgen) - logLik(modnull)))
	anovamat['Genus','PFromNull'] = pchisq(anovamat['Genus', 'ChiFromNull'], 1, lower=F)
	anovamat['Sample', 'ChiFromFull'] = as.numeric(2*(logLik(modfull) - logLik(modnosamp)))
	anovamat['Sample','PFromFull'] = pchisq(anovamat['Sample', 'ChiFromFull'], 1, lower=F)
	anovamat['Tree', 'ChiFromFull'] = as.numeric(2*(logLik(modfull) - logLik(modnotree)))
	anovamat['Tree','PFromFull'] = pchisq(anovamat['Tree', 'ChiFromFull'], 1, lower=F)
	anovamat['Genus', 'ChiFromFull'] = as.numeric(2*(logLik(modfull) - logLik(modnogen)))
	anovamat['Genus','PFromFull'] = pchisq(anovamat['Genus', 'ChiFromFull'], 1, lower=F)

	anovamat['Sample','Variance'] = as.numeric(VarCorr(modfull)$'TreeID:SampID')
	anovamat['Tree','Variance'] = as.numeric(VarCorr(modfull)$TreeID)
	anovamat['Genus','Variance'] = as.numeric(VarCorr(modfull)$Genus)

	residVar = attr(VarCorr(modfull),'sc')

	list(residVar=residVar, anovamat=anovamat)
})

names(trait_aov) = functraits[functraits%in%numtraits]


## Conclusions:
# STA : genus and sample explain significant variation, but lots of residual and genus more so
# Cortex_thickness : genus and sample explain, but more so genus
# Thallus_thickness : only genus explains significant variation
# Tot_chl_DW : genus explains some variation and sample might, but lots of residual variation
# Chla2b : models didn't work?
# Water_capacity : tree and genus explain significant variation, but more so tree and there is lots of residual variance
# Rhizine_length : all variance explained by genus


# Make a data frame of data for modeling
use_traits = c('Water_capacity','Thallus_thickness','STA','Cortex_thickness','Rhizine_length','Tot_chl_DW','Chla2b')
model_data = use_data[,c(use_traits, env_vars, 'Genus','SampID')]

# Transform responses
# Effectively we'll be fitting log-normal models since all responses are constrained to be positive
# Cortex thickness can be 0 so we add 1
model_data$Cortex_thickness = model_data$Cortex_thickness + 1

layout(matrix(1:(2*length(use_traits)), nrow=2, byrow=F))
for(y in use_traits){
	hist(model_data[,y], main=y)
	hist(log(model_data[,y]), main=paste('log',y))
}

for(y in use_traits) model_data[,y] = log(model_data[,y])

# Re-scale predictors
model_data$Light_mean = model_data$Light_mean/10000
model_data$Temp_max = model_data$Temp_max/10

# Create objects for storing models
modlist = vector('list', length(use_traits)*length(env_vars))
dim(modlist) = c(length(use_traits),length(env_vars))
dimnames(modlist) = list(use_traits, env_vars)

parm_ests = array(NA, dim=c(length(use_traits), length(env_vars), 5, 3), 
	dimnames = list(use_traits, env_vars, c('sigma.samp','sigma.genus','sigma.res','b0','b1'), c('est','low95','up95')))

genus_ests = array(NA, dim=c(length(use_traits), length(env_vars), 11, 3), 
	dimnames = list(use_traits, env_vars, levels(factor(model_data$Genus)), c('est','low95','up95')))


## Loop through all traits and all env vars
# Note that number of observations and levels(SampID) will differ across models
for(i in use_traits){
for(j in env_vars){
	mod = lmer(model_data[,i] ~ 1 + model_data[,j] + (1|Genus) + (1|SampID), data=model_data, REML=T)
	modlist[i,j][[1]] = mod
	
	ests = c(data.frame(VarCorr(mod))[,'sdcor'], fixef(mod))
	g_ests = ranef(mod, whichel='Genus')[[1]]$'(Intercept)'
	ints = confint(mod, parm=1:5, method='profile')
	g_sds = sqrt(as.numeric(attr(ranef(mod, whichel='Genus', condVar=T)[[1]], 'postVar')))
	g_ints = g_ests + g_sds%*%t(c(-1.96,1.96))
	
	parm_ests[i,j,,'est'] = ests
	parm_ests[i,j,,c('low95','up95')] = ints
	genus_ests[i,j,,'est'] = g_ests
	genus_ests[i,j,,c('low95','up95')] = g_ints
}}

save(modlist, parm_ests, genus_ests, file='./Analysis/REML single variable FT models.RData')

## Plot estimated fixed effects
library(reshape)

plot_data = melt(parm_ests[,,'b1','est'])
low95_data = melt(parm_ests[,,'b1','est'])

names(plot_data) = c('yvar','xvar','est')
xyplot(xvar~est|yvar, data=plot_data, panel=function(x,y,subscripts,...){
	panel.abline(v=0, col='grey50', lty=2)
	panel.points(x,y, pch=16, col=1)
	#panel.segments(parm_ests[1,y,,'low95'],y,parm_ests[1,y,,'up95'],y)
})

xvar_levels = factor(env_vars)

pdf('./Analysis/Figures/REML single variable model effects zoomed.pdf', height=9, width=4)
layout(matrix(1:length(use_traits), ncol=1))
par(mar=c(2,13,2,.5))
for(i in use_traits){
	this_data = parm_ests[i,,'b1',]
	xrange = c(-2,2)#range(this_data)
	plot(as.numeric(xvar_levels)~this_data[,'est'], xlim=xrange, axes=F, 
		xlab='', ylab='', lwd=1)
	abline(v=0, col='grey50', lty=2)
	segments(this_data[, 'low95'], as.numeric(xvar_levels), 
		this_data[, 'up95'], as.numeric(xvar_levels), lwd=1)
	axis(1)
	axis(2, at=1:length(env_vars), labels=xvarnames[xvar_levels,'displayName'], las=1)
	mtext(traitnames[i,'displayName'], 3, 0)
	box()
}
dev.off()



###########################################
### Fit Baysian Models


library('R2OpenBUGS')
library('coda')
library(MASS) # ginv for calculating inverse matrix

setwd('./Analysis/BUGS')

## See script BUGS_models.R for fitting of all models.

## Fit one model:

# Set trait and find missing observations
i = 'Cortex_thickness'
keeprows = !is.na(use_data[,i])
j = 'Vpd_daysatfreq'

# Set model variables
y = use_data[keeprows,i] + 1# Some cortices have 0 thickness
n = length(y)
samp = as.numeric(factor(use_data[keeprows,'SampID']))
J = length(unique(samp))
genus = as.numeric(factor(use_data[keeprows,'Genus']))
G = length(unique(genus))

# Set model parameters
mod_data = list('y','samp','genus','n','J','G','xvar')
mod_init = function(){list(a=rnorm(J), b=rnorm(G), g0=rnorm(1), g1=rnorm(1), 
	sigma.y=runif(1), sigma.a=runif(1), sigma.b=runif(1))}
mod_parms = c('a','b','g0','g1','sigma.y','sigma.a','sigma.b')
mod_file = paste('bayes_mod_',i,'.txt', sep='')

# Run model
modbug = bugs(mod_data, mod_init, mod_parms, mod_file, 
	n.chains=3, n.iter=500, codaPkg=T, n.burnin=0, n.thin=100,
	saveExec=T, working.directory='./Analysis/BUGS', restart=F)
modbug = read.bugs(modbug)

# Thin results
thinned = window(modbug, thin=1, start=100)

# Plot results of each chain
mylayout=c(2,5,8)
par(ask=T)
xyplot(thinned, layout=mylayout, scales=list(rot=0))
densityplot(thinned, layout=mylayout, scales=list(rot=0), aspect='fill')
acfplot(thinned, layout=mylayout, scales=list(rot=0), aspect='fill')
gelman.diag(thinned, autoburnin=F)



### Load all models and plot probability smear graphs ###
use_xvars = c('Light_mean','Light_high','Temp_max','Vpd_mean','Vpd_daysatfreq')

# Make empty list for storing univariate models
mod_list = list(list(NA,NA,NA,NA,NA),list(NA,NA,NA,NA,NA),list(NA,NA,NA,NA,NA),list(NA,NA,NA,NA,NA),list(NA,NA,NA,NA,NA),list(NA,NA,NA,NA,NA))
names(mod_list) = c('Tot_chl_DW','Chla2b','Water_capacity','Cortex_thickness','Thallus_thickness','Rhizine_length')
for(i in 1:length(mod_list)) names(mod_list[[i]]) = use_xvars

# Read in and stored thinned posterior distributions
for(i in names(mod_list)){
for(j in use_xvars){
	load(paste('./Analysis/BUGS/',i,'/',j,'/',i,'-',j,'-model.RData',sep=''))
	thinned = window(modbug, start=5001, thin=5) # original sampler thinned every 100, so total is every 500
	mod_list[[i]][[j]] = thinned
}}

# Check autocorrelation - do this manually for each model
mylayout=c(3,4)
#par(ask=T)
pdf('./Analysis/Figures/model autocorrelation plots.pdf', height=8, width=10.5)
for(i in names(mod_list)){
for(j in use_xvars){
print(acfplot(mod_list[[i]][[j]], main = paste(i,j, sep=' vs. '), layout=mylayout, scales=list(rot=0), aspect='fill'))
}}
dev.off()

# Check convergence
sink('./Analysis/Figures/Model convergence.txt')
for(i in names(mod_list)){
for(j in use_xvars){
print(paste(i,j,sep=' vs. '))
print(gelman.diag(mod_list[[i]][[j]], autoburnin=F))
}}
sink()

## Thin chains based on observed autocorrelation
thin_ints = read.csv('./Analysis/BUGS/thin_intervals.csv', row.names=1)

mod_thin = mod_list
for(i in names(mod_list)){
for(j in use_xvars){
	thinned = window(mod_list[[i]][[j]], start=1, thin=thin_ints[i,j])
	mod_thin[[i]][[j]] = thinned
}}


# Create list of pooled chains
pooled_list = mod_thin

for(i in names(mod_list)){
for(j in names(mod_list[[1]])){
	chains = mod_list[[i]][[j]]
	c1 = as.matrix(chains[[1]])
	c2 = as.matrix(chains[[2]])
	c3 = as.matrix(chains[[3]])
	comb = rbind(c1,c2,c3)
	comb = as.mcmc(comb)
	
	pooled_list[[i]][[j]] = comb
}}


## Plot slopes


# Define response variable
i = 'Tot_chl_DW'

## Slope estimates
# Define plot limits
for(i in names(pooled_list)){

Xlims = sapply(use_xvars, function(j) max(abs(range(pooled_list[[i]][[j]][,'g1']))))

svg(paste('./Analysis/Figures/',i,' slope parameter posteriors.svg', sep=''), height=6, width=4.5)
layout(matrix(1:5, ncol=1, nrow=5))
par(mar=c(4,1,0.5,1))
for(j in use_xvars){

	post = pooled_list[[i]][[j]][,'g1']	
	d = density(post)
	g1int = HPDinterval(post, 0.95)
	g1int50 = HPDinterval(post, 0.50)	
	plot(d, axes=F, lwd=2, xlab='', ylab='', main='', xlim=c(-Xlims[j],Xlims[j]))
	
	lines(c(g1int[1],g1int[2]),c(0,0), lwd=3, col='black')
	lines(c(g1int50[1],g1int50[2]),c(0,0), lwd=4, col='grey50')
	points(median(post),0, cex=2, pch='|', lwd=2)

	#par(xpd=T)
	abline(v=0, lty=1, col=mycolor[7], lwd=2)
	#par(xpd=F)
	
	mtext(xvarnames[j,'displayName'], 3, -1, adj=0)
	axis(1, line=.5, cex.lab=2)
	
}
dev.off()
}


## Variance Components
# inspect and use to set limit
i = 'Cortex_thickness'
j = 'Vpd_mean'

Xlims = sapply(use_xvars, function(j) max(pooled_list[[i]][[j]][,c('sigma.a','sigma.b','sigma.y')]))
Xlim = c(0,1)

use_cols = c(samp='orangered1',genus='royalblue1',resid='black')

svg(paste('./Analysis/Figures/',i,' variance components posteriors.svg', sep=''), height=5, width=10.5)
layout(matrix(1:15, ncol=5, nrow=3))
for(j in use_xvars){

	post.samp = (pooled_list[[i]][[j]][,'sigma.a'])#^2 # to get variances from stddev	 
	post.genus = (pooled_list[[i]][[j]][,'sigma.b'])#^2
	post.resid = (pooled_list[[i]][[j]][,'sigma.y'])#^2
	
	d.samp = density(post.samp)
	d.genus = density(post.genus)
	d.resid = density(post.resid)

	samp.int = HPDinterval(post.samp, 0.95)
	genus.int = HPDinterval(post.genus, 0.95)
	resid.int = HPDinterval(post.resid, 0.95)
	
	par(mar=c(0,.5,3,.5))
	
	plot(d.samp, axes=F, lwd=2, xlab='', ylab='', main='', xlim=Xlim, col=use_cols['samp'])
	lines(c(samp.int[1],samp.int[2]),c(0,0), lwd=3, col=use_cols['samp'])
	points(median(post.samp), 0, cex=2, pch='|', lwd=2, col=use_cols['samp'])
	mtext(xvarnames[j,'displayName'], 3, 0)
	if(j==use_xvars[1]) mtext('Sample', 2, 0, col=use_cols['samp'])	
	
	par(mar=c(1.5,.5,1.5,.5))
	plot(d.genus, axes=F, lwd=2, xlab='', ylab='', main='', xlim=Xlim, col=use_cols['genus'])
	lines(c(genus.int[1],genus.int[2]),c(0,0), lwd=3, col=use_cols['genus'])
	points(median(post.genus), 0, cex=2, pch='|', lwd=2, col=use_cols['genus'])
	if(j==use_xvars[1]) mtext('Genus', 2, 0, col=use_cols['genus'])
	
	par(mar=c(3,.5,0,.5))

	plot(d.resid,axes=F, lwd=2, xlab='', ylab='', main='', xlim=Xlim, col=use_cols['resid'])
	lines(c(resid.int[1],resid.int[2]),c(0,0), lwd=3, col=use_cols['resid'])
	points(median(post.resid), 0, cex=2, pch='|', lwd=2, col=use_cols['resid'])
	if(j==use_xvars[1]) mtext('Residual', 2, 0, col=use_cols['resid'])

	axis(1, line=.5)
	axis(1, line=.5, at=median(post.samp), labels=F, col=use_cols['samp'], lwd=2)
	axis(1, line=.5, at=median(post.genus), labels=F, col=use_cols['genus'], lwd=2)
	axis(1, line=.5, at=median(post.resid), labels=F, col=use_cols['resid'], lwd=2)
}
dev.off()



## Plot prediction for Tot_chl_DW ~ Light_mean with separate curves for each genus

i = 'Cortex_thickness' #'Tot_chl_DW'
j = 'VPD_mean' # 'Light_mean'


use_mod = mod_thin[[i]][[j]]
mod_sum = summary(use_mod)
means = mod_sum$statistics[,'Mean']

use_mcmc = pooled_list[[i]][[j]]

# Prediction based on parameter means
modfunc = function(x, b) exp( means[paste('b[',b,']',sep='')] + means['g0'] + means['g1']*x )

# Prediction based on MCMC samples
predict_func = function(x, b) use_mcmc[,paste('b[',b,']',sep='')] + use_mcmc[,'g0'] + use_mcmc[,'g1']*x

use_x = seq(min(use_data[,j]), max(use_data[,j]), length.out=50)

pdf('./Analysis/Figures/Tot_chl_DW vs Light_mean.pdf', height=4, width=4)
par(mar=c(5,5,1,1))
plot(Tot_chl_DW~Light_mean, data=use_data, las=1, xlab=xvarnames['Light_mean','displayName'], 
	ylab='Chlorophyll Conc. (umol/g)')
for(b in 1:8){
	use_y = sapply(use_x, function(x) do.call(predict_func, list(x=x,b=b)))
	mean_y = colMeans(use_y)
	lines(use_x, exp(mean_y), col=b)
}
dev.off()

pdf('./Analysis/Figures/Tot_chl_DW vs Light_high.pdf', height=4, width=4)
par(mar=c(5,5,1,1))
plot(Tot_chl_DW~Light_high, data=use_data, las=1, xlab=xvarnames['Light_high','displayName'], 
	ylab='Chlorophyll Conc. (umol/g)')
for(b in 1:8){
	use_y = sapply(use_x, function(x) do.call(predict_func, list(x=x,b=b)))
	mean_y = colMeans(use_y)
	lines(use_x, exp(mean_y), col=b)
}
dev.off()

# Plot for all env vars
for(j in use_xvars){

# Prediction based on MCMC samples
use_mcmc = pooled_list[[i]][[j]]
predict_func = function(x, b) use_mcmc[,paste('b[',b,']',sep='')] + use_mcmc[,'g0'] + use_mcmc[,'g1']*x
use_x = seq(min(use_data[,j]), max(use_data[,j]), length.out=50)

pdf(paste('./Analysis/Figures/', i, ' vs ', j,'.pdf', sep=''), height=4, width=4)
par(mar=c(5,5,1,1))
plot(use_data[,i]~use_data[,j], las=1, xlab=xvarnames[j,'displayName'], 
	ylab=traitnames[i,'displayName'])
for(b in 1:8){
	use_y = sapply(use_x, function(x) do.call(predict_func, list(x=x,b=b)))
	mean_y = colMeans(use_y)
	lines(use_x, exp(mean_y-1), col=b) # Substract one for cortex thickness only due to initialdata transformation
}
dev.off()
}


## NOT DONE YET

svg('./Analysis/Figures/Tot_chl_DW vs Light_mean_sum.svg', height=4, width=4)
par(mar=c(4,4.5,.5,.5))
plot(log(Tot_chl_DW)~Light_mean_sum, data=use_data, las=1, lwd=2,cex=1.5,
	ylab=expression(paste('Log Chlorophyll (',mu,'mol / g DW)', sep='')), 
	xlab='Mean Light Intensity (lumens)'
)
mod = lm(log(Tot_chl_DW)~Light_mean_sum, data=use_data)
modfunc = function(x) {coef(mod)[1] + coef(mod)[2]*x}
curve(modfunc, min(use_data$Light_mean_sum), max(use_data$Light_mean_sum), 
	add=T, lwd=3, col=mycolor[7])
dev.off()

svg('./Analysis/Figures/Tot_chl_DW vs Light_mean_sum.svg', height=4, width=4)
par(mar=c(4,4.5,.5,.5))
plot(log(Tot_chl_DW)~Light_mean_sum, data=use_data, las=1, lwd=2,cex=1.5,
	ylab=expression(paste('Log Chlorophyll (',mu,'mol / g DW)', sep='')), 
	xlab='Mean Light Intensity (lumens)'
)
mod = lm(log(Tot_chl_DW)~Light_mean_sum, data=use_data)
modfunc = function(x) {coef(mod)[1] + coef(mod)[2]*x}
curve(modfunc, min(use_data$Light_mean_sum), max(use_data$Light_mean_sum), 
	add=T, lwd=3, col=mycolor[7])
dev.off()




for(a in 1:18){
for(b in 1:8){
	
	use_y = do.call(modsumfunc,list(a=a,b=b))
	points(samples[paste('S',a,sep=''),'Temp_mean_sum'],use_y, col=b, pch=0, lwd=2)
	
}
	points(samples[paste('S',a,sep=''),'Temp_mean_sum'],means[paste('a[',a,']',sep='')], col=1, pch=16, cex=1.5)
}

par(xpd=T)
points(rep(par('usr')[1],18), means[sapply(1:18, function(a) paste('a[',a,']',sep=''))])
points(par('usr')[1], means['g0']+(par('usr')[1]*means['g1']), pch=4, col=2)
points(par('usr')[1], means['g0']+means['g1'], pch=4, col=2)

mean(means[sapply(1:18, function(a) paste('a[',a,']',sep=''))])

## Visualizing distributions
use_x = seq(-20,20,.1)

plot(use_x,dgamma(use_x, shape=0.001, scale=0.001), type='l')

plot(use_x, dnorm(use_x, 0, 3.6^2), type='l', ylim=c(0,.1))
points(y,rep(0,length(y)))


#####################################################
### Testing fake data

make_fake_data = function(x, samp, genus, g0, beta, sigma.y, sigma.a, sigma.b, do.log){
	
	G = length(unique(genus))
	J = length(unique(samp))
	n = length(x)

	genus = as.numeric(factor(genus))
	samp = as.numeric(factor(samp))

	b.g = g0 + rnorm(G, 0, sigma.b)
	a.j = rnorm(J, 0, sigma.a)	

	if(length(beta)>1){
		y = x%*%beta
	} else {
		y = x*beta
	}
		
	y  = y + b.g[genus] + a.j[samp] + rnorm(n, 0, sigma.y)

	if(do.log) y = exp(y)

	list(y = y, b.g = b.g, a.j = a.j)

}


fake_data = use_data[,c('SampID','Genus')]
fake_data$xvar = use_data$Vpd_mean
fake_data$y1 = make_fake_data(fake_data$xvar, fake_data$SampID, fake_data$Genus, 
	g0 = 0.5, beta = .5, sigma.y = 1, sigma.a = 1, sigma.b = 1, do.log=T)


samp = as.numeric(factor(fake_data$SampID))
J = length(unique(samp))
genus = as.numeric(factor(fake_data$Genus))
G = length(unique(genus))
xvar = fake_data$xvar

try_parms = expand.grid(g0 = .5, beta = c(0, .5, 1), sigma.y = 1, sigma.a = c(.5, 1, 2), sigma.b = c(.5, 1, 2))

mod_list = as.list(rep(NA, nrow(try_parms)))
Y = matrix(NA, nrow=nrow(fake_data), ncol=nrow(try_parms))
B = matrix(NA, nrow=G, ncol=nrow(try_parms))
A = matrix(NA, nrow=J, ncol=nrow(try_parms))


for(i in 1:nrow(try_parms)){
	sim_data = make_fake_data(fake_data$xvar, fake_data$SampID, fake_data$Genus, 
		g0 = try_parms$g0[i], beta = try_parms$beta[i], sigma.y = try_parms$sigma.y[i], 
		sigma.a = try_parms$sigma.a[i], sigma.b = try_parms$sigma.b[i], do.log=T)
	Y[,i] = sim_data$y
	B[,i] = sim_data$b.g
	A[,i] = sim_data$a.j

	# Set model variables
	y = Y[,i] 
	n = length(y)

	# Set model parameters
	mod_data = list('y','samp','genus','n','J','G','xvar')
	mod_init = function(){list(a=rnorm(J), b=rnorm(G), g0=rnorm(1), g1=rnorm(1), 
		sigma.y=runif(1), sigma.a=runif(1), sigma.b=runif(1))}
	mod_parms = c('a','b','g0','g1','sigma.y','sigma.a','sigma.b')
	mod_file = 'bayes_mod_Cortex_thickness.txt'

	# Run model
	modbug = bugs(mod_data, mod_init, mod_parms, mod_file, 
		n.chains=3, n.iter=600, codaPkg=T, n.burnin=100, n.thin=100,
		saveExec=T, working.directory='./Analysis/BUGS', restart=F)
	mod_list[[i]] = read.bugs(modbug)
}

save(mod_list, Y, A, B, try_parms, file='./Analysis/BUGS/simulated data and models.RData')

# Get estimated parameters
parm_ests = array(NA, dim=c(nrow(try_parms), ncol(mod_list[[1]][[1]]), 3), dimnames=list(1:nrow(try_parms), colnames(mod_list[[1]][[1]]), c('med','low95','up95')))
for(i in 1:nrow(try_parms)){
	this_mod = mod_list[[i]]

	# Pool all chains
	c1 = as.matrix(this_mod[[1]])
	c2 = as.matrix(this_mod[[2]])
	c3 = as.matrix(this_mod[[3]])
	pooled = rbind(c1,c2,c3)
	pooled = as.mcmc(pooled)
	
	# Get median
	parm_ests[i,,1] = apply(pooled[,dimnames(parm_ests)[[2]]], 2, median)

	# Get 95% intervals
	parm_ests[i,,2:3] = HPDinterval(pooled, .95)[dimnames(parm_ests)[[2]],]
}

par(ask=T)

# Load simulation that was run on cluster
load('./Analysis/BUGS/simulated data and models.RData')


pdf('./Analysis/Figures/simulated data model performance coefs.pdf', height=15, width=15)
par(mfrow=c(3,3))
par(mar=c(2,2,.5,.5))
for(i in 1:nrow(try_parms)){
	plot(log(Y[,i])~xvar, ylab='', xlab='')

	parm_text = paste('g0=',try_parms[i,1],' g1=',try_parms[i,2],' sigma.y=',try_parms[i,3],' sigma.a=',try_parms[i,4],' sigma.b=',try_parms[i,5], sep='') 
	use_ests = format(parm_ests[i,c('g0','g1','sigma.y','sigma.a','sigma.b'),'med'], digits=2)
	est_text = paste('g0=',use_ests[1],' g1=',use_ests[2],' sigma.y=',use_ests[3],' sigma.a=',use_ests[4],' sigma.b=',use_ests[5], sep='')
	mtext(parm_text, 3, -1)
	mtext(est_text, 3, -2, col=2)	

	abline(lm(log(Y[,i])~xvar), lwd=3, , col='grey30', lty=3) #Line from simple linear model	
	abline(try_parms[i,1], try_parms[i,2], lwd=3) # Line for simulated data
	abline(parm_ests[i,'g0','med'], parm_ests[i,'g1','med'], lwd=3, col=2)
}
dev.off()

# Plot Genus-level means
pdf('./Analysis/Figures/simulated data model performance genus coefs.pdf', height=9, width=9)
par(mfrow=c(3,3))
par(mar=c(2,2,.5,.5))
for(i in 1:nrow(try_parms)){
	use_b = B[,i]
	b_ests = parm_ests[i,paste('b[',1:G,']', sep=''),'med']
	b_ints = parm_ests[i,paste('b[',1:G,']', sep=''),c('low95','up95')]
	b_range = range(c(b_ints, use_b))
	plot(use_b, b_ests, ylim=b_range)
	segments(use_b, b_ints[,'low95'], use_b, b_ints[,'up95'], lwd=2, col='grey50')
	points(use_b, b_ests, pch=16, cex=1)
	abline(0,1,col=2)
}
dev.off()

# Thin results
thinned = window(modbug, thin=1, start=100)

# Plot results of each chain
mylayout=c(2,5,9)
par(ask=T)
xyplot(thinned, layout=mylayout, scales=list(rot=0))
densityplot(thinned, layout=mylayout, scales=list(rot=0), aspect='fill')
acfplot(thinned, layout=mylayout, scales=list(rot=0), aspect='fill')
gelman.diag(thinned, autoburnin=F)

### Test performance of lmer
lmer_modlist = as.list(rep(NA, nrow(try_parms)))
reml_modlist = as.list(rep(NA, nrow(try_parms)))
lmer_ests = matrix(NA, nrow=nrow(try_parms), ncol=87)
reml_ests = matrix(NA, nrow=nrow(try_parms), ncol=87)
colnames(lmer_ests) = dimnames(parm_ests)[[2]]
colnames(reml_ests) = dimnames(parm_ests)[[2]]

# Using same simulated data as before
for(i in 1:nrow(try_parms)){
	y = Y[,i]

	this_mod = lmer(log(y) ~ 1 + xvar + (1|factor(genus)) + (1|factor(samp)), REML=F)
	lmer_modlist[[i]] = this_mod
	
	lmer_ests[i,paste('a[',1:J,']', sep='')] = as.matrix(ranef(this_mod, whichel='factor(samp)')[[1]])
	lmer_ests[i,paste('b[',1:G,']', sep='')] = as.matrix(ranef(this_mod, whichel='factor(genus)')[[1]])
	lmer_ests[i,c('g0','g1')] = fixef(this_mod)
	lmer_ests[i, 'deviance']=  deviance(this_mod)
	lmer_ests[i, c('sigma.a','sigma.b','sigma.y')] = as.data.frame(VarCorr(this_mod))[,'sdcor']
	
	this_mod = lmer(log(y) ~ 1 + xvar + (1|factor(genus)) + (1|factor(samp)), REML=T)
	reml_modlist[[i]] = this_mod
	
	reml_ests[i,paste('a[',1:J,']', sep='')] = as.matrix(ranef(this_mod, whichel='factor(samp)')[[1]])
	reml_ests[i,paste('b[',1:G,']', sep='')] = as.matrix(ranef(this_mod, whichel='factor(genus)')[[1]])
	reml_ests[i,c('g0','g1')] = fixef(this_mod)
	reml_ests[i, 'deviance']=  deviance(this_mod)
	reml_ests[i, c('sigma.a','sigma.b','sigma.y')] = as.data.frame(VarCorr(this_mod))[,'sdcor']
	
}

pdf('./Analysis/Figures/simulated data ml model performance coefs.pdf', height=15, width=15)
par(mfrow=c(3,3))
par(mar=c(2,2,.5,.5))
for(i in 1:nrow(try_parms)){
	plot(log(Y[,i])~xvar, ylab='', xlab='')

	parm_text = paste('g0=',try_parms[i,1],' g1=',try_parms[i,2],' sigma.y=',try_parms[i,3],' sigma.a=',try_parms[i,4],' sigma.b=',try_parms[i,5], sep='') 
	use_ests = format(lmer_ests[i,c('g0','g1','sigma.y','sigma.a','sigma.b')], digits=2)
	lmer_text = paste('g0=',use_ests[1],' g1=',use_ests[2],' sigma.y=',use_ests[3],' sigma.a=',use_ests[4],' sigma.b=',use_ests[5], sep='')
	use_ests = format(reml_ests[i,c('g0','g1','sigma.y','sigma.a','sigma.b')], digits=2)
	reml_text = paste('g0=',use_ests[1],' g1=',use_ests[2],' sigma.y=',use_ests[3],' sigma.a=',use_ests[4],' sigma.b=',use_ests[5], sep='')
	mtext(parm_text, 3, -1)
	mtext(est_text, 3, -2, col=2)	
	mtext(est_text, 3, -3, col='blue')

	abline(lm(log(Y[,i])~xvar), lwd=3, , col='grey30', lty=3) #Line from simple linear model	
	abline(try_parms[i,1], try_parms[i,2], lwd=3) # Line for simulated data
	abline(lmer_ests[i,'g0'], lmer_ests[i,'g1'], lwd=3, col=2)
	abline(reml_ests[i,'g0'], reml_ests[i,'g1'], lwd=3, lty=2, col='blue')
}
dev.off()

pdf('./Analysis/Figures/simulated data ml model performance genus coefs.pdf', height=9, width=9)
par(mfrow=c(3,3))
par(mar=c(2,2,.5,.5))
for(i in 1:nrow(try_parms)){
	use_b = B[,i]
	b_ests = lmer_ests[i,paste('b[',1:G,']', sep='')]
	bs = ranef(lmer_modlist[[i]], whichel='factor(genus)', condVar=T)[[1]]
	b_ints = b_ests+sqrt(as.numeric(attr(bs, 'postVar')))%*%t(c(-1.96, 1.96))
	b_range = range(c(b_ints, use_b))
	plot(use_b, b_ests, ylim=b_range)
	segments(use_b, b_ints[,1], use_b, b_ints[,2], lwd=2, col='grey50')
	points(use_b, b_ests, pch=16, cex=1)
	abline(0,1,col=2)

	b_ests = reml_ests[i,paste('b[',1:G,']', sep='')]
	bs = ranef(reml_modlist[[i]], whichel='factor(genus)', condVar=T)[[1]]
	b_ints = b_ests+sqrt(as.numeric(attr(bs, 'postVar')))%*%t(c(-1.96, 1.96))
	segments(use_b, b_ints[,1], use_b, b_ints[,2], lwd=2, col='lightblue', lty=2)
	points(use_b, b_ests, pch=1, cex=1, col='blue')
	abline(0,1,col=2)

}
dev.off()



################################################################
### Linear Models of Trait Variation at sample scale (trait diversity)
library(reshape) #melt

# Make a data frame of data for modeling
use_traits = c('Water_capacity','Thallus_thickness','STA','Cortex_thickness','Rhizine_length','Tot_chl_DW','Chla2b')
model_data = use_data[,c(use_traits, env_vars, 'Genus','SampID')]

# DONT DO THIS YET
# Transform responses
# Effectively we'll be fitting log-normal models since all responses are constrained to be positive
# Cortex thickness can be 0 so we add 1
#model_data$Cortex_thickness = model_data$Cortex_thickness + 1
#for(y in use_traits) model_data[,y] = log(model_data[,y])

# Re-scale predictors
model_data$Light_mean = model_data$Light_mean/10000
model_data$Temp_max = model_data$Temp_max/10


## Calculate metrics of single trait dispersion for each sample
samps = unique(model_data$SampID)
cft_disp = array(NA, dim=c(length(samps),7,2), dimnames=list(SampID=samps, Trait=use_traits, Metric = c('cv','mpd')))
for(i in use_traits){
	cft_disp[,i,'cv'] = tapply(model_data[,i], model_data$SampID, function(x) sqrt(var(x, na.rm=T))/mean(x, na.rm=T))
	cft_disp[,i,'mpd'] = tapply(model_data[,i], model_data$SampID, function(x) mean(dist(x[!is.na(x)])))
}

# Bootstrap null distribution for trait dispersion
N=10000
cft_disp_null = array(NA, dim=c(length(samps),7,2,N), dimnames=list(SampID=samps, Trait=use_traits, Metric= c('cv','mpd'), 1:N))

for(i in use_traits){
	NAinds = which(is.na(model_data[,i]))
	x = model_data[-NAinds, i]
	fact = model_data[-NAinds, 'SampID']
	
	for(j in 1:N){
		use_order = sample(fact)
		cv = tapply(x, use_order, function(y) sqrt(var(y, na.rm=T))/mean(y, na.rm=T))
		mpd = tapply(x, use_order, function(y) mean(dist(y[!is.na(y)])))
	
		cft_disp_null[names(cv),i,'cv',j] = cv
		cft_disp_null[names(mpd),i,'mpd',j] = mpd
	}
}	


# Calculate z-scores for actual trait dispersions based on null distribution
null_mean = apply(cft_disp_null, c(1,2,3), mean)
null_sd = apply(cft_disp_null, c(1,2,3), function(x) sqrt(var(x)))

cft_disp_z = (cft_disp - null_mean)/null_sd

plot(cft_disp_z[,'Water_capacity','mpd']~ env[paste('S',samps, sep=''),'Vpd_mean'])
plot(cft_disp_z[,'Cortex_thickness','mpd']~ env[paste('S',samps, sep=''),'Light_mean'])
plot(cft_disp_z[,'STA','mpd']~ env[paste('S',samps, sep=''),'Vpd_mean'])

cft_disp_df = melt(cft_disp)
cft_disp_df$z = melt(cft_disp_z)[,'value']

cft_disp_df = merge(cft_disp_df, unique(model_data[,c('SampID',env_vars)]))

pdf('./Analysis/Figures/FT mpd vs env.pdf', height=9, width=9)
par(mfrow=c(length(use_traits), length(env_vars)))
par(mar=c(4,4,1,1))
for(i in use_traits){
for(j in env_vars){

	this_data = subset(cft_disp_df, Trait==i&Metric=='mpd')
	plot(z~this_data[,j], data=this_data, xlab=j, ylab=i)
	abline(h=c(-2,2))

}}
dev.off()


#################################################################
### Assessing variation within species

library(DTK) # multiple comparison test for unequal variances

# Check trait variances within Genera
sapply(rownames(subset(traitnames, mode=='N')), function(x){
	tapply(use_data[,x], use_data$Genus, var, na.rm=T)
})

# Compare trait means across Genera (T3- Dunnett modified Tukey-Kramer test)
DTK_tests = sapply(rownames(subset(traitnames, mode=='N')), function(x){
	DTK = DTK.test(use_data[,x], use_data$Genus, a=0.05)
	DTK[[2]]
}, simplify=F)

# Pair-wise tests (DTK)
diff_pairs = sapply(rownames(subset(traitnames, mode=='N')), function(x){
	names(which(apply(DTK_tests[[x]][,2:3], 1, prod)>0))
})

# Number of observations of each trait for each genus
Nobs = sapply(rownames(traitnames), function(x){
	tapply(use_data[,x], use_data$Genus, function(y) sum(!is.na(y)))
})

# Wilcox rank-sum test for differences among pairs in ordinal data
diff_pairs_ord = sapply(rownames(subset(traitnames, mode=='O')), function(x){
	WT = pairwise.wilcox.test(as.numeric(use_data[,x]), use_data$Genus, p.adjust.method='hochberg', exact=F)$p.value
	apply(which(WT <0.05, arr.ind=T), 1, function(y) paste(rownames(WT)[y[1]],colnames(WT)[y[2]], sep='-'))
})

# Boxplots of trait variation across genera
traitnames = traitnames[order(traitnames$type),]

# Numeric traits
numtraits = rownames(subset(traitnames, mode=='N'))
for( i in numtraits) ){
	svg(paste('./Analysis/Figures/',i,' by Genus boxplot.svg', sep=''), height=5, width=4)
	par(mar=c(8,5,1,3))
	meds = tapply(use_data[,i], use_data$Genus, median, na.rm=T)
	plotorder = names(meds)[order(meds)]
	boxplot(use_data[,i]~factor(use_data$Genus, levels=plotorder), 
		las=3, ylab=traitnames[i,'displayName'], 
		lwd=1, pt.lwd=1, varwidth=T)
	axis(4)

	# ANOVA
	ftest = anova(lm(use_data[,i]~factor(use_data$Genus, levels=plotorder)))
	fval = round(ftest$'F value'[1],2)
	pval = round(ftest$'Pr(>F)'[1],3)
	usr = par('usr')
	text(usr[1], usr[4], paste('F =', fval, 'P =',pval), srt=90, adj=c(1.1,1.1))

	dev.off()
}

# Categorical traits
cattraits = rownames(subset(traitnames, mode=='O'))
for(i in cattraits){
	counts = sapply(tapply(use_data[,i], use_data$Genus, function(y) as.numeric(table(y))), 
		function(x) x)
	rownames(counts) = levels(use_data[,i])
	freqs = t(t(counts) / colSums(counts))
	means = tapply(as.numeric(use_data[,i]), use_data$Genus, mean, na.rm=T)
	plotorder = names(means)[order(means)]
	
	svg(paste('./Analysis/Figures/',i,' by Genus barplot.svg', sep=''), height=5, width=5)
	par(mar=c(4,8,3,4))
	use_col = colorRampPalette(mycol[2:9])(nrow(freqs))
	barplot(freqs[,plotorder], horiz=T, las=1, col=use_col, xlab='',
		width=log(colSums(counts)), main=traitnames[i, 'displayName'])
	par(xpd=NA)

	usr=par('usr')
	legend(0.5, usr[4], rownames(freqs), fill=use_col, horiz=T,
		xjust=.5, yjust=.5, bty='n')
	par(xpd=F)
	mtext('Proportion of Individuals', 1, 2.2)
	
	# Kruskal-Wallis One-way ANOVA
	krus = format(kruskal.test(use_data[,i], factor(use_data$Genus))$p.value, scientific=T, digits=T)

	mtexti(paste('Kruskal-Wallis: P =', krus), 4)

	dev.off()
}


###################################################################
### Distribution of species along environmental gradients

# Calculate site X genus matrix

counts  = data.frame(table(use_data$SampID, use_data$Genus))
names(counts) = c('SampID','Genus','Abun')

counts = merge(counts, env, all.x=T)

i='Cortex_thickness'
meds = tapply(use_data[,i], use_data$Genus, median, na.rm=T)
plotorder = names(meds)[order(meds)]
plotorder = plotorder[!(plotorder %in% c('Parmelinopsis','Myelochroa','Phaeophyscia'))]

pdf('./Analysis/Figures/genus frequency vs Vpd_mean.pdf',height=7, width=4)
layout(matrix(1:8, nrow=4, ncol=2, byrow=F))
par(mar=c(2,2.5,1,.5))
for(g in plotorder){
	this_data = subset(counts, Genus==g)
	plot(Abun/24~Vpd_mean, data=this_data, xlab='', ylab='', pch=16, ylim=c(0,0.5), las=1)
	mtext(g,3,0,adj=0, cex=.8)
}
dev.off()




svg('./Analysis/Figures/plot 3d light mean range temp.svg', height=4.8, width=6)
par(las=1)
scatterplot3d(samples$Light_range_sum,samples$Light_mean_sum,  samples$Temp_mean_sum, 
	type='h', pch=(15:17)[samples$TreeID], box=F, cex.symbols=1.5, 
	color=mycolor[c(1,8)][factor(samples$Pair)],
	xlab='Hours of full sunlight',ylab='Mean light intensity (lumens)', 
	zlab=expression(paste('Mean temp. ',degree,'C')), y.margin.add=.2
)
usr = par('usr')
uxlim = usr[2]-usr[1]
uylim = usr[4]-usr[3]
legB = legend(usr[1]+uxlim*.02,usr[4],rep('',3), col=mycolor[1], pch=15:17, pt.cex=1.5, bty='n')
legR = legend(usr[1]+uxlim*.06,usr[4],paste('Tree',1:3), col=mycolor[8], pch=15:17, pt.cex=1.5, bty='n')
par(xpd=T)
text(legB$rect$left+uxlim*0.05, legB$rect$top+uylim*0.03, 'Side', srt=60)
text(legR$rect$left+uxlim*0.06, legR$rect$top+uylim*0.02, 'Top', srt=60)
dev.off()


### OLD MODELS ###

# info for OpenBUGS
mod_data = list('y','samp','genus','n','J','G','temp','lmean','lrange')
mod_init = function(){list(a=rnorm(J), b=rnorm(G), g1=rnorm(1), g2=rnorm(1), g3=rnorm(1), 
	sigma.y=runif(1), sigma.a=runif(1), sigma.b=runif(1), mu.a=rnorm(1), mu.b=rnorm(1))}
mod_parms = c('a','b','sigma.y','sigma.a','sigma.b','mu.a','mu.b','g1','g2','g3')

# Run OpenBUGS
#modbug = bugs(mod_data, mod_init, mod_parms, "bayes_mod_genus_temp_light_mean_range.txt", 
#	n.chains=3, n.iter=100000, codaPkg=T, n.burnin=5000, n.thin=1)
#modbug = read.bugs(modbug)

#save(modbug, file=paste(i,'modbug1.RData', sep='_'))
modbug = load('Tot_chl_modbug1.RData')

# model 1b
mod_data = list('y','samp','genus','n','J','G','temp','lmean','lrange')
mod_init = function(){list(a=rnorm(J), b=rnorm(G), g1=rnorm(1), g2=rnorm(1), g3=rnorm(1), 
	sigma.y=runif(1), sigma.a=runif(1), sigma.b=runif(1), mu.a=rnorm(1))}
mod_parms = c('a','b','sigma.y','sigma.a','sigma.b','mu.a','g1','g2','g3')

# Run OpenBUGS
modbug1b = bugs(mod_data, mod_init, mod_parms, "bayes_mod_genus_temp_light_mean_range_1b.txt", 
	n.chains=3, n.iter=100000, codaPkg=T, n.burnin=10000, n.thin=1)
modbug1b = read.bugs(modbug1b)

save(modbug1b, file=paste(i,'modbug1b.RData', sep='_'))
modbug1b = load('Tot_chl_modbug1b.RData')


mod_data = list('y','samp','genus','n','J','G','temp','lmean','lrange')
mod_init = function(){list(a=rnorm(J), b=rnorm(G), g1=rnorm(1), g2=rnorm(1), g3=rnorm(1), 
	sigma.y=runif(1), sigma.a=runif(1), sigma.b=runif(1), mu.y=rnorm(1), mu.a=rnorm(1), mu.b=rnorm(1))}
mod_parms = c('a','b','sigma.y','sigma.a','sigma.b','mu.y','mu.a','mu.b','g1','g2','g3')
#modbug2 = bugs(mod_data, mod_init, mod_parms, "bayes_mod_genus_temp_light_mean_range_2.txt", 
#	n.chains=3, n.iter=100000, codaPkg=T, n.burnin=5000, n.thin=1)
#modbug2 = read.bugs(modbug2)

#save(modbug2, file=paste(i,'modbug2.RData', sep='_'))
modbug2 = load('Tot_chl_modbug2.RData')


# model with covariance between environmental predictors
prior_cor = 0.1
prior_var = 10000
prior_prec = matrix(prior_cor*prior_var, 3,3)
diag(prior_prec) <- prior_var
prior_phi = ginv(prior_prec)

mod_data = list('y','samp','genus','n','J','G','temp','lmean','lrange')
mod_init = function(){list(a=rnorm(J), b=rnorm(G), g=rnorm(3), 
	sigma.y=runif(1), sigma.a=runif(1), sigma.b=runif(1), 
	mu.y=rnorm(1), mu.g=rnorm(3), Phi=rWishart(1, 5, prior_phi),
	mu.g0=rep(0,3), prec=prior_prec, Phi0=prior_phi
)}
mod_parms = c('a','b','sigma.y','sigma.a','sigma.b','mu.y','g','mu.g','Sigma')

modbug3 = bugs(mod_data, mod_init, mod_parms, "bayes_mod_genus_temp_light_mean_range_3.txt", 
	n.chains=3, n.iter=1000, debug=T, n.burnin=0, n.thin=1)
modbug3 = read.bugs(modbug3)


# Questions:
# Clearly not converging b/c a and b dependent on one another
# However, parameters of interes, g, sigma are converging.
# Can I use these model results?


## Focal models with Mean temp and light range
focal_traits = c('Tot_chl','Chla2b','Water_capacity','Cortex_width')

i=focal_traits[1]
mod= lmer(use_data[,i]~ 1+Temp_mean_sum+Light_range_sum+Light_mean_sum+(1|Genus)+(1|SampID), data=use_data, REML=F)
anova(mod)








## Plots:
qqnorm(resid(mod_temp))
plot(fitted(mod_temp), resid(mod_temp))
abline(h=0)
plot(na.omit(use_data[,i]),fitted(mod_temp), xlim=c(50,350), ylim=c(50, 350))
abline(0,1)


plot(use_data[,i]~use_data$Temp_mean_sum)

######################################################
### Old Code

# Use orthogonal predictors from PCA of mean temp, light and light range
orth = prcomp(use_data[,c('Temp_mean_sum','Light_mean_sum','Light_range_sum')], center=T, scale=T)
pcscores = predict(orth)

pc1 = pcscores[keeprows,1]
pc2 = pcscores[keeprows,2]
pc3 = pcscores[keeprows,3]

mod_data = list('y','samp','genus','n','J','G','pc1','pc2','pc3')
mod_init = function(){list(a=rnorm(J), b=rnorm(G), g0=rnorm(1), g1=rnorm(1), 
	g2=rnorm(1), g3=rnorm(1), sigma.y=runif(1), sigma.a=runif(1), sigma.b=runif(1))}
mod_parms = c('a','b','g0','g1','g2','g3','sigma.y','sigma.a','sigma.b')

modbug_full = bugs(mod_data, mod_init, mod_parms, "bayes_mod_genus_full.txt", 
	n.chains=3, n.iter=200000, codaPkg=T, n.burnin=0,
	saveExec=T, working.directory='./WD5', restart=F)
modbug_full = read.bugs(modbug_full)

save(modbug_full, file=paste(i,'modbug_full.RData', sep='_'))

setwd('./Analysis/BUGS/WDtest/')
modbug_full = read.openbugs('')

### Plot model results ###
m='lmean'
use_mod = modbug_lmean
thinned = window(use_mod, thin=100, start=10000)
mylayout = c(2,5,9) # depends on number of parameters
par(ask=T)

pdf(paste(i,'modbug',m,'trace.pdf', sep='_'),height=11, width=8.5)
xyplot(thinned, layout=mylayout, scales=list(rot=0))
dev.off()

pdf(paste(i,'modbug',m,'density.pdf', sep='_'),height=11, width=8.5)
densityplot(thinned, layout=mylayout, scales=list(rot=0), aspect='fill')
dev.off()

pdf(paste(i,'modbug',m,'ac.pdf', sep='_'),height=11, width=8.5)
acfplot(thinned, layout=mylayout, scales=list(rot=0), aspect='fill')
dev.off()

gelman = gelman.diag(thinned) # near 1 implies convergence
gelman.diag(use_mod)

pdf(paste(i,'modbug',m,'gelman.pdf', sep='_'),height=11, width=8.5)
par(mfrow=c(3,2))
gelman.plot(use_mod)
dev.off()

## Plot slopes
Xlims = data.frame(Tot_chl_DW=c(-.5,.5), Cortex_width=c(-.3,.3), Water_capacity=c(-.3,.3))
Ylims = data.frame(temp=c(0,21), lmean=c(0,14))
displayNames=c(temp='Mean temp.',lmean='Mean light', Tot_chl_DW='Chlorophyll concentration',
	Cortex_width='Cortex thickness', Water_capacity='Water holding capacity')

svg('./Analysis/Figures/univariate model slope parameter posteriors.svg', height=3, width=8)
layout(matrix(1:6, ncol=3, nrow=2))
for(i in names(pooled_list)){
for(m in names(pooled_list[[1]])){
	
	if(m=='temp') par(mar=c(0,2.5,2.5,0))
	if(m=='lmean') par(mar=c(2.5,2.5,0,0))
	
	post = pooled_list[[i]][[m]][,'g1']	
	d = density(post)
	g1int = HPDinterval(post, 0.95)
	g1int50 = HPDinterval(post, 0.50)	
	plot(d, axes=F, lwd=2, xlab='', ylab='', main='', xlim=Xlims[,i],
		ylim=Ylims[,m])
	
	lines(c(g1int[1],g1int[2]),c(0,0), lwd=3, col='black')
	lines(c(g1int50[1],g1int50[2]),c(0,0), lwd=4, col='grey50')
	points(median(post),0, cex=2, pch='|', lwd=2)

	#par(xpd=T)
	abline(v=0, lty=1, col=mycolor[7], lwd=2)
	#par(xpd=F)
	
	if(m=='lmean') axis(1, line=.5)
	if(i=='Tot_chl_DW') mtext(displayNames[m], 2, 1)
	if(m=='temp') mtext(displayNames[i], 3, 1)
	
}}
dev.off()





