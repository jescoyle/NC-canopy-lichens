## This script models canopy lichen functional traits
options(stringsAsFactors=F)

working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/'
sql_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Data/SQLite Tables/'
data_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Data/Derived Tables/'
script_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Analysis/GitHub/NC-canopy-lichens/Code/'

setwd(working_dir)

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

xvarnames = c('Mean Light Intensity (Lux)', 'Freq. Full Sun', 
	'Max. Temperature (C*degree)','Mean VPD','Freq. Daytime Saturated Air')
names(xvarnames) = env_vars
xvarnames_short = c('Mean Light','Sun Freq.','Max Temp.','Mean VPD','Air Sat. Freq.')
names(xvarnames_short) = env_vars

######################################################
### Read in Data

traits = read.csv(paste(sql_dir, 'derived_traits.csv', sep=''))
env = read.csv(paste(data_dir, 'loggerdata.csv', sep=''))
samples = read.csv(paste(sql_dir, 'samples.csv', sep=''))
trees = read.csv(paste(sql_dir, 'trees.csv', sep=''))

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

# Remove thalli from samples on T11-T12
use_samps = subset(samples, TreeID <=10)$SampID
use_data = subset(use_data, SampID %in% use_samps)

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
use_data[use_data$Genus=='Usnea','Rhizine_length'] = NA # Usnea don't have rhizines

# Remove fruticose thalli (for some analyses)
#use_data = subset(use_data, Genus!='Usnea')

##########################################################################
### Make a data frame of data for modeling

# Create area normalized Water Holding Capacity and STM variables 
# WHC = Mass water / Area
# STM = Mass thallus / Area
use_data$WHC = use_data$Water_capacity/use_data$STA
use_data$STM = 1/use_data$STA

# Plot WHC vs STM
use_genera = names(table(use_data$Genus))[table(use_data$Genus)>5]
use_col = colorRampPalette(mycolor)(length(use_genera))

pdf('./Analysis/Figures/WHC vs STM by genus.pdf', height=5, width=5)
par(mar=c(4,4,1,1))
plot(WHC~STM, type='n', data=use_data, las=1, ylab=expression(WHC[A]), xlim=c(0,75), ylim=c(0,75))
abline(0,1)
for(i in 1:length(use_genera)){
	g = use_genera[i]
	these_data = subset(use_data, Genus==g)
	points(WHC~STM, data=these_data, bg=use_col[i], pch=21)
	use_mod = lm(WHC~STM, data=these_data)
	mod_func = function(x) predict(use_mod, data.frame(STM=x))
	curve(mod_func, from=min(these_data$STM, na.rm=T), to=max(these_data$STM, na.rm=T), add=T, col=use_col[i], lwd=2)
}
legend('bottomright', use_genera, lwd=2, col=use_col, bty='n')
dev.off()

colorby = cut(use_data$Area, breaks=100, include.lowest=T)
use_col = rainbow(100)
plot(WHC~STM, data=use_data, las=1, ylab=expression(WHC[A]), xlim=c(0,75), ylim=c(0,75), pch=16, col=use_col[colorby])

plot(Water_capacity ~ WHC, data=use_data)


use_traits = c('Water_capacity','WHC','Thallus_thickness','STA','STM','Cortex_thickness','Rhizine_length','Tot_chl_DW','Chla2b')
model_data = use_data[,c(use_traits, env_vars, 'Genus','SampID','Year')]

## Decision to Transform responses
# Effectively we'll be fitting log-normal models since all responses are constrained to be positive
# Cortex thickness can be 0 so we add 1
model_data$Cortex_thickness = model_data$Cortex_thickness + 1

# Distributions
layout(matrix(1:(2*length(use_traits)), nrow=2, byrow=F))
for(y in use_traits){
	hist(model_data[,y], main=y)
	hist(log(model_data[,y]), main=paste('log',y))
}

# Plot var ~ mean
layout(matrix(1:length(use_traits), nrow=1, byrow=F))
par(mar=c(3,3,3,1))
for(y in use_traits){
	#means = tapply(use_data[,y], use_data$SampID, mean, na.rm=T)
	#vars = tapply(use_data[,y], use_data$SampID, var, na.rm=T)
	#plot(vars~means, xlab='mean', ylab='var', main=y)
	means = tapply(log(use_data[,y]), use_data$SampID, mean, na.rm=T)
	vars = tapply(use_data[,y], use_data$SampID, var, na.rm=T)
	plot(vars~means, xlab='mean', ylab='var', main=paste('log',y))
}

# Transform data
for(y in use_traits) model_data[,y] = log(model_data[,y])

# Re-scale predictors
model_data$Light_high = model_data$Light_high*100
model_data$Vpd_daysatfreq = model_data$Vpd_daysatfreq*100
model_data$Light_mean = model_data$Light_mean/1000

##########################################################################
### Quick examination of the distribution of Genera/Species across samples

## See examine_community_structure.R for more involved analyses of community structure

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
taxon_tab = merge(taxon_tab, unique(use_data[,c('Genus','TaxonID')]))

#########################################################
## Correlations among traits
library(corrplot)

trait_cor = cor(model_data[,use_traits], use='pairwise.complete.obs')

pdf('./Analysis/Figures/trait correlations.pdf', height=7, width=10)
corrplot(trait_cor, method='ellipse', type='lower', order='hclust', p.mat=trait_cor, insig='p-value', sig.level=-1)
dev.off()

## Ordination of individuals
library(vegan)

noNA_data = na.omit(model_data)
traitord = rda(noNA_data[,use_traits], center=T, scale=T)

mycol=rainbow(11)

plot(traitord, type = "n")
points(traitord, display = "sites", cex = 0.8, pch=1, col=mycol[factor(noNA_data$Genus)])
text(traitord, display = "spec", cex=0.7, col="blue")

# Log-transformed STA and STM should just be the same but negative
# So, it is unnecessary to model both.
with(model_data, plot(STA~STM))

######################################################
## Spatial Variation in Environment

env$BranchPos = factor(env$BranchPos, levels=c('low','mid','high'))

# Calculate distance from top of tree
env$Height_top = sapply(1:nrow(env), function(i){
	trees[trees$TreeID==env[i,'TreeID'],'Height']-env[i,'Height']
}) 


## Models of environmental variation by height, branch position and site

env_mods = lapply(env_vars, function(x){
	lm(env[,x] ~ factor(Year)*Angle + factor(Year)*I(-1*Height_top) , data=env)
})
names(env_mods) = env_vars

modline = function(x, i, year, angle){
	this_mod = coef(env_mods[[i]])
	y = this_mod[1] + this_mod[2]*ifelse(year==2013, 0,1) + this_mod[3]*angle + 
		this_mod[4]*x + this_mod[5]*ifelse(year==2013, 0,1)*angle + 
		this_mod[6]*ifelse(year==2013, 0,1)*x
	y
}

# Plot effects
mycol=c('white','grey80')


svg('./Analysis/Figures/environmental variation by height.svg', height=5, width=12)
layout(matrix(1:(2*length(env_vars)), nrow=2))
for(i in env_vars){
	par(mar=c(3,3,4,1))
	plot(env[,i], -env$Height_top, xlab='', ylab='', main=xvarnames[i], las=1, pch=21, bg=mycol[factor(env$Year)])
	
	xvals = seq(min(-subset(env, Year==2013)$Height_top), max(-subset(env, Year==2013)$Height_top), length.out=100)
	yvals = modline(xvals, i, 2013, 0)	
	lines(yvals, xvals, lwd=4, col=1)
	lines(yvals, xvals, lwd=3, col=mycol[1])
	xvals = seq(min(-subset(env, Year==2014)$Height_top), max(-subset(env, Year==2014)$Height_top), length.out=100)
	yvals = modline(xvals, i, 2014, 0)	
	lines(yvals, xvals, lwd=4, col=1)
	lines(yvals, xvals, lwd=3, col=mycol[2])
	
	par(mar=c(5,3,2,1))
	plot(env[,i], -env$Height_top, xlab='', ylab='', las=1, pch=21, bg=mycol[factor(env$Year)])
	
	xvals = seq(min(-subset(env, Year==2013)$Height_top), max(-subset(env, Year==2013)$Height_top), length.out=100)
	yvals = modline(xvals, i, 2013, 90)	
	lines(yvals, xvals, lwd=4, col=1)
	lines(yvals, xvals, lwd=3, col=mycol[1])
	xvals = seq(min(-subset(env, Year==2014)$Height_top), max(-subset(env, Year==2014)$Height_top), length.out=100)
	yvals = modline(xvals, i, 2014, 90)	
	lines(yvals, xvals, lwd=4, col=1)
	lines(yvals, xvals, lwd=3, col=mycol[2])
}	
dev.off()

# Extract effect of height controling for other factors

height_effects = t(sapply(env_vars, function(i){
	this_mod = env_mods[[i]]
	new_mod = update(this_mod, .~.-I(-1*Height_top)-factor(Year):I(-1*Height_top))
	est = coef(summary(this_mod))[4,1:2]
	ftest = anova(new_mod, this_mod)[2,5:6]
	as.numeric(c(est,ftest))
}))
colnames(height_effects) = c('Est.','SE','F','P')
height_effects = data.frame(height_effects)

angle_effects = t(sapply(env_vars, function(i){
	this_mod = env_mods[[i]]
	new_mod = update(this_mod, .~.-Angle-factor(Year):Angle)
	est =	coef(summary(this_mod))[3,1:2]
	ftest = anova(new_mod, this_mod)[2,5:6]
	as.numeric(c(est,ftest))
}))
colnames(angle_effects) = c('Est.','SE','F','P')
angle_effects = data.frame(angle_effects)

height_effects$Predictor = 'Canopy Depth'
height_effects$Response = xvarnames[rownames(height_effects)]
angle_effects$Predictor = 'Surface Angle'
angle_effects$Response = xvarnames[rownames(angle_effects)]

env_tab = rbind(height_effects, angle_effects)
env_tab = env_tab[order(env_tab$Response, env_tab$Predictor), c(6,5,1:4)]

write.csv(env_tab, './Analysis/Figures/effect of height and angle on env.csv', row.names=F)


# Plot correlations among env variables
plot(env$Light_mean~env$Light_high)
plot(env$Light_mean~env$Light_p90)
plot(env$Vpd_mean~env$Vpd_satfreq)
plot(env$Vpd_mean~env$Vpd_daysatfreq)

# PCA of chosen variables
use_vars = c('Light_mean','Light_high','Temp_max','Vpd_mean','Vpd_daysatfreq')

envpca = prcomp(env[,use_vars], scale=T)
summary(envpca);envpca

cor(env[,env_vars])

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
### ANOVA of Individual Traits

library(lme4)
library(LMERConvenienceFunctions)
library(DTK) # multiple comparison test for unequal variances

## Test for trait differences across sites - definitely differences, but these may be due to environment
trait_aov_site = sapply(use_traits, function(i){
	mod = summary(aov(use_data[,i]~Year, data=use_data))[[1]]
	data.frame(SS=mod[1,2], MS=mod[1,3], F=mod[1,4], P=mod[1,5], N=sum(mod[,1])+1)
})
trait_aov_site = t(trait_aov_site)
trait_kw_site = sapply(use_traits, function(i){
	kwt = kruskal.test(use_data[,i], factor(use_data$Year))
	data.frame(Chi.sq=kwt$statistic, P=kwt$p.value)
})
trait_kw_site = t(trait_kw_site)

### ANOVA of traits among genera - using raw trait values, not transformed traits used for modeling

# Check trait variances within Genera - heterogeneous
sapply(use_traits, function(x){
	tapply(use_data[,x], use_data$Genus, var, na.rm=T)
})
# Check trait variances within Sites - heterogeneous
sapply(use_traits, function(x){
	tapply(use_data[,x], use_data$Year, var, na.rm=T)
})


# ANOVA among genera
trait_aov_gen = sapply(use_traits, function(i){
	mod = summary(aov(use_data[,i]~Genus, data=use_data))[[1]]
	c(TSS=sum(mod[,2]), SS=mod[1,2], MS=mod[1,3], F=mod[1,4], P=mod[1,5], N=sum(mod[,1])+1)
})
trait_aov_gen = t(trait_aov_gen)


# Kruskal-Wallis H-test in case variances are not equal
trait_kw = sapply(use_traits, function(i){
	kwt = kruskal.test(use_data[,i], factor(use_data$Genus))
	data.frame(Chi.sq=kwt$statistic, P=kwt$p.value)
})
trait_kw = t(trait_kw)

# using log-transformed traits
trait_aov_gen = sapply(use_traits, function(i){
	mod = summary(aov(model_data[,i]~Genus, data=model_data))[[1]]
	c(TSS=sum(mod[,2]), SS=mod[1,2], MS=mod[1,3], F=mod[1,4], P=mod[1,5], N=sum(mod[,1])+1)
})
trait_aov_gen = t(trait_aov_gen)

trait_kw = sapply(use_traits, function(i){
	kwt = kruskal.test(model_data[,i], factor(model_data$Genus))
	data.frame(Chi.sq=kwt$statistic, P=kwt$p.value)
})
trait_kw = t(trait_kw)

cbind(trait_aov_gen, trait_kw)

write.csv(cbind(trait_aov_gen, trait_kw), './Analysis/Figures/ANOVA and K-W test of traits by genus.csv')

# Plot variance explained
aov_props_gen = trait_aov_gen[,'SS']/trait_aov_gen[,'TSS']

barplot(aov_props_gen, ylim=c(0,1), col='black', las=2, ylab='Proportion SS Explained')



# Compare trait means across Genera (T3- Dunnett modified Tukey-Kramer multiple comparisons test)
DTK_tests = sapply(use_traits, function(x){
	this_data = subset(use_data, !is.na(use_data[,x])) # remove missing values
	keep_genera = names(which(table(this_data$Genus)>2)) # remove genera where variance can't be estimated
	this_data = subset(this_data, Genus %in% keep_genera)
	DTK = DTK.test(this_data[,x], this_data$Genus, a=0.05)
	DTK[[2]]
}, simplify=F)

# Pair-wise tests (DTK)
diff_pairs = sapply(use_traits, function(x){
	names(which(apply(DTK_tests[[x]][,2:3], 1, prod)>0))
})

# Number of observations of each trait for each genus
Nobs = sapply(rownames(traitnames), function(x){
	tapply(use_data[,x], use_data$Genus, function(y) sum(!is.na(y)))
})

write.csv(Nobs, './Analysis/Figures/number of trait measures for each genus.csv')

# Number of observation of each trait for each species
Nobs_sp = sapply(use_traits[1:7], function(x){
	tapply(use_data[,x], use_data$TaxonID, function(y) sum(!is.na(y)))
})

rownames(Nobs_sp) = paste(taxa[rownames(Nobs_sp),'Genus'], taxa[rownames(Nobs_sp),'Species'])

write.csv(Nobs_sp, './Analysis/Figures/number of trait measured for each species.csv')

# NOT DONE YET:  Wilcox rank-sum test for differences among pairs in ordinal data
diff_pairs_ord = sapply(rownames(subset(traitnames, mode=='O')), function(x){
	WT = pairwise.wilcox.test(as.numeric(use_data[,x]), use_data$Genus, p.adjust.method='hochberg', exact=F)$p.value
	apply(which(WT <0.05, arr.ind=T), 1, function(y) paste(rownames(WT)[y[1]],colnames(WT)[y[2]], sep='-'))
})


# Boxplots of trait variation across genera
traitnames = traitnames[order(traitnames$type),]

# Numeric traits
pdf('./Analysis/Figures/trait distributions by genus boxplot.pdf', height=5, width=4)
for( i in use_traits) {
	
	par(mar=c(8,5,1,3))
	meds = tapply(use_data[,i], use_data$Genus, median, na.rm=T)
	plotorder = names(meds)[order(meds)]
	
	plotlabel1 = traitnames[i,'exprName']
	plotlabel2 = paste('(',traitnames[i,'units'],')', sep='')
	ylab = ifelse(plotlabel2=='()', parse(text=plotlabel1), parse(text=paste(plotlabel1, plotlabel2, sep='~~')))
	
	boxplot(use_data[,i]~factor(use_data$Genus, levels=plotorder), 
		las=3, ylab=ylab, 
		lwd=1, pt.lwd=1, varwidth=T)
	axis(4)
	
	# Kruskal-Wallis Test
	kwtest = kruskal.test(use_data[,i], factor(use_data$Genus))
	chi2 = kwtest$statistic
	pval = kwtest$p.value

	usr = par('usr')
	this_text = substitute(Chi^2==c~~P==p, list(c = format(chi2, digits=3), p = format(pval, scientific=T, digits=2)))
	text(usr[1], usr[4], this_text, srt=90, adj=c(1.1,1.1))

	# Find groups that are different
	#this_diff = diff_pairs[[i]]
	#pairs = sapply(strsplit(this_diff, '-'), function(x) c(which(plotorder==x[1]), which(plotorder==x[2])))
	#pairs = apply(pairs, 2, function(x) x[order(x)])
	#pairs = pairs[,order(pairs[1,], pairs[2,])]
}
dev.off()

## Test trait differences just within Parmotrema
parm_data = subset(use_data, Genus=='Parmotrema')

trait_kw_parm = sapply(use_traits, function(i){
	kwt = kruskal.test(parm_data[,i], factor(parm_data$TaxonID))
	data.frame(Chi.sq=kwt$statistic, P=kwt$p.value)
})
trait_kw_parm = t(trait_kw_parm)

# Compare trait means across Parmotrema species
DTK_tests_parm = sapply(use_traits, function(x){
	this_data = subset(parm_data, !is.na(parm_data[,x])) # remove missing values
	keep_sp = names(which(table(this_data$TaxonID)>2)) # remove genera where variance can't be estimated
	this_data = subset(this_data, TaxonID %in% keep_sp)
	DTK = DTK.test(this_data[,x], this_data$TaxonID, a=0.05)
	DTK[[2]]
}, simplify=F)

# Pair-wise tests (DTK)
diff_pairs_parm = sapply(use_traits, function(x){
	names(which(apply(DTK_tests_parm[[x]][,2:3], 1, prod)>0))
})

pdf('./Analysis/Figures/trait distributions by Parmotrema species boxplot.pdf', height=5.5, width=4)
for( i in use_traits) {
	
	par(mar=c(5,5,1,3))
	meds = tapply(parm_data[,i], parm_data$TaxonID, median, na.rm=T)
	plotorder = names(meds)[order(meds)]
	
	plotlabel1 = traitnames[i,'exprName']
	plotlabel2 = paste('(',traitnames[i,'units'],')', sep='')
	ylab = ifelse(plotlabel2=='()', parse(text=plotlabel1), parse(text=paste(plotlabel1, plotlabel2, sep='~~')))
	
	boxplot(parm_data[,i]~factor(parm_data$TaxonID, levels=plotorder), 
		las=3, ylab=ylab, 
		lwd=1, pt.lwd=1, varwidth=T)
	axis(4)
	
	# Kruskal-Wallis Test
	kwtest = kruskal.test(parm_data[,i], factor(parm_data$TaxonID))
	chi2 = kwtest$statistic
	pval = kwtest$p.value

	usr = par('usr')
	this_text = substitute(Chi^2==c~~P==p, list(c = format(chi2, digits=3), p = format(pval, scientific=T, digits=2)))
	text(usr[1], usr[4], this_text, srt=90, adj=c(1.1,1.1))

}
dev.off()

## Boxplots of traits for all species within each genus
use_g = unique(use_data$Genus)

pdf('./Analysis/Figures/trait distributions by species.pdf', height=7, width=14)
for(g in use_g){
	this_data = subset(use_data, Genus==g)

	par(mar=c(7,3,3,1))
	par(mfrow=c(2,4))
	for(i in use_traits){
	if(i=='Rhizine_length'&g=='Usnea'){
		plot.new()
	} else {
		bp = boxplot(this_data[,i]~this_data$TaxonID, las=2, main=traitnames[i,'displayName'])
		text(1:length(bp$n), bp$stats[4,], bp$n, adj=c(-.25,-.25))
	
	}}
	plot.new()
}
dev.off()




# Number of observations of each trait for each genus
Nobs_parm = sapply(rownames(traitnames), function(x){
	tapply(parm_data[,x], parm_data$TaxonID, function(y) sum(!is.na(y)))
})

# Save significant differences
sink('./Analysis/Figures/pairwise trait differences among taxa.txt')
writeLines('All Genera\n')
diff_pairs
writeLines('\nParmotrema Species\n')
diff_pairs_parm
sink()

######
# OLD PLOTS WHICH I MAY RE-USE
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




######################################################
### Linear Models of Individual Traits

library(lme4)
library(LMERConvenienceFunctions)
library(boot) # for boot.ci

use_traits = c('Water_capacity','Thallus_thickness','STA','Cortex_thickness','Rhizine_length','Tot_chl_DW','Chla2b','WHC','STM')

# Create objects for storing models
modlist = vector('list', length(use_traits)*length(env_vars))
dim(modlist) = c(length(use_traits),length(env_vars))
dimnames(modlist) = list(use_traits, env_vars)

parm_ests = array(NA, dim=c(length(use_traits), length(env_vars), 6, 3), 
	dimnames = list(use_traits, env_vars, c('sigma.env', 'sigma.samp','sigma.genus','sigma.res','b0','b1'), c('est','low95','up95')))

genus_ests = array(NA, dim=c(length(use_traits), length(env_vars), length(unique(model_data$Genus)), 3), 
	dimnames = list(use_traits, env_vars, levels(factor(model_data$Genus)), c('est','low95','up95')))

# A function that calculated the fixed effects variance of a lmm model
calc_fixedvar = function(x){
	var(as.numeric(as.vector(fixef(x)) %*% t(model.matrix(x))))
}
# A function that extracts all variance components from a lmm model
get_allvar = function(x){
	var_df = data.frame(VarCorr(x))
	varcomps = c(calc_fixedvar(x), var_df$vcov)
	names(varcomps) = c('Fixed',var_df$grp)
	varcomps
}

## Loop through all traits and all env vars
# Note that number of observations and levels(SampID) will differ across models
for(i in use_traits){
for(j in env_vars){
	
	mod = lmer(model_data[,i] ~ 1 + model_data[,j]  + (1|Genus) + (1|SampID), data=model_data, REML=T)
	
	modlist[i,j][[1]] = mod
	
	# Calculate variance of fixed effects
	fixed_sigma2 = calc_fixedvar(mod)
	
	ests = c(fixed_sigma2, data.frame(VarCorr(mod))[,'vcov'], fixef(mod)) # Recording sigma (std dev not variance)
	g_ests = ranef(mod, whichel='Genus')[[1]]
	
	## CI METHOD 1 ##
	# Bootstrap estimate of confidence intervals for variance
	#varboot= bootMer(mod, get_allvar, nsim=1000, use.u=T, type='parametric')
	#ints_boot = sapply(1:4, function(x) boot.ci(varboot, index=x, type='perc')$percent[4:5])
	
	# Profile likelihood confidence intervals of fixed effects
	#ints = confint(mod, parm=4:5, method='profile')

	## CI METHOD 2 ##
	# Bootstrap fixed effects variance estimate only
	varboot = bootMer(mod, calc_fixedvar, nsim=1000, use.u=T, type='parametric')
	ints_boot = boot.ci(varboot, index=1, type='perc')$percent[4:5]
	
	# Profile likelihood confidence intervals for random and fixed effects
	pf = profile(mod, 1:5)
	vpf = varianceProf(pf) # converts to variance scale (by squaring)
	ints = confint(vpf, parm=1:5, method='profile')

	# 95% prediction intervals for random genus effects
	# May want to use bootstrap estimates later
	g_sds = sqrt(as.numeric(attr(ranef(mod, whichel='Genus', condVar=T)[[1]], 'postVar')))
	g_ints = g_ests$'(Intercept)' + g_sds%*%t(c(-1.96,1.96))
	
	# Save estimates
	parm_ests[i,j,,c('low95','up95')] = rbind(t(ints_boot), ints)
	parm_ests[i,j,,'est'] = ests
	genus_ests[i,j,rownames(g_ests),'est'] = g_ests$'(Intercept)'
	genus_ests[i,j,rownames(g_ests),c('low95','up95')] = g_ints
}}

# Save models 
#save(modlist, parm_ests, genus_ests, file='./Analysis/REML single variable FT models no Usnea.RData')
#save(modlist, parm_ests, genus_ests, file='./Analysis/REML single variable FT models with WHC STM.RData')


# Load previously saved models
load('./Analysis/REML single variable FT models.RData')
load('./Analysis/REML single variable FT models no Usnea.RData')
load('./Analysis/REML single variable FT models with WHC STM.RData')

## Plot estimated fixed effects
library(reshape)
library(lattice)

xvar_levels = factor(env_vars)

pdf('./Analysis/Figures/REML single variable model effects site effects no Usnea.pdf', height=9, width=4)
layout(matrix(1:length(use_traits), ncol=1))
par(mar=c(2,13,2,.5))
for(i in use_traits){
	this_data = parm_ests[i,,'b1',]
	xrange = range(this_data)
	plot(as.numeric(xvar_levels)~this_data[,'est'], xlim=xrange, axes=F, 
		xlab='', ylab='', lwd=1)
	abline(v=0, col='grey50', lty=2)
	segments(this_data[, 'low95'], as.numeric(xvar_levels), 
		this_data[, 'up95'], as.numeric(xvar_levels), lwd=1)
	axis(1)
	axis(2, at=1:length(env_vars), labels=xvarnames[xvar_levels], las=1)
	mtext(traitnames[i,'displayName'], 3, 0)
	box()
}
dev.off()

# Compare Mass and Area specific WHC
pdf('./Analysis/Figures/REML single variable model effects compare water capacity.pdf', height=7, width=4)
layout(matrix(1:4, ncol=1))
par(mar=c(2,13,2,.5))
for(i in c('Water_capacity','WHC','STA','STM')){
	this_data = parm_ests[i,,'b1',]
	xrange = range(this_data)
	plot(as.numeric(xvar_levels)~this_data[,'est'], xlim=xrange, axes=F, 
		xlab='', ylab='', lwd=1)
	abline(v=0, col='grey50', lty=2)
	segments(this_data[, 'low95'], as.numeric(xvar_levels), 
		this_data[, 'up95'], as.numeric(xvar_levels), lwd=1)
	axis(1)
	axis(2, at=1:length(env_vars), labels=xvarnames[xvar_levels], las=1)
	mtext(traitnames[i,'displayName'], 3, 0)
	box()
}
dev.off()


## Plot variance components (one chart for each predictor variable)


## Make figure for ESA presentation with cortex thickness
pdf(paste('./Analysis/Figures/variance components cortex models.pdf'), height=4, width=6)
varcomp = parm_ests['Cortex_thickness',c('Vpd_mean','Light_high'),c('sigma.env','sigma.samp','sigma.genus','sigma.res'),]
plot_data = t(varcomp[,,'est'])
colnames(plot_data) = c('VPD','Sun Freq.')

dotchart(plot_data, las=2, labels=c('Environment','Sample','Genus','Residual'), xlim=c(0, max(varcomp)), pch=16)
xvals = as.numeric(matrix(1:12, nrow=6)[1:4,])
segments(t(varcomp[2:1,,'low95']),xvals,t(varcomp[2:1,,'up95']),xvals)
mtext(expression(sigma^2), 1, 3, cex=1.5)
dev.off()

parm_ests[,,'sigma.env','est'] = parm_ests[,,'sigma.env','est']^2

## Plot Variance components: same scale across traits
for(i in env_vars){

varcomp = parm_ests[dim(parm_ests)[1]:1, i,c('sigma.res','sigma.samp','sigma.env','sigma.genus'),]
ests = varcomp[,,'est']
varcomp[ests==0] = NA # removes data for variance components that could not be estimated
ests = t(varcomp[,,'est'])
low95 = t(varcomp[,,'low95'])
up95 = t(varcomp[,,'up95'])
ests = rbind(ests, rep(NA, ncol(ests)))
low95 = rbind(low95, rep(NA, ncol(low95)))
up95 = rbind(up95, rep(NA, ncol(up95)))
yvals = (1:length(ests))-0.5

pdf(paste('./Analysis/Figures/variance components', i,'models.pdf'), height=7, width=6)
layout(matrix(1:2, nrow=1), widths=c(0.6, 0.4))
par(oma=c(0,2,0,2))
par(mar=c(3,0,3,0))
plot.new()
plot.window(ylim=c(0,length(ests)), xlim=c(0,1))
axis(2, at=yvals, labels=rep(c('Residual','Sample','Environment','Genus',''), ncol(ests)), las=1,
	tick=F, line=0, col=0, pos=1)
trait_text = rep(traitnames[colnames(ests),'exprName'], each=nrow(ests))
trait_text[!(1:length(ests) %in% seq(nrow(ests),length(ests),nrow(ests)))] = NA
for(i in 1:length(trait_text)) text(0.5, i-1.5, parse(text=trait_text[i]), pos=2, font=2)
abline(h=seq(0,length(ests),nrow(ests))-.5)

plot.new()
plot.window(ylim=c(0,length(ests)), xlim=c(0,max(up95, na.rm=T)))
abline(h=yvals, col='grey50', lty=3)
abline(h=seq(0,length(ests),nrow(ests))-.5)
segments(low95,yvals,up95,yvals, lend=1, lwd=3)
points(ests, yvals, pch=21, bg='white')
axis(1, pos=-.5, las=1)
mtext(expression(sigma^2), 1, 2, cex=1.5)
dev.off()

}

# Plot Variance components figure used in manuscript (Fig 2): different scales for each trait
for(i in env_vars){

varcomp = parm_ests[dim(parm_ests)[1]:1, i,c('sigma.res','sigma.samp','sigma.env','sigma.genus'),]
ests = varcomp[,,'est']
varcomp[ests==0] = NA # removes data for variance components that could not be estimated
ests = t(varcomp[,,'est'])
low95 = t(varcomp[,,'low95'])
up95 = t(varcomp[,,'up95'])
ests = rbind(ests, rep(NA, ncol(ests)))
low95 = rbind(low95, rep(NA, ncol(low95)))
up95 = rbind(up95, rep(NA, ncol(up95)))
yvals = (1:nrow(ests))-0.5

pdf(paste('./Analysis/Figures/variance components', i,'models no Usnea v2.pdf'), height=7, width=6)
layout(matrix((2*ncol(ests)):1, nrow=ncol(ests), byrow=T)[,2:1], widths=c(0.6, 0.4))
par(oma=c(2,0,0,2))
par(mar=c(2,0,0,0))

for(j in colnames(ests)){
plot.new()
plot.window(ylim=c(-.5,nrow(ests)), xlim=c(0,1))
axis(2, at=yvals, labels=c('Residual','Sample','Environment','Genus',''), las=1,
	tick=F, line=0, col=0, pos=1.05)
trait_text = traitnames[j,'exprName']
text(0.7, nrow(ests)-1.5, parse(text=trait_text), pos=2, font=2)
segments(.3, c(0,nrow(ests))-.5, 1.1, c(0,nrow(ests))-.5)

plot.new()
plot.window(ylim=c(-.5,nrow(ests)), xlim=c(0,max(up95[,j], na.rm=T)))
abline(h=yvals, col='grey50', lty=3)
abline(h=c(0,nrow(ests))-.5)
segments(low95[,j], yvals, up95[,j], yvals, lend=1, lwd=3)
points(ests[,j], yvals, pch=21, bg='white',cex=1.5)
axis(1, pos=-.5, las=1)
if(j==colnames(ests)[1]) mtext(expression(paste('Estimated Variance ',(sigma^2))), 1, 2.5, cex=.8)
}

dev.off()

}

# Compare water capacity measured
wtraits = c('Water_capacity','WHC')#,'STA','STM')

varcomp = parm_ests[wtraits, 'Vpd_mean',c('sigma.res','sigma.samp','sigma.env','sigma.genus'),]
ests = varcomp[,,'est']
varcomp[ests==0] = NA # removes data for variance components that could not be estimated
ests = t(varcomp[,,'est'])
low95 = t(varcomp[,,'low95'])
up95 = t(varcomp[,,'up95'])
ests = rbind(ests, rep(NA, ncol(ests)))
low95 = rbind(low95, rep(NA, ncol(low95)))
up95 = rbind(up95, rep(NA, ncol(up95)))
yvals = (1:nrow(ests))-0.5

pdf(paste('./Analysis/Figures/variance components Vpd_mean models compare water capacity.pdf'), height=3, width=6)
layout(matrix((2*ncol(ests)):1, nrow=ncol(ests), byrow=T)[,2:1], widths=c(0.6, 0.4))
par(oma=c(2,0,0,2))
par(mar=c(2,0,0,0))

for(j in colnames(ests)){
plot.new()
plot.window(ylim=c(-.5,nrow(ests)), xlim=c(0,1))
axis(2, at=yvals, labels=c('Residual','Sample','Environment','Genus',''), las=1,
	tick=F, line=0, col=0, pos=1.05)
trait_text = traitnames[j,'exprName']
text(0.7, nrow(ests)-1.5, parse(text=trait_text), pos=2, font=2)
segments(.3, c(0,nrow(ests))-.5, 1.1, c(0,nrow(ests))-.5)

plot.new()
plot.window(ylim=c(-.5,nrow(ests)), xlim=c(0,max(up95[,j], na.rm=T)))
abline(h=yvals, col='grey50', lty=3)
abline(h=c(0,nrow(ests))-.5)
segments(low95[,j], yvals, up95[,j], yvals, lend=1, lwd=3)
points(ests[,j], yvals, pch=21, bg='white',cex=1.5)
axis(1, pos=-.5, las=1)
if(j==colnames(ests)[1]) mtext(expression(paste('Estimated Variance ',(sigma^2))), 1, 2.5, cex=.8)
}

dev.off()



## Plot estimated random effects for each genus
# From models of Light_mean
j = 'Light_mean'

xvar_levels = factor(dimnames(genus_ests)[[3]])

pdf('./Analysis/Figures/REML genus random effects from Light_mean models.pdf', height=12, width=4)
layout(matrix(1:length(use_traits), ncol=1))
par(mar=c(2,13,2,.5))
for(i in use_traits){
	this_data = genus_ests[i,j,,]
	plot_order = rank(this_data[,'est'])
	xrange = range(this_data, na.rm=T)
	plot(plot_order~this_data[,'est'], xlim=xrange, axes=F, 
		xlab='', ylab='', lwd=1)
	#abline(v=0, col='grey50', lty=2)
	segments(this_data[, 'low95'], plot_order, 
		this_data[, 'up95'], plot_order, lwd=1)
	axis(1)
	axis(2, at=as.numeric(xvar_levels), labels=xvar_levels[order(plot_order)], las=1)
	mtext(traitnames[i,'displayName'], 3, 0)
	box()
}
dev.off()

# CHLA2B MODELS CAN'T ESTIMATE GENUS EFFECTS 
this_mod = modlist['Chla2b','Vpd_mean'][[1]]
table(this_mod@frame$Genus) # only 1 Phaeophyscia and 2 Myelochroa

mod = lmer(Chla2b~Light_mean + (1|Genus) + (1|SampID), data=subset(model_data, !(Genus %in% c('Phaeophyscia','Myelochroa','Parmelinopsis','Physcia'))), REML=T)
mod1 = lm(Chla2b~Vpd_mean, data=model_data)

mod = lmer(Chla2b~Light_mean + (1|Genus) + (1|SampID), data=subset(model_data, Genus %in% c('Punctelia','Parmotrema')), REML=T)
mod = lmer(Chla2b~Light_mean + (1|SampID), data=model_data, REML=T)


### Run models that control for environmental differences across sites
modlist_site = vector('list', length(use_traits)*length(env_vars))
dim(modlist_site) = c(length(use_traits),length(env_vars))
dimnames(modlist_site) = list(use_traits, env_vars)

parm_ests_site = array(NA, dim=c(length(use_traits), length(env_vars), 6, 3), 
	dimnames = list(use_traits, env_vars, c('sigma.samp','sigma.genus','sigma.res','b0','site','b1'), c('est','low95','up95')))

genus_ests_site = array(NA, dim=c(length(use_traits), length(env_vars), length(unique(model_data$Genus)), 3), 
	dimnames = list(use_traits, env_vars, levels(factor(model_data$Genus)), c('est','low95','up95')))

for(i in use_traits){
for(j in env_vars){
	
	mod = lmer(model_data[,i] ~ 1 + factor(Year) + model_data[,j]  + (1|Genus) + (1|SampID), data=model_data, REML=T)
	
	modlist_site[i,j][[1]] = mod
	
	ests = c(data.frame(VarCorr(mod))[,'sdcor'], fixef(mod))
	g_ests = ranef(mod, whichel='Genus')[[1]]
	ints = confint(mod, parm=1:6, method='profile')
	g_sds = sqrt(as.numeric(attr(ranef(mod, whichel='Genus', condVar=T)[[1]], 'postVar')))
	g_ints = g_ests$'(Intercept)' + g_sds%*%t(c(-1.96,1.96))
	
	parm_ests_site[i,j,,'est'] = ests
	parm_ests_site[i,j,,c('low95','up95')] = ints
	genus_ests_site[i,j,rownames(g_ests),'est'] = g_ests$'(Intercept)'
	genus_ests_site[i,j,rownames(g_ests),c('low95','up95')] = g_ints
}}

# Save models
#save(modlist_site, parm_ests_site, genus_ests_site, file='./Analysis/REML single variable FT models with site effects.RData')
#save(modlist_site, parm_ests_site, genus_ests_site, file='./Analysis/REML single variable FT models with site effects with WHC STM.RData')


# Load previously saved models
#load('./Analysis/REML single variable FT models with site effects.RData')
load('./Analysis/REML single variable FT models with site effects with WHC STM.RData')

# Compare site effects across models
parm_ests_site[,,'site',]

## Make plot of estimated effects with and without controlling for site

xvar_levels = factor(env_vars)
ynudge = 0.4

# Figure grouping variables by trait
pdf('./Analysis/Figures/REML single variable model effects compare site.pdf', height=7, width=7)
layout(matrix(1:(length(use_traits)+1), ncol=2))
par(oma=c(2,0,0,0))
par(mar=c(2,7,2,.5))
par(lend=1)
count = 1
for(i in use_traits){
	
	this_data = parm_ests[i,,'b1',]
	this_site = parm_ests_site[i,,'b1',]
	xrange = range(c(this_data, this_site))

	plot.new()
	plot.window(xlim=xrange, ylim=c(0,length(env_vars)))
	abline(v=0, col='black', lty=1, lwd=.8)
	
	segments(this_data[, 'low95'], as.numeric(xvar_levels), 
		this_data[, 'up95'], as.numeric(xvar_levels), lwd=2)
	points(this_data[,'est'], as.numeric(xvar_levels), pch=21, bg='white')

	segments(this_site[, 'low95'], as.numeric(xvar_levels)-ynudge, 
		this_site[, 'up95'], as.numeric(xvar_levels)-ynudge, lwd=2, col='grey50')
	points(this_site[,'est'], as.numeric(xvar_levels)-ynudge, pch=21, col='grey50', bg='white')
	

	axis(1)
	axis(2, at=1:length(env_vars), labels=xvarnames_short[xvar_levels], las=1)
	mtext(parse(text=traitnames[i,'exprName']), 3, 0, cex=0.8)
	box()
	
	if(count %in% c(4,7)) mtext('Estimated Effect', 1, 2.5, cex=.8)
	count = count+1
}

dev.off()

# Figure grouping models by environmental predictor
N = length(use_traits)
trait_order = use_traits[N:1]

# Factor by which to multiply estimates so that modeled parameter estimates match units in manuscript
xfact = c(1,1,1,100,1)
names(xfact) = env_vars

# Units for env variables
envunits = c('klx','%','degree*C','Pa','%')
names(envunits) = env_vars

pdf('./Analysis/Figures/REML single variable model effects compare site by env.pdf', height=7, width=7)
layout(matrix(1:(length(env_vars)+1), ncol=2))
par(mar=c(4,8,1,.5))
par(lend=1)

for(j in env_vars){
	
	this_data = parm_ests[trait_order,j,'b1',]*xfact[j]
	this_site = parm_ests_site[trait_order,j,'b1',]*xfact[j]
	xrange = range(c(this_data, this_site))

	plot.new()
	plot.window(xlim=xrange, ylim=c(0,N))
	abline(v=0, col='black', lty=1, lwd=.8)
	
	segments(this_data[, 'low95'], 1:N, 
		this_data[, 'up95'], 1:N, lwd=2)
	points(this_data[,'est'], 1:N, pch=21, bg='white')

	segments(this_site[, 'low95'], (1:N)-ynudge, 
		this_site[, 'up95'], (1:N)-ynudge, lwd=2, col='grey50')
	points(this_site[,'est'], (1:N)-ynudge, pch=21, col='grey50', bg='white')
	
	axis(1)
	axis(2, at=1:N, labels=parse(text=traitnames[trait_order,'exprName']), las=1)
	mtext(parse(text=traitnames[i,'exprName']), 3, 0, cex=0.8)

	
	this_xname = paste('Effect of',xvarnames_short[j])
	#this_unit = parse(text=envunits[j])
	#this_lab = paste(this_xname,'~~(',this_unit,')', sep='')
	#this_lab = bquote(.(this_xname)~~(.(this_unit))) # Can't get this to work

	mtext(this_xname, 1, 2.5, cex=0.8)
	box()
}

dev.off()




## Multi-panel plot of model predictions with and without site effects
## Manuscript Figure 1.

# Make vector for converting predictors in model data back to original scale for plotting
xvar_factor = c(0.001, 100, 1, 1, 100); names(xvar_factor) = env_vars

# A function that fits the mean line for the models without site effects
modline = function(x, i, j){
	y = parm_ests[i,j,'b0','est'] + parm_ests[i,j,'b1','est']*(x*xvar_factor[j])
	y = exp(y)
	if(i=='Cortex_thickness') y = y-1
	y
}

# A function that fits separate lines for each site
modline_site = function(x, i, j, year){
	y = parm_ests_site[i,j,'b0','est'] +  ifelse(year==2013, 0, 1)*parm_ests_site[i,j,'site','est'] + parm_ests_site[i,j,'b1','est']*(x*xvar_factor[j])
	y = exp(y)
	if(i=='Cortex_thickness') y = y-1
	y
}

use_col = c('white','grey70')

svg('./Analysis/Figures/compare site effects vpd_mean.svg', height=5.5, width=4.5)
par(mfrow=c(3,2))
par(mar=c(2,4.4,1,1))
par(lend=2)
for(i in c('Water_capacity','STA','Rhizine_length','Cortex_thickness','Tot_chl_DW','Chla2b')){

	this_unit = parse(text=traitnames[i,'units'])
	this_traitname = parse(text=traitnames[i,'exprName'])
	this_lab = ifelse(length(this_unit)>0, parse(text=paste(this_traitname,'~~(',this_unit,')', sep='')), this_traitname)
	
	plot(use_data[,i]~use_data$Vpd_mean, pch=21, bg=use_col[factor(use_data$Year)],las=1, 
		ylab=this_lab, xlab='', axes=F)
	
	# Significance of fixed effects slopes
	this_mod = modlist[i,'Vpd_mean'][[1]]
	drop_mod = update(this_mod, .~.-model_data[, j])
	pval1 = anova(this_mod, drop_mod)[2,'Pr(>Chisq)']

	this_mod = 	modlist_site[i,'Vpd_mean'][[1]]
	drop_mod = update(this_mod, .~.-model_data[, j])
	pval2 = anova(this_mod, drop_mod)[2,'Pr(>Chisq)']

	use_lty = c(1,6)[1+(c(pval1, pval2)>=0.05)]

	curve(modline(x, i, 'Vpd_mean'), from=min(use_data$Vpd_mean), to = max(use_data$Vpd_mean), add=T, lwd=4, lty=use_lty[1])
	curve(modline_site(x, i, 'Vpd_mean', 2013), from=min(subset(use_data, Year==2013)$Vpd_mean), to = max(subset(use_data, Year==2013)$Vpd_mean),
		add=T, lwd=3, col='grey30')
	curve(modline_site(x, i, 'Vpd_mean', 2013), from=min(subset(use_data, Year==2013)$Vpd_mean), to = max(subset(use_data, Year==2013)$Vpd_mean),
		add=T, lwd=2, col=use_col[1], lty=use_lty[2])
	curve(modline_site(x, i, 'Vpd_mean', 2014), from=min(subset(use_data, Year==2014)$Vpd_mean), to = max(subset(use_data, Year==2014)$Vpd_mean),
		add=T, lwd=3, col='grey30')
	curve(modline_site(x, i, 'Vpd_mean', 2014), from=min(subset(use_data, Year==2014)$Vpd_mean), to = max(subset(use_data, Year==2014)$Vpd_mean),
		add=T, lwd=2, col=use_col[2], lty=use_lty[2])

	axis(1, at=5:8, labels=(5:8)*100)
	axis(2, las=1)
	box()

}
dev.off()

# Compare water capacity
svg('./Analysis/Figures/compare site effects vpd_mean compare water capacity.svg', height=4.5, width=5)
par(mfrow=c(2,2))
par(mar=c(2,4.4,1,1))
par(lend=2)
for(i in c('Water_capacity','WHC','STA','STM')){

	this_unit = parse(text=traitnames[i,'units'])
	this_traitname = parse(text=traitnames[i,'exprName'])
	this_lab = ifelse(length(this_unit)>0, parse(text=paste(this_traitname,'~~(',this_unit,')', sep='')), this_traitname)
	
	plot(use_data[,i]~use_data$Vpd_mean, pch=21, bg=use_col[factor(use_data$Year)],las=1, 
		ylab=this_lab, xlab='', axes=F)
	
	# Significance of fixed effects slopes
	this_mod = modlist[i,'Vpd_mean'][[1]]
	drop_mod = update(this_mod, .~.-model_data[, j])
	pval1 = anova(this_mod, drop_mod)[2,'Pr(>Chisq)']

	this_mod = 	modlist_site[i,'Vpd_mean'][[1]]
	drop_mod = update(this_mod, .~.-model_data[, j])
	pval2 = anova(this_mod, drop_mod)[2,'Pr(>Chisq)']

	use_lty = c(1,6)[1+(c(pval1, pval2)>=0.05)]

	curve(modline(x, i, 'Vpd_mean'), from=min(use_data$Vpd_mean), to = max(use_data$Vpd_mean), add=T, lwd=4, lty=use_lty[1])
	curve(modline_site(x, i, 'Vpd_mean', 2013), from=min(subset(use_data, Year==2013)$Vpd_mean), to = max(subset(use_data, Year==2013)$Vpd_mean),
		add=T, lwd=3, col='grey30')
	curve(modline_site(x, i, 'Vpd_mean', 2013), from=min(subset(use_data, Year==2013)$Vpd_mean), to = max(subset(use_data, Year==2013)$Vpd_mean),
		add=T, lwd=2, col=use_col[1], lty=use_lty[2])
	curve(modline_site(x, i, 'Vpd_mean', 2014), from=min(subset(use_data, Year==2014)$Vpd_mean), to = max(subset(use_data, Year==2014)$Vpd_mean),
		add=T, lwd=3, col='grey30')
	curve(modline_site(x, i, 'Vpd_mean', 2014), from=min(subset(use_data, Year==2014)$Vpd_mean), to = max(subset(use_data, Year==2014)$Vpd_mean),
		add=T, lwd=2, col=use_col[2], lty=use_lty[2])

	axis(1, at=5:8, labels=(5:8)*100)
	axis(2, las=1)
	box()

}
dev.off()


## Plot cortex thickness models with different lines for each genus
genera = levels(factor(use_data$Genus))
use_col = colorRampPalette(mycolor)(length(genera)); names(use_col) = genera

modline_g = function(x, i, j, g){
	y = parm_ests[i,j,'b0','est'] + parm_ests[i,j,'b1','est']*(x*xvar_factor[j]) + genus_ests[i,j,g,'est']
	y = exp(y)
	if(i=='Cortex_thickness') y = y-1
	y
}

svg('./Analysis/Figures/model predictions cortex thickness light_high vpd_mean.svg', height=4, width=8)
par(mfrow=c(1,2))
par(mar=c(4,3,1,1))
for(j in c('Vpd_mean','Light_high')){
	plot(use_data$Cortex_thickness ~ use_data[,j], pch=21, bg = use_col[use_data$Genus], las=1,
		ylab='', xlab=xvarnames[j])

	for(g in genera){
		curve(modline_g(x, 'Cortex_thickness', j, g), from=min(use_data[,j]), to = max(use_data[,j]), 
			add=T, col=use_col[g])
	}

	curve(modline(x, 'Cortex_thickness',j), from=min(use_data[,j]), to = max(use_data[,j]), 
		add=T, lwd=3)

}
dev.off()




### Examine random slope models: effects of env change with genus

# Plot trait ~ env relationships for different genera

# Loop through genera
genera = levels(factor(model_data$Genus))
use_col = colorRampPalette(mycolor)(length(genera)); names(use_col) = genera


pdf('./Analysis/Figures/trait-env relationships by genus.pdf', height=6, width=8)
par(mfrow=c(2,3))
#layout(matrix(1:3, nrow=1), widths=c(.4, .4, .2))

for(i in use_traits){
par(mar=c(4,4,1,1))
for(j in env_vars){
	
	plot(use_data[,c(j,i)], las=1, xlab=xvarnames[j], ylab=traitnames[i,'displayName'],
		pch=21, bg=use_col[factor(use_data$Genus)], cex=.8)
	
	for(g in genera){
		this_data = subset(model_data, Genus==g)
		xvar = this_data[,j]
		yvar = this_data[,i] 
		if(length(yvar[!is.na(yvar)]) > 3){
			this_mod = lm(yvar~xvar)
			this_func = function(x){
				new_x = x*xvar_factor[j]
				y = exp(predict(this_mod, data.frame(xvar=new_x)))
				if(i=='Cortex_thickness') y = y+1
				y
			}
	
			sig = anova(this_mod)[1,'Pr(>F)']<0.05			
			use_lwd = ifelse(sig, 2, 1)
			use_lty = ifelse(sig, 1, 2)

			curve(this_func, from=min(xvar)/xvar_factor[j], to=max(xvar)/xvar_factor[j], add=T,
				lwd=use_lwd, lty=use_lty, col=use_col[g])
		}
	}
}
plot.new()
par(mar=c(4,1,1,1))
legend('right', genera, col=use_col, lwd=2, bty='n')
}
dev.off()


# Only plot the two that differ: WHC ~ Mean light and CHLTOT ~ Mean light
pdf('./Analysis/Figures/trait-env relationship by genus only sig.pdf', height=3, width=6)
layout(t(matrix(1:3)), widths=c(0.4,0.4,0.2))
par(mar=c(4,4,1,1))
for(i in c('Water_capacity', 'Tot_chl_DW')){

	plotlabel1 = traitnames[i,'exprName']
	plotlabel2 = paste('(',traitnames[i,'units'],')', sep='')
	ylab = ifelse(plotlabel2=='()', parse(text=plotlabel1), parse(text=paste(plotlabel1, plotlabel2, sep='~~')))
	

	plot(use_data[,c('Light_mean',i)], las=1, xlab=xvarnames['Light_mean'], ylab=ylab,
		pch=21, bg=use_col[factor(use_data$Genus)], cex=.8)
	
	for(g in genera){
		this_data = subset(model_data, Genus==g)
		xvar = this_data[,'Light_mean']
		yvar = this_data[,i] 
		if(length(yvar[!is.na(yvar)]) > 3){
			this_mod = lm(yvar~xvar)
			this_func = function(x){
				new_x = x*xvar_factor['Light_mean']
				y = exp(predict(this_mod, data.frame(xvar=new_x)))
				if(i=='Cortex_thickness') y = y+1
				y
			}
	
			sig = anova(this_mod)[1,'Pr(>F)']<0.05			
			use_lwd = ifelse(sig, 2, 1)
			use_lty = ifelse(sig, 1, 2)

			curve(this_func, from=min(xvar)/xvar_factor['Light_mean'], to=max(xvar)/xvar_factor['Light_mean'], add=T,
				lwd=use_lwd, lty=use_lty, col=use_col[g])
		}
	}
}
plot.new()
par(mar=c(4,1,1,1))
legend('right', genera, col=use_col, lwd=2, bty='n')
dev.off()

## Make table of separate means model coefficients
sepgen_mods = array(NA, dim=c(length(use_traits), length(env_vars), length(genera)+1, 4), 
	dimnames=list(Trait=use_traits, Predictor=env_vars, Genus=c('All',genera), Statistic=c('b0','b1','F','P')))

for(i in use_traits){
for(j in env_vars){
for(g in genera){
	this_data = subset(model_data, Genus==g)
	xvar = scale(this_data[,j], center=T, scale=F) # so we can compare coefs
	yvar = this_data[,i] 
	
	if(length(yvar[!is.na(yvar)]) > 3){
		this_mod = lm(yvar~xvar)
		sepgen_mods[i,j,g,'F'] = anova(this_mod)[1,'F value']	
		sepgen_mods[i,j,g,'P'] = anova(this_mod)[1,'Pr(>F)']
		sepgen_mods[i,j,g,c('b0','b1')] = coef(this_mod)			
	}
}
	
	yvar = scale(model_data[,j], center=T, scale=F)
	xvar = model_data[,i]

	this_mod = lm(yvar~xvar)
	sepgen_mods[i,j,'All','F'] = anova(this_mod)[1,'F value']	
	sepgen_mods[i,j,'All','P'] = anova(this_mod)[1,'Pr(>F)']
	sepgen_mods[i,j,'All',c('b0','b1')] = coef(this_mod)

}}

# b0 gives genus means, b1 gives slopes
sepgen_melt = melt(sepgen_mods)
gen_pvals = cast(sepgen_melt, Trait+Predictor~Genus, subset=Statistic=='P')
gen_slopes = cast(sepgen_melt, Trait+Predictor~Genus, subset=Statistic=='b1')
write.csv(cbind(gen_slopes, gen_pvals), './Analysis/Figures/trait env models by genus slopes and pvals.csv', row.names=F)

# Write out parameter estimates from models
names(dimnames(parm_ests)) = c('Trait','Predictor','Statistic','Value')
parm_melt = melt(parm_ests)
parm_df = cast(parm_melt, Trait+Predictor~Statistic, subset=Value=='est')
write.csv(parm_df, './Analysis/Figures/individual trait models parameter estimates.csv' ,row.names=F)

# Based on these figure it is probably not necessary to fit a random slopes model.
# Also, there isn't enough data across levels to fit a random slopes model (we get cor=-1 between slope and intercept)
i = 'WHC'
j = 'Vpd_mean'

# Compare different model formulations
yvar = scale(model_data[,i], center=T, scale=F) # center response so that we can estimate random slopes
lme0 = lmer(yvar ~ model_data[,j] + (1|Genus), data=model_data, REML=T)
lme1 = lmer(yvar ~ model_data[,j] + (model_data[,j]|Genus), data=model_data, REML=T)
lme2 = lmer(yvar ~ model_data[,j] + (1|Genus) + (1|SampID), data=model_data, REML=T)
AIC(lme0, lme1, lme2)
VarCorr(lme1)



#### Run models for Duke Forest data only in case VPD relationships are driven by site-level differences
duke_data = subset(model_data, Year==2014)

# Create objects for storing models
modlist_duke = vector('list', length(use_traits)*length(env_vars))
dim(modlist_duke) = c(length(use_traits),length(env_vars))
dimnames(modlist_duke) = list(use_traits, env_vars)

parm_ests_duke = array(NA, dim=c(length(use_traits), length(env_vars), 5, 3), 
	dimnames = list(use_traits, env_vars, c('sigma.samp','sigma.genus','sigma.res','b0','b1'), c('est','low95','up95')))

genus_ests_duke = array(NA, dim=c(length(use_traits), length(env_vars), length(unique(model_data$Genus)), 3), 
	dimnames = list(use_traits, env_vars, levels(factor(model_data$Genus)), c('est','low95','up95')))

# Loop through all traits and all env vars
for(i in use_traits){
for(j in env_vars){
	mod = lmer(duke_data[,i] ~ 1 + duke_data[,j] + (1|Genus) + (1|SampID), data=duke_data, REML=T)
	modlist_duke[i,j][[1]] = mod
	
	ests = c(data.frame(VarCorr(mod))[,'sdcor'], fixef(mod))
	g_ests = ranef(mod, whichel='Genus')[[1]]
	ints = confint(mod, parm=1:5, method='profile')
	g_sds = sqrt(as.numeric(attr(ranef(mod, whichel='Genus', condVar=T)[[1]], 'postVar')))
	g_ints = g_ests$'(Intercept)' + g_sds%*%t(c(-1.96,1.96))
	
	parm_ests_duke[i,j,,'est'] = ests
	parm_ests_duke[i,j,,c('low95','up95')] = ints
	genus_ests_duke[i,j,rownames(g_ests),'est'] = g_ests$'(Intercept)'
	genus_ests_duke[i,j,rownames(g_ests),c('low95','up95')] = g_ints
}}

save(modlist_duke, parm_ests_duke, genus_ests_duke, file='./Analysis/REML single variable FT models Duke Forest.RData')

pdf('./Analysis/Figures/REML single variable model effects Duke Forest.pdf', height=9, width=4)
layout(matrix(1:length(use_traits), ncol=1))
par(mar=c(2,13,2,.5))
for(i in use_traits){
	this_data = parm_ests_duke[i,,'b1',]
	xrange = range(this_data)
	plot(as.numeric(xvar_levels)~this_data[,'est'], xlim=xrange, axes=F, 
		xlab='', ylab='', lwd=1)
	abline(v=0, col='grey50', lty=2)
	segments(this_data[, 'low95'], as.numeric(xvar_levels), 
		this_data[, 'up95'], as.numeric(xvar_levels), lwd=1)
	axis(1)
	axis(2, at=1:length(env_vars), labels=xvarnames[xvar_levels], las=1)
	mtext(traitnames[i,'displayName'], 3, 0)
	box()
}
dev.off()




##############################################
### Decomposition of community trait variation into turnover versus intraspecific variation
## From Leps et al 2011, Ecography
source('./Analysis/Leps2011_trait_functions.R')

## Run code at the beginning of the previous section on individual trait models to get model_data with transformed traits
## note: we are not used scaled traits- those are only for multi-trait RDA 

## If performing analyses on foliose thalli only, make sure to subset out Usnea at the top

# Define seven focal traits for analysis
use_traits = c('Water_capacity','Thallus_thickness','STA','Cortex_thickness','Rhizine_length','Tot_chl_DW','Chla2b')

## Calculate genus-level mean traits
gen_means = aggregate(model_data[,use_traits], list(Genus=model_data$Genus), FUN=function(x) mean(x, na.rm=T))
gen_means[is.na(gen_means)] = NA
rownames(gen_means) = gen_means$Genus

gen_vars = aggregate(model_data[,use_traits], list(Genus=model_data$Genus), FUN=function(x) var(x, na.rm=T))

## Calculate species-level mean traits
sp_means = aggregate(model_data[,use_traits], list(Species=use_data$TaxonID), FUN=function (x) mean(x, na.rm=T))
sp_means[is.na(sp_means)] = NA
rownames(sp_means) = sp_means$Species

# Substitute genus-averaged traits for genus-only taxa
gen_only = rownames(sp_means)[which(taxa[rownames(sp_means),'TaxonConcept']=='genus')]
sp_means[gen_only,use_traits] = gen_means[taxa[gen_only,'Genus'],use_traits]

## Calculate specific average traits for each sample
SA_traits = aggregate(model_data[,use_traits], list(SampID=model_data$SampID), FUN=function(x) mean(x, na.rm=T))
SA_traits[is.na(SA_traits)] = NA
rownames(SA_traits) = paste('S', SA_traits$SampID, sep='')
SA_traits = SA_traits[,-1]
SA_traits = as.matrix(SA_traits)

## Calculate fixed average traits for each sample
sampXsp_abun = xtabs(~SampID+TaxonID, data=use_data)
rownames(sampXsp_abun) = paste('S',rownames(sampXsp_abun), sep='')
sampXgen_abun = xtabs(~SampID+Genus, data=use_data)
rownames(sampXgen_abun) = paste('S', rownames(sampXgen_abun), sep='')

# FA traits based on genus-level means
FA_traits_gen = matrix(NA, nrow=nrow(sampXgen_abun), ncol=length(use_traits))
rownames(FA_traits_gen) = rownames(sampXgen_abun); colnames(FA_traits_gen) = use_traits

for(j in use_traits){
	these_means = gen_means[,j]

	# Re-scale community data omitting thalli that have missing trait values
	use_comm = sampXgen_abun[,!is.na(these_means)]
	use_comm = use_comm/rowSums(use_comm)
	these_means = these_means[!is.na(these_means)]
	
	# Calculate community weighted averages using taxon means
	FA_traits_gen[,j] = use_comm %*% these_means
}

# Set fixed averages to NA whenever specific averages are NA
FA_traits_gen[is.na(SA_traits)] = NA

# FA trait based on species-level means
FA_traits = matrix(NA, nrow=nrow(sampXsp_abun), ncol=length(use_traits))
rownames(FA_traits) = rownames(sampXsp_abun); colnames(FA_traits) = use_traits

for(j in use_traits){
	these_means = sp_means[,j]

	# Re-scale community data omitting thalli that have missing trait values
	use_comm = sampXsp_abun[,!is.na(these_means)]
	use_comm = use_comm/rowSums(use_comm)
	these_means = these_means[!is.na(these_means)]
	
	# Calculate community weighted averages using taxon means
	FA_traits[,j] = use_comm %*% these_means
}

# Set fixed averages to NA whenever specific averages are NA
FA_traits[is.na(SA_traits)] = NA

# Calculate intraspecific variability effect for each sample
ISV_traits = SA_traits - FA_traits
ISV_traits_gen = SA_traits - FA_traits_gen

# Save trait means
#save(SA_traits, FA_traits, ISV_traits, FA_traits_gen, ISV_traits_gen, file='./Analysis/sample mean traits no Usnea.RData')
#save(SA_traits, FA_traits, ISV_traits, FA_traits_gen, ISV_traits_gen, file='./Analysis/sample mean traits with WHC STM.RData')


# Load previously saved matrices of sample-mean traits
load('./Analysis/sample mean traits no Usnea.RData')

# Check distributions
layout(matrix(1:(3*length(use_traits)), nrow=3, byrow=F))
par(mar=c(4, 4, 3.5, 1))
for(j in use_traits){
	hist(SA_traits[,j], main=j, xlab='SA')
	hist(FA_traits[,j], main='', xlab='FA')
	hist(ISV_traits[,j], main='', xlab='ISV')
}
for(j in use_traits){
	hist(SA_traits[,j], main=j, xlab='SA')
	hist(FA_traits_gen[,j], main='', xlab='FA')
	hist(ISV_traits_gen[,j], main='', xlab='ISV')
}

## Are species turnover and  intraspecific variability correlated?
cor(FA_traits, ISV_traits, use='complete.obs')
# Cortex_thickness moderately positively correlated (r = 0.37)
# Rhizine length strongly negatively correlated (r = -0.72)
# Chla2b moderately negatively correlted (r = -0.36)
cor(FA_traits_gen, ISV_traits_gen, use='complete.obs')
# Same correlation structure as for species-level averages

# Plot correlations
par(mfrow=c(2,4))
par(mar=c(4,4,1,1))
for(j in use_traits){
	plot(FA_traits[,j], ISV_traits[,j], xlab='Interspecific', ylab='Intraspecific', 
		main=traitnames[j,'displayName'])
}

### Trait decomposition analysis

# Define environmental data
Xdata = unique(model_data[,c('SampID', env_vars)])
rownames(Xdata) = paste('S', Xdata$SampID, sep='')
Xdata = Xdata[,-1]
Xdata = Xdata[rownames(SA_traits),]

# Define site covariate
year = unique(model_data[,c('SampID','Year')])
rownames(year) = paste('S', year$SampID, sep='')
year = year[rownames(SA_traits),]
year = factor(year[,-1])


## Decomposition of total trait variance
aov0 = matrix(NA, nrow=length(use_traits), ncol=3)
rownames(aov0) = use_traits; colnames(aov0) = c('Interspecific','Intraspecific','Covariation')

for(j in use_traits){
	decomp = trait.flex.anova(~1, specif.avg=SA_traits[,j], const.avg=FA_traits[,j])
	aov0[j,] = as.numeric(decomp$RelSumSq[1:3])
}

aov0_gen = matrix(NA, nrow=length(use_traits), ncol=3)
rownames(aov0_gen) = use_traits; colnames(aov0_gen) = c('Interspecific','Intraspecific','Covariation')

for(j in use_traits){
	decomp = trait.flex.anova(~1, specif.avg=SA_traits[,j], const.avg=FA_traits_gen[,j])
	aov0_gen[j,] = as.numeric(decomp$RelSumSq[1:3])
}

# Decomposition of total trait variance controlling for site
aov0_site = array(NA, dim=c(length(use_traits), 4, 5), 
	dimnames=list(Trait=use_traits, Component=c('Interspecific','Intraspecific','Covariation','Total'), Statistic=c('Var','Res','Tot','F','P')))

for(j in use_traits){
	this_aov = trait.flex.anova(~year, specif.avg=SA_traits[,j], const.avg=FA_traits[,j])

	aov0_site[j,,c('Var','Res','Tot')] = t(this_aov$SumSq)
	aov0_site[j,'Interspecific',c('F','P')] = as.numeric(this_aov$anova.turnover[1,4:5])
	aov0_site[j,'Intraspecific',c('F','P')] = as.numeric(this_aov$anova.diff[1,4:5])
	aov0_site[j,'Total',c('F','P')] = as.numeric(this_aov$anova.total[1,4:5])
} 

# Save results to table: array and file names change depending on which object is being saved.
site0_melt = melt(aov0_site[,c('Interspecific','Intraspecific','Total'),c('Var','F','P')])
site0_tab = cast(site0_melt, Trait~Component+Statistic)

write.csv(site0_tab, './Analysis/Figures/trait variation decomposition by site.csv', row.names=F)


## Plot bar chart of tot trait variance decomposition
bar_width=1
space_width=0.25*bar_width

bar_col = c('white','black','grey50')

pdf('./Analysis/Figures/Trait SS Decomposition no Usnea.pdf', height=6, width=7)
par(mar=c(0, 4, 8, 7))
plot.new()
plot.window(xlim=c(0, length(use_traits)*(bar_width+space_width)), ylim=c(min(aov0), 1))
abline(h=0)
for(i in 1:length(use_traits)){
	these_y = aov_heights(aov0[i,])
	rect(i*space_width+(i-1)*bar_width, these_y[1,], i*(space_width+bar_width), these_y[2,],
		col=bar_col)
}
axis(2, las=1)
mtext('Proportion Sum of Squares', 2, 2.8)
par(xpd=T)
text((1:length(use_traits))*(space_width+bar_width)-0.5*bar_width, 
	1, parse(text=traitnames[rownames(aov0),'exprName']), adj=c(-.05,0), srt=60)
legend('right',colnames(aov0), fill=bar_col, bty='n', inset=-.25)
par(xpd=F)
dev.off()

pdf('./Analysis/Figures/Trait SS Decomposition site residuals.pdf', height=6, width=7)
par(mar=c(0, 4, 8, 8))
plot.new()
plot.window(xlim=c(0, length(use_traits)*(bar_width+space_width)), ylim=c(-.4, 1))
abline(h=0)
for(i in 1:length(use_traits)){
	varcomp = aov0_site[i,c('Interspecific','Intraspecific','Covariation'),'Res']
	varcomp = varcomp / aov0_site[i,'Total','Tot']
	these_y = aov_heights(varcomp)
	rect(i*space_width+(i-1)*bar_width, these_y[1,], i*(space_width+bar_width), these_y[2,],
		col=bar_col)
}
axis(2, las=1)
mtext('Proportion Sum of Squares', 2, 2.8)
par(xpd=T)
text((1:length(use_traits))*(space_width+bar_width)-0.5*bar_width, 
	1, parse(text=traitnames[rownames(aov0),'exprName']), adj=c(-.05,0), srt=60)
legend('right',colnames(aov0), fill=bar_col, bty='n', inset=-.35)
par(xpd=F)
dev.off()


## Decomposition of environmental model
# All env variables
j='Water_capacity' # Define trait
aov_env = trait.flex.anova(~., specif.avg=SA_traits[,j], const.avg=FA_traits[,j], data=Xdata)

# Each variable separately
aov_env = vector('list', length(use_traits)*length(env_vars))
dim(aov_env) = c(length(use_traits),length(env_vars))
dimnames(aov_env) = list(use_traits, env_vars)

aovSS_array = array(NA, dim=c(length(use_traits), length(env_vars), 4, 5), 
	dimnames=list(Trait=use_traits, Predictor=env_vars, Component=c('Interspecific','Intraspecific','Covariation','Total'), Statistic=c('Var','Res','Tot','F','P')))

for(i in use_traits){
for(j in env_vars){
	this_aov = trait.flex.anova(~Xdata[,j], specif.avg = SA_traits[,i], const.avg=FA_traits[,i])
	aov_env[i,j][[1]] = this_aov
	aovSS_array[i,j,,c('Var','Res','Tot')] = t(this_aov$SumSq)
	aovSS_array[i,j,'Interspecific',c('F','P')] = as.numeric(this_aov$anova.turnover[1,4:5])
	aovSS_array[i,j,'Intraspecific',c('F','P')] = as.numeric(this_aov$anova.diff[1,4:5])
	aovSS_array[i,j,'Total',c('F','P')] = as.numeric(this_aov$anova.total[1,4:5])
}}

# Direction of effects of each env var on each component
effects_array = array(NA, dim=c(length(use_traits), length(env_vars), 3), 
	dimnames=list(Trait=use_traits, Predictor=env_vars, Component=c('Interspecific','Intraspecific','Total')))

for(i in use_traits){
for(j in env_vars){
	these_mods = lapply(list(ISV_traits[,i], FA_traits[,i], SA_traits[,i]), function(Y){
		lm(Y ~ Xdata[,j])
	})

	effects_array[i,j,] = sapply(these_mods, function(x) coef(x)[2])

}}


aovSS_melt = melt(aovSS_array[,,c('Interspecific','Intraspecific','Total'),c('Var','F','P')])
aovSS_tab = cast(aovSS_melt, Trait+Predictor~Component+Statistic)

write.csv(aovSS_tab, './Analysis/Figures/Environmental predictors of trait variation decomposition no Usnea.csv', row.names=F)
write.csv(aovSS_tab, './Analysis/Figures/Environmental predictors of trait variation decomposition with WHC STM.csv', row.names=F)

totSS_tab = cast(melt(aovSS_array[,'Light_high',c('Interspecific','Intraspecific','Total'),'Tot']), Trait~Component)
write.csv(totSS_tab, './Analysis/Figures/Total SS trait decomposition no Usnea.csv', row.names=F)
write.csv(totSS_tab, './Analysis/Figures/Total SS trait decomposition with WHC STM.csv', row.names=F)

## Each variable separately controlling for site
aov_env_site = vector('list', length(use_traits)*length(env_vars))
dim(aov_env_site) = c(length(use_traits),length(env_vars))
dimnames(aov_env_site) = list(use_traits, env_vars)

aovSS_site = array(NA, dim=c(length(use_traits), length(env_vars), 4, 5), 
	dimnames=list(Trait=use_traits, Predictor=env_vars, Component=c('Interspecific','Intraspecific','Covariation','Total'), Statistic=c('Var','Res','Tot','F','P')))

for(i in use_traits){
for(j in env_vars){
	this_aov = trait.flex.anova(~year+Xdata[,j], specif.avg = SA_traits[,i], const.avg=FA_traits[,i])
	
	aov_env_site[i,j][[1]] = this_aov
	aovSS_site[i,j,,c('Var','Res','Tot')] = t(this_aov$SumSq[2:4,])
	aovSS_site[i,j,'Interspecific',c('F','P')] = as.numeric(this_aov$anova.turnover[2,4:5])
	aovSS_site[i,j,'Intraspecific',c('F','P')] = as.numeric(this_aov$anova.diff[2,4:5])
	aovSS_site[i,j,'Total',c('F','P')] = as.numeric(this_aov$anova.total[2,4:5])
}}


# Save results in table
aovSS_site_melt = melt(aovSS_site[,,c('Interspecific','Intraspecific','Total'),c('Var','F','P')])
aovSS_site_tab = cast(aovSS_site_melt, Trait+Predictor~Component+Statistic)

site0_tab$Predictor = 'Site'
site0_tab = site0_tab[,colnames(aovSS_site_tab)]

aovSS_site_tab = rbind(aovSS_site_tab, site0_tab)

aovSS_site_tab = aovSS_site_tab[order(aovSS_site_tab$Trait, factor(aovSS_site_tab$Predictor, levels=c('Site',env_vars))),c(1,2,5,3,4,8,6,7,11,9,10)]

write.csv(aovSS_site_tab, './Analysis/Figures/Environmental predictors of trait variation decomposition control site.csv', row.names=F)
write.csv(aovSS_site_tab, './Analysis/Figures/Environmental predictors of trait variation decomposition control site with WHC STM.csv', row.names=F)


## Direction of effects of each env var on each component
effects_array_site = array(NA, dim=c(length(use_traits), length(env_vars), 3), 
	dimnames=list(Trait=use_traits, Predictor=env_vars, Component=c('Interspecific','Intraspecific','Total')))

for(i in use_traits){
for(j in env_vars){
	these_mods = lapply(list(ISV_traits[,i], FA_traits[,i], SA_traits[,i]), function(Y){
		lm(Y ~ year + Xdata[,j])
	})

	effects_array_site[i,j,] = sapply(these_mods, function(x) coef(x)[3])

}}


## Plot proportion of Total SS explained by each environmental variable
# use aovSS_array and effects_array for models without site
# use aovSS_site and effects_array_site for models with site

use_tab = cast(aovSS_melt, Trait+Predictor~Statistic, subset=Component=='Total'&Statistic=='Var')
props = aovSS_array[,,c('Interspecific','Intraspecific','Total'),'Var'] / aovSS_array[,,c('Interspecific','Intraspecific','Total'),'Tot']
props_mat = cast(melt(props), Trait+Predictor~Component)

props_sig = aovSS_array[,,c('Interspecific','Intraspecific','Total'),'P']
props_sig_mat = cast(melt(props_sig), Trait+Predictor~Component)
props_sig_mat[,3:5] = props_sig_mat[,3:5] < 0.05

effsign = cast(melt(effects_array), Trait+Predictor~Component)
effsign[,3:5] = effsign[,3:5]>0

image(t(as.matrix(props_mat[,c('Interspecific','Intraspecific','Total')])))

barcols = c('grey70','black')
signsymbs = c('-','+')

plot_props = props_mat

pdf('./Analysis/Figures/Sample mean trait variance explained by env predictors site.pdf', height=7, width=7)
layout(matrix(1:4, nrow=1), widths=c(0.4, 0.2, 0.2, 0.2))
par(oma=c(0,2,0,2))
par(mar=c(3,0,3,0.5))
plot.new()
plot.window(ylim=c(0,nrow(plot_props)), xlim=c(0,1))
axis(2, at=(1:nrow(plot_props))-0.5, labels=xvarnames_short[as.character(plot_props$Predictor)], las=1,
	tick=F, line=0, col=0, pos=1)
trait_text = traitnames[as.character(plot_props$Trait),'exprName']
trait_text[duplicated(trait_text)] = NA
for(i in 1:length(trait_text)) text(0.6, i+length(env_vars)-1.5, parse(text=trait_text[i]), pos=2, font=2)

plot.new()
plot.window(ylim=c(0,nrow(plot_props)), xlim=c(0,1))
for(i in 1:nrow(plot_props)){
	rect(0,i-1, plot_props[i,'Total'],i, col=barcols[props_sig_mat[i,'Total']+1])
	if(props_sig_mat[i,'Total']) text(plot_props[i,'Total']+.1, i-.5, signsymbs[effsign[i,'Total']+1])
}
abline(h=0:nrow(plot_props), col='grey50',lty=3)
abline(h=seq(0,nrow(plot_props), length(env_vars)))
axis(1, line=-1, las=3)
mtext('Total Trait Variation', 3, 0, cex=.8)

plot.new()
plot.window(ylim=c(0,nrow(plot_props)), xlim=c(0,1))
for(i in 1:nrow(plot_props)){
	rect(0,i-1, plot_props[i,'Intraspecific'],i, col=barcols[props_sig_mat[i,'Intraspecific']+1])
	if(props_sig_mat[i,'Intraspecific']) text(plot_props[i,'Intraspecific']+.1, i-.5, signsymbs[effsign[i,'Intraspecific']+1])
}
abline(h=0:nrow(plot_props), col='grey50',lty=3)
abline(h=seq(0,nrow(plot_props), length(env_vars)))
axis(1, line=-1, las=3)
mtext('Intraspecific', 3, 0, cex=.8)

plot.new()
plot.window(ylim=c(0,nrow(plot_props)), xlim=c(0,1))
for(i in 1:nrow(plot_props)){
	rect(0,i-1, plot_props[i,'Interspecific'],i, col=barcols[props_sig_mat[i,'Interspecific']+1])
	if(props_sig_mat[i,'Interspecific']) text(plot_props[i,'Interspecific']+.1, i-.5, signsymbs[effsign[i,'Interspecific']+1])
}
abline(h=0:nrow(plot_props), col='grey50',lty=3)
abline(h=seq(0,nrow(plot_props), length(env_vars)))
axis(1, line=-1, las=3)
mtext('Interspecific', 3, 0, cex=.8)
dev.off()

## Make smaller version of figure used in manuscript (Figure 4)
library(R.utils)

bestvars = apply(aovSS_array[,,'Total','Var'], 1, function(x) which(x==max(x)))
aovSS_best = sapply(use_traits, function(i) aovSS_array[i,bestvars[i],,], simplify='array')
names(dimnames(aovSS_best))[3] = 'Trait'
props = aovSS_best[c('Interspecific','Intraspecific','Total'),'Var',] / aovSS_best[c('Interspecific','Intraspecific','Total'),'Tot',]
props_mat = cast(melt(props), Trait~Component)
props_mat$Predictor = dimnames(aovSS_array)$Predictor[bestvars[as.character(props_mat$Trait)]]
rownames(props_mat) = props_mat$Trait

props_sig = aovSS_best[c('Interspecific','Intraspecific','Total'),'P',]
props_sig_mat = cast(melt(props_sig), Trait~Component)
props_sig_mat[,2:4] = props_sig_mat[,2:4] < 0.05
rownames(props_sig_mat) = props_sig_mat$Trait

effects_best = sapply(use_traits, function(i) effects_array[i,bestvars[i],], simplify='array')
effsign = t(effects_best)>0
effsign = effsign[as.character(props_mat$Trait),]

plot_props = props_mat[use_traits[7:1],]
props_sig_mat = props_sig_mat[use_traits[7:1],]
effsign = effsign[use_traits[7:1],]

pdf('./Analysis/Figures/Sample mean trait variance explained by best env predictors site.pdf', height=2, width=7)

# Set up plot layout: 4 columns
layout(matrix(1:4, nrow=1), widths=c(0.4, 0.2, 0.2, 0.2))
par(oma=c(0,2,0,2))
par(mar=c(3,0,3,0.5))
plot.new()
plot.window(ylim=c(0,nrow(plot_props)), xlim=c(0,1))

# Add y-axis labels in 1st plot
lab1 = parse(text=traitnames[as.character(plot_props$Trait),'exprName']) # Trait names
lab2 = xvarnames_short[as.character(plot_props$Predictor)] # Predictors
N = nrow(plot_props)
yvals = (1:N)-0.5 # locations to put labels
mtext('Predictor', 3, 0, adj=1, cex=0.8)
axis(2, at=yvals, labels=lab2, las=1,
	tick=F, col=0, pos=1.1)
mtext('Trait', 3, 0, adj=0.4, cex=0.8)
text(0.65, yvals, parse(text=lab1), pos=2, font=2)

for(comp in c('Total','Intraspecific','Interspecific')){
	# Set up plot window
	plot.new()
	plot.window(ylim=c(0,nrow(plot_props)), xlim=c(0,1))
	
	# Add lines dividing bars
	abline(h=0:N, col='grey70',lty=3, lwd=.8, lend=2)
	
	# Add rectangles for variance explained
	rect(0, 0:(N-1), plot_props[,comp], 1:N, col=barcols[props_sig_mat[,comp]+1], border=NA)
	
	# Add symbol for direction of effect
	signtext = signsymbs[effsign[,comp]+1]
	signtext[!props_sig_mat[,comp]] = NA
	text(plot_props[,comp]+.1, yvals, signtext)

	# Add x-axis and title
	axis(1, line=0, las=2, srt=180)
	mtext(comp, 3, 0, cex=.8)
}

dev.off()

# Compare water capacity
wtraits = c('Water_capacity','WHC','STA','STM')
plot_props = props_mat[wtraits,]
props_sig_mat = props_sig_mat[wtraits,]
effsign = effsign[wtraits,]

pdf('./Analysis/Figures/Sample mean trait variance explained by best env predictors site compare water capacity.pdf', height=1.75, width=7)

# Set up plot layout: 4 columns
layout(matrix(1:4, nrow=1), widths=c(0.4, 0.2, 0.2, 0.2))
par(oma=c(0,2,0,2))
par(mar=c(3,0,3,0.5))
plot.new()
plot.window(ylim=c(0,nrow(plot_props)), xlim=c(0,1))

# Add y-axis labels in 1st plot
lab1 = parse(text=traitnames[as.character(plot_props$Trait),'exprName']) # Trait names
lab2 = xvarnames_short[as.character(plot_props$Predictor)] # Predictors
N = nrow(plot_props)
yvals = (1:N)-0.5 # locations to put labels
mtext('Predictor', 3, 0, adj=1, cex=0.8)
axis(2, at=yvals, labels=lab2, las=1,
	tick=F, col=0, pos=1.1)
mtext('Trait', 3, 0, adj=0.4, cex=0.8)
text(0.65, yvals, parse(text=lab1), pos=2, font=2)

for(comp in c('Total','Intraspecific','Interspecific')){
	# Set up plot window
	plot.new()
	plot.window(ylim=c(0,nrow(plot_props)), xlim=c(0,1))
	
	# Add lines dividing bars
	abline(h=0:N, col='grey70',lty=3, lwd=.8, lend=2)
	
	# Add rectangles for variance explained
	rect(0, 0:(N-1), plot_props[,comp], 1:N, col=barcols[props_sig_mat[,comp]+1], border=NA)
	
	# Add symbol for direction of effect
	signtext = signsymbs[effsign[,comp]+1]
	signtext[!props_sig_mat[,comp]] = NA
	text(plot_props[,comp]+.1, yvals, signtext)

	# Add x-axis and title
	axis(1, line=0, las=2, srt=180)
	mtext(comp, 3, 0, cex=.8)
}

dev.off()


################################################################################
### RDA

library(vegan)

# Define focal traits for analysis
use_traits = c('Water_capacity','Thallus_thickness','STA','Cortex_thickness','Rhizine_length','Tot_chl_DW','Chla2b')


# Scale trait data before calculating sample level means so that some traits are not more stongly weighted than others in the ordination 
scaled_data = model_data

scaled_data[,use_traits] = scale(scaled_data[,use_traits], center=T, scale=T)

# Re-calculate species mean traits
sp_means = aggregate(scaled_data[,use_traits], list(Species=use_data$TaxonID), FUN=function (x) mean(x, na.rm=T))
sp_means[is.na(sp_means)] = NA
rownames(sp_means) = sp_means$Species
gen_means = aggregate(scaled_data[,use_traits], list(Genus=scaled_data$Genus), FUN=function(x) mean(x, na.rm=T))
gen_means[is.na(gen_means)] = NA
rownames(gen_means) = gen_means$Genus
gen_only = rownames(sp_means)[which(taxa[rownames(sp_means),'TaxonConcept']=='genus')]
sp_means[gen_only,use_traits] = gen_means[taxa[gen_only,'Genus'],use_traits]


# This code relies on matrices of trait means across samples calculated in the previous section:
# Xdata = sample environment data that matches order  or trait mean matrices
# year = factor indicating which site samples come from

### RDA on thallus traits

## Unconstrained ordination of thalli to identify primary covarition of traits
## Followed by fitting env variables

# Impute missing trait values from species-level means (this may bias toward and effect of genus)
Ydata = scaled_data[,use_traits]
for(i in use_traits){
	missing = is.na(Ydata[,i]) # Find missing data
	Ydata[missing,i] = sp_means[use_data$TaxonID[missing],i] # Look up trait values for the species that are missing

	# If data are still missing, get them from genus-level means
	missing = is.na(Ydata[,i])
	if(sum(missing)>0){
		Ydata[missing,i] = gen_means[scaled_data$Genus[missing],i]
	}
}
# note: rhizine length will still be NA for Usnea
# two analyses: only foliose taxa with all traits, all taxa without rhizines


thalli_ord0 = rda(Ydata[,colnames(Ydata)!='Rhizine_length'])
foliose_ord0 = rda(Ydata[scaled_data$Genus!='Usnea',])

# Save trait loadings and eigenvalues
score_mat = scores(thalli_ord0, 1:6)$species
eigs = eigenvals(thalli_ord0)
eigs_pct = eigs/sum(eigs)
score_mat = rbind(score_mat, eigs, eigs_pct)
write.csv(score_mat, './Analysis/Figures/unconstrained ind thalli PCA trait loadings.csv')
score_mat = scores(foliose_ord0, 1:7)$species
eigs = eigenvals(foliose_ord0)
eigs_pct = eigs/sum(eigs)
score_mat = rbind(score_mat, eigs, eigs_pct)
write.csv(score_mat, './Analysis/Figures/unconstrained ind thalli PCA trait loadings foliose only.csv')


Xdata_foliose = scaled_data[scaled_data$Genus!='Usnea',env_vars]
Xdata_thalli = scaled_data[,env_vars]

year_thalli  = factor(scaled_data$Year)
year_foliose = factor(scaled_data[scaled_data$Genus!='Usnea','Year'])


use_ord = foliose_ord0
use_Xdata = Xdata_foliose
use_year = year_foliose
ord_sum = summary(use_ord)$cont$importance

scores(use_ord, 1:4, display='species')
ord_sum

# Colored by site
use_col = mycolor[c(2,8)]

pdf('./Analysis/Figures/ordination of thalli foliose only color by site.pdf', height=5, width=10)
layout(matrix(1:3, nrow=1), widths=c(0.4,0.4,0.1))
par(mar=c(4,4,1,1))

ev = envfit(use_ord, use_Xdata, choices=1:2)
pcts = paste(format(ord_sum[2,1:2]*100, digits=2), '%', sep='')
op = ordiplot(use_ord, c(1,2), type='n', xlab=paste('PC1 (',pcts[1], ')', sep=''), ylab=paste('PC2 (',pcts[2], ')', sep=''))
points(op, 'sites', pch=21, bg=use_col[factor(use_year)], cex=0.8)
ordiellipse(op, groups=use_year, kind='sd', draw='polygon', col=use_col[1], show.groups=2013)
ordiellipse(op, groups=use_year, kind='sd', draw='polygon', col=use_col[2], show.groups=2014)
par(xpd=T)
text(op, 'species', labels=traitnames[names(use_ord$colsum), 'shortName'])
plot(ev, labels=xvarnames_short[colnames(Xdata)], col='grey20', add=T, p.max=0.05, arrow.mul=3, cex=0.9)
par(xpd=F)

ev = envfit(use_ord, use_Xdata, choices=3:4)
pcts = paste(format(ord_sum[2,3:4]*100, digits=2), '%', sep='')
op = ordiplot(use_ord, c(3,4), type='n', xlab=paste('PC3 (',pcts[1], ')', sep=''), ylab=paste('PC4 (',pcts[2], ')', sep=''))
points(op, 'sites', pch=21, bg=use_col[factor(use_year)], cex=0.8)
ordiellipse(op, groups=use_year, kind='sd', draw='polygon', col=use_col[1], show.groups=2013)
ordiellipse(op, groups=use_year, kind='sd', draw='polygon', col=use_col[2], show.groups=2014)
par(xpd=T)
text(op, 'species', labels=traitnames[names(use_ord$colsum), 'shortName'])
plot(ev, labels=xvarnames_short[colnames(Xdata)], col='grey20', add=T, p.max=0.05, arrow.mul=7, cex=0.9)
par(xpd=F)

plot.new()
par(mar=c(0,0,0,0))
legend('right', c('Water Dog','Duke'), col=use_col, pch=16, bty='n')

dev.off()

# Colored by genus
genera = levels(factor(model_data$Genus)) 
genera = genera[genera!='Usnea']
use_col = colorRampPalette(mycolor)(length(genera)); names(use_col) = genera
color_fact = factor(model_data$Genus)
color_fact = color_fact[model_data$Genus!='Usnea']

pdf('./Analysis/Figures/ordination of thalli foliose only color by Genus.pdf', height=5, width=10)
layout(matrix(1:3, nrow=1), widths=c(0.4,0.4,0.1))
par(mar=c(4,4,1,1))

pcts = paste(format(ord_sum[2,1:2]*100, digits=2), '%', sep='')
op = ordiplot(use_ord, c(1,2), type='n', xlab=paste('PC1 (',pcts[1], ')', sep=''), ylab=paste('PC2 (',pcts[2], ')', sep=''))
points(op, 'sites', pch=21, bg='white', cex=0.8) #use_col[color_fact]
for(g in genera){
	ordihull(op, groups=color_fact, kind='sd', draw='polygon', col=use_col[g], show.groups=g, alpha=50)
	ordihull(op, groups=color_fact, kind='sd', draw='line',lwd=2, col=use_col[g], show.groups=g)
}
par(xpd=T)
text(op, 'species', labels=traitnames[names(use_ord$colsum), 'shortName'])
par(xpd=F)

pcts = paste(format(ord_sum[2,3:4]*100, digits=2), '%', sep='')
op = ordiplot(use_ord, c(3,4), type='n', xlab=paste('PC3 (',pcts[1], ')', sep=''), ylab=paste('PC4 (',pcts[2], ')', sep=''))
points(op, 'sites', pch=21, bg='white', cex=0.8) #use_col[color_fact]
for(g in genera){
	ordihull(op, groups=color_fact, kind='sd', draw='polygon', col=use_col[g], show.groups=g, alpha=50)
	ordihull(op, groups=color_fact, kind='sd', draw='line',lwd=2, col=use_col[g], show.groups=g)
}
par(xpd=T)
text(op, 'species', labels=traitnames[names(use_ord$colsum), 'shortName'])
par(xpd=F)
plot.new()
par(mar=c(0,0,0,0))
legend('right', genera, lwd=2, col=use_col)

dev.off()


## Variation partitioning of traits of thalli among ENV, GENUS, and SITE

# All thalli but no rhizine length
vp_thalli = varpart(Ydata[,colnames(Ydata)!='Rhizine_length'], Xdata_thalli, ~Genus, ~factor(Year), data=scaled_data)

# Foliose thalli only, inclduing rhizine length
data_foliose = subset(scaled_data, Genus!='Usnea')
vp_foliose = varpart(Ydata[model_data$Genus!='Usnea',], Xdata_foliose, ~Genus, ~factor(Year),data=data_foliose)

# X1=ENV, X2=Genus, X3=Site
vp_thalli$part
vp_foliose$part

# Fraction explained by site uniquely and shared between site and genus is very small so we can ignore site.
svg('./Analysis/Figures/individual traits variance partition.svg', height=4, width=4)
par(mar=c(0,0,0,0))
plot(vp_thalli, digits=2)
usr = par('usr')
xlen = usr[2]-usr[1]
ylen = usr[4]-usr[3]
text(usr[1]+(.32*xlen),usr[4]-(.1*ylen), 'Environment')
text(usr[1]+(.67*xlen),usr[4]-(.1*ylen), 'Genus')
text(usr[1]+(.5*xlen),usr[3]+(.1*ylen), 'Site')
dev.off()

svg('./Analysis/Figures/individual traits variance partition no Usnea.svg', height=4, width=4)
par(mar=c(0,0,0,0))
plot(vp_foliose, digits=2)
usr = par('usr')
xlen = usr[2]-usr[1]
ylen = usr[4]-usr[3]
text(usr[1]+(.32*xlen),usr[4]-(.1*ylen), 'Environment')
text(usr[1]+(.67*xlen),usr[4]-(.1*ylen), 'Genus')
text(usr[1]+(.5*xlen),usr[3]+(.1*ylen), 'Site')
dev.off()

# Partition only between Env and Genus
vp_thalli = varpart(Ydata[,colnames(Ydata)!='Rhizine_length'], Xdata_thalli, ~Genus, data=model_data)
vp_foliose = varpart(Ydata[model_data$Genus!='Usnea',], Xdata_foliose, ~Genus, data=data_foliose)

vp_df = data.frame(Thalli = rep(c('All','Foliose'), each=3), Predictor=rep(c('Env','Shared','Genus'), 2))
vp_df$R2 = c(vp_thalli$part$indfract[1:3,'Adj.R.squared'],vp_foliose$part$indfract[1:3,'Adj.R.squared']) 

write.csv(vp_df, './Analysis/Figures/thalli RDA varpart.csv', row.names=F)


## Constrained ordination of thallus traits using all env variables
rda_thalli_env = rda(Ydata[colnames(Ydata)!='Rhizine_length'], Xdata_thalli, scale=F)
anova(rda_thalli_env)

year_thalli = scaled_data$Year
year_foliose = year_thalli[scaled_data$Genus!='Usnea']

rda_thalli_env_site = rda(Ydata[colnames(Ydata)!='Rhizine_length'], Xdata_thalli, year_thalli, scale=F)
anova(rda_thalli_env_site)

RsquareAdj(rda_thalli_env)
RsquareAdj(rda_thalli_env_site)


## Constrained ordination of thallus traits: make table of effects
rda_thalli_array = array(NA, dim=c(length(env_vars), 6, 2, 2), 
	dimnames=list(Predictor=env_vars, Statistic=c('R2','Var','Res','Tot','F','P'), Thalli=c('All','Foliose'), Condition=c('None','Site')))
rda_thalli_scores = array(NA, dim=c(length(env_vars), length(use_traits),2,2), 
	dimnames=list(Predictor=env_vars, Trait=use_traits, Thalli=c('All','Foliose'), Condition=c('None','Site')))


for(i in env_vars){
	ord_thalli = rda(Ydata[,colnames(Ydata)!='Rhizine_length'], Xdata_thalli[,i], scale=F)
	ord_thalli_site = rda(Ydata[,colnames(Ydata)!='Rhizine_length'], Xdata_thalli[,i], year_thalli, scale=F)
	ord_foliose = rda(Ydata[model_data$Genus!='Usnea',], Xdata_foliose[,i], scale=F)
	ord_foliose_site = rda(Ydata[model_data$Genus!='Usnea',], Xdata_foliose[,i], year_foliose, scale=F)

	rda_thalli_scores[i,rownames(scores(ord_thalli, 1, 'species')),'All','None'] = scores(ord_thalli, 1, 'species')
	rda_thalli_scores[i,rownames(scores(ord_thalli_site, 1, 'species')),'All','Site'] = scores(ord_thalli_site, 1, 'species')
	rda_thalli_scores[i,rownames(scores(ord_foliose, 1, 'species')),'Foliose','None']  = scores(ord_foliose, 1, 'species')
	rda_thalli_scores[i,rownames(scores(ord_foliose_site, 1, 'species')),'Foliose','Site'] = scores(ord_foliose_site, 1, 'species')

	# Variance components from eigenvalues
	rda_thalli_array[i,c('Var','Res','Tot'),,'None'] = sapply(list(ord_thalli, ord_foliose), function(x){
		eigs = eigenvals(x)
		c(eigs[1], sum(eigs)-eigs[1], sum(eigs))
	})
	rda_thalli_array[i,c('Var','Res','Tot'),,'Site'] = sapply(list(ord_thalli_site, ord_foliose_site), function(x){
		eigs = eigenvals(x)
		c(eigs[1], sum(eigs)-eigs[1], sum(eigs))
	})

	# Adjusted R2 of environmental effect
	rda_thalli_array[i,'R2',,] = 	sapply(list(ord_thalli, ord_foliose, ord_thalli_site, ord_foliose_site), function(x) RsquareAdj(x)$adj.r.squared)
	
	# Significance from F-statistic
	aov_thalli = anova(ord_thalli)
	aov_thalli_site = anova(ord_thalli_site)
	aov_foliose = anova(ord_foliose)
	aov_foliose_site = anova(ord_foliose_site)

	rda_thalli_array[i,c('F','P'),,'None'] = sapply(list(aov_thalli, aov_foliose), function(x){
		as.numeric(x[1,3:4])
	})
	rda_thalli_array[i,c('F','P'),,'Site'] = sapply(list(aov_thalli_site, aov_foliose_site), function(x){
		as.numeric(x[1,3:4])
	})

}

rda_thalli_array[,'Tot',,]
rda_thalli_array[,'R2',,]
rda_thalli_array[,c('R2','P'),'All',]

thalli_melt = melt(rda_thalli_array)
thalli_tab = cast(thalli_melt, Thalli+Predictor+Condition~Statistic)

write.csv(thalli_tab, './Analysis/Figures/Evironmental predictors of multi-trait thalli.csv', row.names=F)

thalli_scores_melt = melt(rda_thalli_scores)
thalli_scores_tab = cast(thalli_scores_melt,Thalli+Predictor+Condition~Trait)
thalli_scores_tab$'% Variation' = thalli_tab$R2

write.csv(thalli_scores_tab, './Analysis/Figures/Evironmental predictors of multi-trait thalli scores.csv', row.names=F)

## Plot variance explained by each env variable
## Manuscript Figure 3
bp_tab = t(rda_thalli_array[,'R2','All',])
bp_tab[bp_tab <0] = 0
colnames(bp_tab) = xvarnames_short[colnames(bp_tab)]
All = lapply(list(rda_thalli_env, rda_thalli_env_site), function(x) RsquareAdj(x)$adj.r.squared)
bp_tab = cbind(bp_tab, All)

svg('./Analysis/Figures/individual traits RDA by env and site.svg', height=4, width=3)
par(mar=c(7,4,1,1))
barplot(bp_tab, las=3, ylim=c(0, .1), col=c('grey80','black'), ylab='Explained Variance')
dev.off()



### RDA on sample-level means

## We re-calculate these matrices here because we are using scaled traits and we will use the species averages to impute trait values for samples where these are missing
SA_traits = aggregate(scaled_data[,use_traits], list(SampID=scaled_data$SampID), FUN=function(x) mean(x, na.rm=T))
SA_traits[is.na(SA_traits)] = NA
rownames(SA_traits) = paste('S', SA_traits$SampID, sep='')
SA_traits = SA_traits[,-1]
SA_traits = as.matrix(SA_traits)

# will calculate FA traits based on abundance of each species
# can't use sampXsp_thalli calculated in examine_community_structure.R b/c it is presence/absence
sampXsp_abun = xtabs(~SampID+TaxonID, data=use_data)
rownames(sampXsp_abun) = paste('S',rownames(sampXsp_abun), sep='')

FA_traits = matrix(NA, nrow=nrow(sampXsp_abun), ncol=length(use_traits))
rownames(FA_traits) = rownames(sampXsp_abun); colnames(FA_traits) = use_traits

for(j in use_traits){
	these_means = sp_means[,j]

	# Re-scale community data omitting thalli that have missing trait values
	use_comm = sampXsp_abun[,!is.na(these_means)]
	use_comm = use_comm/rowSums(use_comm)
	these_means = these_means[!is.na(these_means)]
	
	# Calculate community weighted averages using taxon means
	FA_traits[,j] = use_comm %*% these_means
}

# Impute missing trait values from species averages (because RDA does not allow missing values)
SA_traits[is.na(SA_traits)] = FA_traits[is.na(SA_traits)] # 3 cases

# Calculate intraspecific variability effect for each sample
ISV_traits = SA_traits - FA_traits # Note that this will be 0 for cases when SA was imputed from FA

# Save sample mean traits
save(SA_traits, FA_traits, ISV_traits, file='./Analysis/sample mean traits.RData')

## Unconstrained ordination of sample mean traits to identify primary covariation of traits
## Followed by environmental fitting

SA_ord = rda(SA_traits)
FA_ord = rda(FA_traits)
ISV_ord = rda(ISV_traits)

use_col = mycolor[c(2,8)]

summary(SA_ord)

use_ord = FA_ord
ord_sum = summary(use_ord)$cont$importance


## Define environmental data matrix to match mean trait matrices
Xdata = unique(model_data[,c('SampID', env_vars)])
rownames(Xdata) = paste('S', Xdata$SampID, sep='')
Xdata = Xdata[,-1]
Xdata = Xdata[rownames(SA_traits),]

year = unique(model_data[,c('SampID','Year')])
rownames(year) = paste('S', year$SampID, sep='')
year = year[rownames(SA_traits),]
year = factor(year[,-1])


## Fit env variables to unconstrained ordination
ev = envfit(use_ord, Xdata, choices=1:4)

pdf('./Analysis/Figures/ordination of FA.pdf', height=5, width=10)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))

ev = envfit(use_ord, Xdata, choices=1:2)
pcts = paste(format(ord_sum[2,1:2]*100, digits=2), '%', sep='')
op = ordiplot(use_ord, c(1,2), type='n', xlab=paste('PC1 (',pcts[1], ')', sep=''), ylab=paste('PC2 (',pcts[2], ')', sep=''))
ordiellipse(op, groups=year, kind='sd', draw='polygon', col=use_col[1], show.groups=2013)
ordiellipse(op, groups=year, kind='sd', draw='polygon', col=use_col[2], show.groups=2014)
points(op, 'sites', pch=21, bg=use_col[factor(year)], cex=0.8)
text(op, 'species', labels=traitnames[use_traits, 'shortName'])
plot(ev, labels=xvarnames_short[colnames(Xdata)], col='grey20', add=T, p.max=0.05, arrow.mul=1, cex=0.9)

ev = envfit(use_ord, Xdata, choices=3:4)
pcts = paste(format(ord_sum[2,3:4]*100, digits=2), '%', sep='')
op = ordiplot(use_ord, c(3,4), type='n', xlab=paste('PC3 (',pcts[1], ')', sep=''), ylab=paste('PC4 (',pcts[2], ')', sep=''))
ordiellipse(op, groups=year, kind='sd', draw='polygon', col=use_col[1], show.groups=2013)
ordiellipse(op, groups=year, kind='sd', draw='polygon', col=use_col[2], show.groups=2014)
points(op, 'sites', pch=21, bg=use_col[factor(year)], cex=0.8)
text(op, 'species', labels=traitnames[use_traits, 'shortName'])
plot(ev, labels=xvarnames_short[colnames(Xdata)], col='grey20', add=T, p.max=0.05, arrow.mul=1, cex=0.9)

dev.off()

## RDA using all variables
SA_ord_env = rda(SA_traits, Xdata, scale=F)
anova(SA_ord_env)
RsquareAdj(SA_ord_env)

SA_ord_env_site = rda(SA_traits, Xdata, as.numeric(year)-1, scale=F)
anova(SA_ord_env_site)
RsquareAdj(SA_ord_env_site)

## Regress on each environmental variables separately

rda_env = vector('list', length(env_vars)*3)
dim(rda_env) = c(length(env_vars), 3)
dimnames(rda_env) = list(env_vars, c('Intraspecific','Interspecific','Total'))

rda_array = array(NA, dim=c(length(env_vars), 4, 6), 
	dimnames=list(Predictor=env_vars, Component=c('Interspecific','Intraspecific','Covariation','Total'), Statistic=c('R2','Var','Res','Tot','F','P')))
rda_scores = array(NA, dim=c(length(env_vars), length(use_traits), 3), 
	dimnames=list(Predictor=env_vars, Trait=use_traits, Component=c('Interspecific','Intraspecific','Total')))

for(i in env_vars){
	ord_SA = rda(SA_traits, Xdata[,i], scale=F)
	ord_FA = rda(FA_traits, Xdata[,i], scale=F)
	ord_ISV = rda(ISV_traits, Xdata[,i], scale=F)

	rda_env[i,] = list(ord_FA, ord_ISV, ord_SA)
	rda_scores[i,,'Total'] = scores(ord_SA, 1, 'species') 
	rda_scores[i,,'Interspecific'] = scores(ord_FA, 1, 'species')
	rda_scores[i,,'Intraspecific'] = scores(ord_ISV, 1, 'species')

	aov_SA = anova(ord_SA)
	aov_FA = anova(ord_FA)
	aov_ISV = anova(ord_ISV)

	# Variance components from eigenvalues
	rda_array[i,c('Intraspecific', 'Interspecific','Total'),c('Var','Res','Tot')] = t(sapply(list(ord_ISV,ord_FA,ord_SA), function(x){
		eigs = eigenvals(x)
		c(eigs[1], sum(eigs)-eigs[1], sum(eigs))
	}))

	# Variance explained using adjusted R2
	rda_array[i,c('Intraspecific', 'Interspecific','Total'),'R2'] = sapply(list(ord_ISV, ord_FA, ord_SA), function(x) RsquareAdj(x)$adj.r.squared)

	# Covariation as the difference
	rda_array[i,'Covariation',] = rda_array[i,'Total',] - rda_array[i,'Interspecific',] - rda_array[i,'Intraspecific',]

	rda_array[i,c('Intraspecific','Interspecific','Total'), c('F','P')] = t(sapply(list(aov_ISV, aov_FA, aov_SA), function(x){
		as.numeric(x[1,3:4])
	}))
}


rda_melt = melt(rda_array[,c('Interspecific','Intraspecific','Total'),c('R2','Var','F','P')])
rda_tab = cast(rda_melt, Component+Predictor~Statistic)

write.csv(rda_tab, './Analysis/Figures/Environmental predictors of multi-trait variance decomposition.csv', row.names=F )

# Variance explained based on adjusted R2
rda_array[,'Total',c('R2','P')]

# Proportion of total trait variation that can be explained by environmentally constrained axis
# Based on eigenvalue decomposition
rda_total = melt(rda_array[,'Total',])
rda_total_tab = cast(rda_total, Predictor~Statistic)
rda_total_tab$Proportion = rda_total_tab$Var / rda_total_tab$Tot

bar_width=1
space_width=0.25*bar_width

bar_col = c('white','black','grey50')

pdf('./Analysis/Figures/Multi-trait Eigenvalue Decomposition.pdf', height=6, width=7)
par(mar=c(11, 4, 0, 7))
plot.new()
plot.window(xlim=c(0, (length(env_vars) + 1)*(bar_width+space_width)), ylim=c(-.15, 1))
abline(h=0)
for(i in 1:length(env_vars)){
	# Explained variance for each component
	varcomp = rda_array[i,c('Interspecific','Intraspecific','Covariation'),'Var'] 
	varcomp = varcomp / rda_array[i,'Total','Tot'] # Scale by total variance
	
	these_y = aov_heights(varcomp)
	rect(i*space_width+(i-1)*bar_width, these_y[1,], i*(space_width+bar_width), these_y[2,],
		col=bar_col)
}

varcomp = rda_array[1,c('Interspecific','Intraspecific','Covariation'),'Tot'] / rda_array[1,'Total','Tot']
these_y = aov_heights(varcomp)
rect((i+1)*space_width+(i)*bar_width, these_y[1,], (i+1)*(space_width+bar_width), these_y[2,], col=bar_col)

axis(2, las=3)
mtext('Proportion Variance', 2, 2.8)
par(xpd=T)
text((1:length(env_vars))*(space_width+bar_width)-0.5*bar_width, 
	-.15, xvarnames[env_vars], adj=c(1,0.5), srt=90)
legend('right',colnames(aov0), fill=bar_col, bty='n', inset=-.25)
text((i+1)*(space_width+bar_width)-0.5*bar_width, -.15, 'Total', adj=c(1,0.5), srt=90)
par(xpd=F)
dev.off()

svg('./Analysis/Figures/Multi-trait RDA R2.svg', height=5, width=4)
par(mar=c(8, 4, 0.5,1))
plot.new()
plot.window(xlim=c(0, length(env_vars)*(bar_width+space_width)), ylim=c(-.15, 0.5))
abline(h=0)
for(i in 1:length(env_vars)){
	# Explained variance for each component
	varcomp = rda_array[i,c('Interspecific','Intraspecific','Covariation'),'R2'] 
	
	these_y = aov_heights(varcomp)
	rect(i*space_width+(i-1)*bar_width, these_y[1,], i*(space_width+bar_width), these_y[2,],
		col=bar_col)
}

axis(2, las=3)
mtext('Proportion Variance', 2, 2.8)
par(xpd=T)
text((1:length(env_vars))*(space_width+bar_width)-0.5*bar_width, 
	-.15, xvarnames_short[env_vars], adj=c(1,0.5), srt=90)
legend('top',c('Interspecific','Intraspecific','Covariation'), fill=bar_col, bty='n')
par(xpd=F)
dev.off()



## Control for site
rda_site = vector('list', length(env_vars)*3)
dim(rda_site) = c(length(env_vars), 3)
dimnames(rda_site) = list(env_vars, c('Intraspecific','Interspecific','Total'))

rda_site_array = array(NA, dim=c(length(env_vars)+1, 4, 6), 
	dimnames=list(Predictor=c('Site',env_vars), Component=c('Interspecific','Intraspecific','Covariation','Total'), Statistic=c('R2','Var','Res','Tot','F','P')))

for(i in env_vars){
	ord_SA = rda(SA_traits, Xdata[,i], as.numeric(year)-1, scale=F)
	ord_FA = rda(FA_traits, Xdata[,i], as.numeric(year)-1, scale=F)
	ord_ISV = rda(ISV_traits, Xdata[,i], as.numeric(year)-1, scale=F)

	rda_site[i,] = list(ord_FA, ord_ISV, ord_SA)

	aov_SA = anova(ord_SA)
	aov_FA = anova(ord_FA)
	aov_ISV = anova(ord_ISV)

	# Variance components from eigenvalues
	rda_site_array[i,c('Intraspecific', 'Interspecific','Total'),c('Var','Res','Tot')] = t(sapply(list(ord_ISV,ord_FA,ord_SA), function(x){
		eigs = eigenvals(x)
		c(eigs[1], sum(eigs)-eigs[1], sum(eigs))
	}))

	# Variance explained using adjusted R2
	rda_site_array[i,c('Intraspecific', 'Interspecific','Total'),'R2'] = sapply(list(ord_ISV, ord_FA, ord_SA), function(x) RsquareAdj(x)$adj.r.squared)

	# Covariation as the difference
	rda_site_array[i,'Covariation',] = rda_site_array[i,'Total',] - rda_site_array[i,'Interspecific',] - rda_site_array[i,'Intraspecific',]

	rda_site_array[i,c('Intraspecific','Interspecific','Total'), c('F','P')] = t(sapply(list(aov_ISV, aov_FA, aov_SA), function(x){
		as.numeric(x[1,3:4])
	}))
}

# Do site by itself
ord_SA = rda(SA_traits, as.numeric(year)-1, scale=F)
ord_FA = rda(FA_traits, as.numeric(year)-1, scale=F)
ord_ISV = rda(ISV_traits, as.numeric(year)-1, scale=F)
aov_SA = anova(ord_SA)
aov_FA = anova(ord_FA)
aov_ISV = anova(ord_ISV)
rda_site_array['Site',c('Intraspecific', 'Interspecific','Total'),c('Var','Res','Tot')] = t(sapply(list(ord_ISV,ord_FA,ord_SA), function(x){
	eigs = eigenvals(x)
	c(eigs[1], sum(eigs)-eigs[1], sum(eigs))
}))
rda_site_array['Site','Covariation',] = rda_site_array[i,'Total',] - rda_site_array[i,'Interspecific',] - rda_site_array[i,'Intraspecific',]
rda_site_array['Site',c('Intraspecific','Interspecific','Total'), c('F','P')] = t(sapply(list(aov_ISV, aov_FA, aov_SA), function(x){
	as.numeric(x[1,3:4])
}))
rda_site_array['Site',c('Intraspecific','Interspecific','Total'),'R2'] = sapply(list(ord_ISV, ord_FA, ord_SA), function(x) RsquareAdj(x)$adj.r.squared)

rda_site_melt = melt(rda_site_array[,c('Interspecific','Intraspecific','Total'),c('Var','F','P')])
rda_site_tab = cast(rda_site_melt, Predictor~Component+Statistic)
rda_site_tab = rda_site_tab[order(factor(rda_site_tab$Predictor, levels=c('Site',env_vars))),c(1,4,2,3,7,5,6,10,8,9)]

write.csv(rda_site_tab, './Analysis/Figures/Environmental predictors of multi-trait variance decomposition control site.csv', row.names=F )

## Save arrays
save(rda_array, rda_site_array, file='./Analysis/RDA sample mean traits.RData')

# Variance explained based on adjusted R2
rda_site_array[,'Total',c('R2','P')]

# Proportion of total trait variation that can be exaplained by environmentally constrained axis
# based on eigenvalue decomposition
rda_total_site = melt(rda_site_array[,'Total',])
rda_total_site_tab = cast(rda_total_site, Predictor~Statistic)
rda_total_site_tab$Proportion = rda_total_site_tab$Var / rda_total_site_tab$Tot


# Figure including  column for effect of site
pdf('./Analysis/Figures/Multi-trait Eigenval Decomposition control site.pdf', height=6, width=7)
par(mar=c(11, 4, 0, 7))
plot.new()
plot.window(xlim=c(0, (length(env_vars) + 2)*(bar_width+space_width)), ylim=c(-.15, 1))
abline(h=0)
for(i in 1:dim(rda_site_array)[1]){
	varcomp = rda_site_array[i,c('Interspecific','Intraspecific','Covariation'),'Var'] # Explained variance for each component
	varcomp = varcomp / rda_site_array[i,'Total','Tot'] # Scale by total variance
	
	these_y = aov_heights(varcomp)
	rect(i*space_width+(i-1)*bar_width, these_y[1,], i*(space_width+bar_width), these_y[2,],
		col=bar_col)
}

varcomp = rda_site_array[1,c('Interspecific','Intraspecific','Covariation'),'Tot'] / rda_site_array[1,'Total','Tot']
these_y = aov_heights(varcomp)
rect((i+1)*space_width+(i)*bar_width, these_y[1,], (i+1)*(space_width+bar_width), these_y[2,], col=bar_col)

axis(2, las=3)
mtext('Proportion Variance', 2, 2.8)
par(xpd=T)
text((1:(length(env_vars)+1))*(space_width+bar_width)-0.5*bar_width, 
	-.15, c('Site',xvarnames[env_vars]), adj=c(1,0.5), srt=90)
legend('right',colnames(aov0), fill=bar_col, bty='n', inset=-.25)
text((i+1)*(space_width+bar_width)-0.5*bar_width, -.15, 'Total', adj=c(1,0.5), srt=90)
par(xpd=F)
dev.off()




################################################################
### Linear Models of Trait Variation at sample scale (trait diversity)
library(reshape) #melt

# We will use model_data where traits are log-transformed for consistency with mean models

## Calculate metrics of single trait dispersion for each sample
samps = unique(model_data$SampID)
cft_disp = array(NA, dim=c(length(samps),7,2), dimnames=list(SampID=samps, Trait=use_traits, Metric = c('sd','mpd')))
for(i in use_traits){
	cft_disp[,i,'sd'] = tapply(model_data[,i], model_data$SampID, function(x) sqrt(var(x, na.rm=T)))
	cft_disp[,i,'mpd'] = tapply(model_data[,i], model_data$SampID, function(x) mean(dist(x[!is.na(x)])))
}


# Bootstrap null distribution for trait dispersion
N=1000
cft_disp_null = array(NA, dim=c(length(samps),7,2,N), dimnames=list(SampID=samps, Trait=use_traits, Metric= c('sd','mpd'), 1:N))

for(i in use_traits){
	NAinds = which(is.na(model_data[,i]))
	x = model_data[-NAinds, i]
	fact = model_data[-NAinds, 'SampID']
	
	for(j in 1:N){
		use_order = sample(fact)
		sdev = tapply(x, use_order, function(y) sqrt(var(y, na.rm=T)))
		mpd = tapply(x, use_order, function(y) mean(dist(y[!is.na(y)])))
	
		cft_disp_null[names(sdev),i,'sd',j] = sdev
		cft_disp_null[names(mpd),i,'mpd',j] = mpd
	}
}	


# Calculate z-scores for actual trait dispersions based on null distribution
null_mean = apply(cft_disp_null, c(1,2,3), mean)
null_sd = apply(cft_disp_null, c(1,2,3), function(x) sqrt(var(x)))
cft_disp_z = (cft_disp - null_mean)/null_sd


cft_disp_df = melt(cft_disp)
cft_disp_df$z = melt(cft_disp_z)[,'value']

cft_disp_df = merge(cft_disp_df, env[,c('SampID',env_vars)])
cft_disp_df = merge(cft_disp_df, samples[,c('SampID','Year')])

use_col = mycolor[c(2,7)]

pdf('./Analysis/Figures/FT sd vs env.pdf', height=15, width=15)
par(mfrow=c(length(use_traits), length(env_vars)))
par(mar=c(4,4,1,1))
for(i in use_traits){
for(j in env_vars){

	this_data = subset(cft_disp_df, Trait==i&Metric=='sd')
	plot(z~this_data[,j], data=this_data, xlab=j, ylab=i, pch=21, bg=use_col[factor(this_data$Year)], las=1)
	abline(h=c(-2,0,2), lty=2)

	mod = lm(z~this_data[,j], data=this_data)
	sig = coef(summary(mod))[2,4] <0.05

	if(sig){
		modfunc = function(x) coef(mod)%*%rbind(1,x)
		curve(modfunc(x), from=min(this_data[,j]), to=max(this_data[,j]), lwd=2, add=T)
	}
}}
dev.off()










