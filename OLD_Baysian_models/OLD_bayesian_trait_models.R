###########################################
### Fit Baysian Models

samples$Height_top

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

