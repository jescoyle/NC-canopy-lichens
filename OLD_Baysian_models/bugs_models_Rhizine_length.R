## This script fits a set of single predictor models for a given trait using OpenBUGS

library('R2OpenBUGS')
library('coda')
library('lattice')		

# Read in Data
use_data = read.csv('../CFT_model_data.csv')

# Set trait and find missing observations
i = 'Rhizine_length'
keeprows = !is.na(use_data[,i])

# Set model variables
y = use_data[keeprows,i]
n = length(y)
samp = as.numeric(factor(use_data[keeprows,'SampID']))
J = length(unique(samp))
genus = as.numeric(factor(use_data[keeprows,'Genus']))
G = length(unique(genus))
use_xvars = c('Light_mean','Light_high','Temp_max','Vpd_mean','Vpd_daysatfreq')

# Save levels corresponding to genera and samples
write.csv(levels(factor(use_data[keeprows,'Genus'])), 'genus_factor_levels.csv', row.names=T)
write.csv(levels(factor(use_data[keeprows,'SampID'])), 'sample_factor_levels.csv', row.names=T)

# Loop through predictor variables
for(j in use_xvars){

	# Define predictor variable
	xvar = use_data[keeprows,j]

	# Set up model
	mod_data = list('y','samp','genus','n','J','G','xvar')
	mod_init = function(){list(a=rnorm(J), b=rnorm(G), g0=rnorm(1), g1=rnorm(1), 
		sigma.y=runif(1), sigma.a=runif(1), sigma.b=runif(1))}
	mod_parms = c('a','b','g0','g1','sigma.y','sigma.a','sigma.b')
	mod_file = paste('bayes_mod_',i,'.txt', sep='')
	
	# Create a new working directory and copy model file into it
	if(!file.exists(j)){
		dir.create(j)
		file.copy(mod_file, j)
	}	

	# Run model
	modbug = bugs(mod_data, mod_init, mod_parms, mod_file, 
		n.chains=3, n.iter=55000, codaPkg=T, n.burnin=0, n.thin=100,
		saveExec=T, working.directory=j, restart=F)
	
	# Load model results
	modbug = read.bugs(modbug)

	# Save as RData file which can be accessed later
	setwd(j)
	save(modbug, file=paste(i,j,'model.RData', sep='-'))	

	# Thin results
	thinned = window(modbug, thin=1, start=5001)

	## Plot results of each chain
	mylayout=c(2,5,7)

	# Trace
	pdf(paste(i,j,'trace.pdf', sep='-'), height=11, width=8.5)
	print(xyplot(thinned, layout=mylayout, scales=list(rot=0)))
	dev.off()

	# Estimates
	pdf(paste(i,j,'density.pdf', sep='-'), height=11, width=8.5)
	print(densityplot(thinned, layout=mylayout, scales=list(rot=0), aspect='fill'))
	dev.off()

	# Autocorrelation
	pdf(paste(i,j,'ac.pdf', sep='-'), height=11, width=8.5)
	print(acfplot(thinned, layout=mylayout, scales=list(rot=0), aspect='fill'))
	dev.off()

	# Convergence
	sink(paste(i,j,'gelman_convergence.txt', sep='-'))
	try(gelman.diag(thinned, autoburnin=F),silent=T)
	sink()

	print(paste('Done with', j))
	setwd('../')
}

	
