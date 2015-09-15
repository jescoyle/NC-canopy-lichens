## This script was used to compare Maximum Liklihood mixed models to Baysian MCMC estimation for the detection of sources of variance in traits




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


