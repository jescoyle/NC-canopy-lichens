## This script assesses variation in reproductive mode for canopy lichens

options(stringsAsFactors=F)

### Read in Data

working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/'
sql_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Data/SQLite Tables/'

setwd(working_dir)

thalli = read.csv(paste(sql_dir, 'derived_traits.csv', sep=''))
samples = read.csv(paste(sql_dir, 'samples.csv', sep=''))

# Create data frame for analysis
lichens = merge(thalli, samples, all.x=T)


####################################################################
### GLMs

library(lme4)

# Subset lichens to avoid samples on T10-12 where traits were not measured
keep_samps = subset(samples, TreeID<10)$SampID
lichens = subset(lichens, SampID %in% keep_samps)

# Convert reproductive abundance to factor
lichens$Sexual_abun = factor(lichens$Sexual_abun, levels=c('none','few','several','many'))
lichens$Asexual_abun = factor(lichens$Asexual_abun, levels=c('none','few','several','many'))
lichens$Sexual_abun = as.numeric(lichens$Sexual_abun)-1
lichens$Asexual_abun = as.numeric(lichens$Asexual_abun)-1

# Drop 2 lichens where reproduction was not quantified
sum(is.na(lichens$Asexual_abun))
sum(is.na(lichens$Sexual_abun))
lichens = subset(lichens, !is.na(Asexual_abun)&!is.na(Sexual_abun))

# Define variable for asexual vs sexual presence
lichens$Asco = lichens$Sexual_abun!=0
lichens$Asex = lichens$Asexual_abun!=0

# Define a variable for reproductive mode
lichens$RepMode = factor(colSums(t(lichens[,c('Asco','Asex')])*c(1,10)), levels=c(0,1,10,11))
levels(lichens$RepMode) = c('none','sexual','asexual','both')

# Calculate prevalence of each mode within samples
rep_pres = xtabs(~SampID+RepMode, data=lichens)
tot_pres = rowSums(rep_pres)

# Calculate response variables
asco_pres = rowSums(rep_pres[,c('sexual','both')])
asex_pres = rowSums(rep_pres[,c('asexual','both')])
norep_pres = rep_pres[,'none']

# Check for missing data - Sample 30 appears to not have any lichens
missing_lich = subset(samples, !(SampID %in% names(asco_pres)))$SampID
subset(lichens, SampID %in% missing_lich)

# Make a dataframe of responses
reproduction = data.frame(SampID = rownames(rep_pres), Asco_pres=asco_pres, 
	Asex_pres=asex_pres, Norep_pres=norep_pres, N=tot_pres)

# Save data
#write.csv(reproduction, './Data/Derived Tables/reproduction.csv', row.names=F)
reproduction = read.csv('./Data/Derived Tables/reproduction.csv')
rownames(reproduction) = reproduction$SampID


### Logistic regression on reproductive mode presence

reproduction = merge(reproduction, samples, all.x=T)

# Define data set and covariates
use_data = reproduction

x = cbind(use_data$Asco_pres, use_data$N)
asco_mod = glm(x~Diam, data=use_data, family=binomial(link='logit'))
asco_mod = glmer(x~Diam+(1|TreeID), data=use_data, family=binomial(link='logit'))
asco_mod = glm(x~Height, data=use_data, family=binomial(link='logit'))

x = cbind(use_data$Asex_pres, use_data$N)
asex_mod = glm(x~Diam, data=use_data, family=binomial(link='logit'))
asex_mod = glmer(x~Diam+(1|TreeID), data=use_data, family=binomial(link='logit'))
asex_mod = glm(x~Height, data=use_data, family=binomial(link='logit'))



### Regression on average reprocutive mode effort

# Calculate average reproductive effort
asex_eff = aggregate(lichens$Asexual_abun, by=list(SampID=lichens$SampID), FUN=mean)
names(asex_eff)[2] = 'Asex_eff'
asco_eff = aggregate(lichens$Sexual_abun, by=list(SampID=lichens$SampID), FUN=mean)
names(asco_eff)[2] = 'Asco_eff'
rep_eff = merge(asex_eff,  asco_eff)

reproduction = merge(reproduction, rep_eff)

plot(Asex_eff~Asco_eff, data=reproduction)

plot(Asex_eff~Diam, data=reproduction)
plot(Asco_eff~Diam, data=reproduction)


### Hierarchical model of individual reproductive effort



