## This script summarizes and examines environmental data from data loggers

working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/'
sql_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Data/SQLite Tables/'

setwd(working_dir)

library(PerformanceAnalytics) # chart.Correlation
library(sp) # spDistsN1

mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = rgb(mycol, maxColorValue=255)[10:1]

##################################################################################
### Read in data
samples = read.csv(paste(sql_dir, 'samples.csv', sep=''))
loggers = read.csv(paste(sql_dir, 'loggers.csv', sep=''))


##############################################################################
### Data Stats

## Convert textual timepoints to actual dates and times
# Put date and time columns together into string
timepoints = paste(loggers$Date, loggers$Time)

# Convert times to UTC using timezone during which events were recorded
datetimes_utc = as.POSIXlt(timepoints, tz='UTC')
datetimes_utc = datetimes_utc + loggers$Time_zone*(60*60)

# Convert UTC times back to U.S. Eastern without DST
datetimes = as.POSIXlt(datetimes_utc, tz='Etc/GMT+5')
times = datetimes$hour + datetimes$min/60
loggers$datetimes = datetimes
loggers$times = times

## Convert branch positions to factors
samples$BranchPos = factor(samples$BranchPos, levels=c('low','mid','high'))

## Examine all loggers time series to identify outliers
temprange = range(loggers$Temp)
lightrange = range(loggers$Lux)

plotorder = samples[order(samples$TreeID, samples$Pair, samples$BranchPos),'SampID']

pdf('./Analysis/Figures/logger temp time series.pdf', height=8.5, width=11)
par(mfrow=c(3,1))
par(mar=c(4,6,3,1))
for(i in plotorder){
	use_loggers = subset(loggers, SampID==i)
	this_samp = subset(samples, SampID==i)
	plotlabel = paste('T',this_samp$TreeID,'-S',i,', ',this_samp$Pair, ' ',this_samp$BranchPos, ' branch', sep='')
	plot(use_loggers$datetimes, use_loggers$Temp, type='l', ylim=temprange, 
		main=plotlabel, xlab='', ylab='Temperature', las=1)
}
dev.off()

pdf('./Analysis/Figures/logger light time series.pdf', height=8.5, width=11)
par(mfrow=c(3,1))
par(mar=c(4,6,3,1))
for(i in plotorder){
	use_loggers = subset(loggers, SampID==i)
	this_samp = subset(samples, SampID==i)
	plotlabel = paste('T',this_samp$TreeID,'-S',i,', ',this_samp$Pair, ' ',this_samp$BranchPos, ' branch', sep='')
	plot(use_loggers$datetimes, use_loggers$Lux, type='l', ylim=lightrange, 
		main=plotlabel, xlab='', ylab='', las=1)
	mtext('Light Intensity',2,4)
	
}
dev.off()

# What is the distribution of light and temperature
hist(loggers$Lux)
quantile(loggers$Lux, probs=c(0,0.01,0.05,0.1,0.25,0.5,.75,.9,.95,.99,1), na.rm=T)
hist(loggers$Temp)

### Calculate environmental data summaries for each sample

## Calculate humidity from measured temp and weather station humidity
weather_df = read.table('./Data/Weather_stations/NOAA_NCDC_NC-Durham-11W_subhourly.txt', header=T)
weather_wdf = read.table('./Data/Weather_stations/Person_County_Airport_NCDC_Surface_Data_Hourly_Global_edited.txt',
	header=F, sep=',', skip=2, na.strings=c('','999','9999.9','999.9'))
weather_kor = read.csv('./Data/Weather_stations/Korstian/Korstian-High-TRH.csv')


wdf_header1 = read.table('./Data/Weather_stations/Person_County_Airport_NCDC_Surface_Data_Hourly_Global_edited.txt', nrows=1)
wdf_header2 = read.table('./Data/Weather_stations/Person_County_Airport_NCDC_Surface_Data_Hourly_Global_edited.txt', nrows=1, skip=1)
wdf_header2[10:11] = c('Dir_Q', 'Dir_I')
wdf_header2[13] = 'Spd_Q'
#wdf_header2[15] = 'Dewpt_Q'
#wdf_header2[17] = 'Slp_Q'
#wdf_header2[20:21] = c('Pr_I','Pr_Q')

wdf_header = unlist(c(wdf_header2[1,1:8], paste(wdf_header1[1,2], wdf_header2[1,9:13], sep='_'), 
	paste(wdf_header1[1,3], wdf_header2[1,14:15], sep='_'),
	paste(wdf_header1[1,4], wdf_header2[1,16:17], sep='_'),
	paste(wdf_header1[1,5], wdf_header2[1,18:21], sep='_')
))
names(weather_wdf) = wdf_header

## Convert weather station times to actual times

# Add initial 0 to time column when time before 10am
weather_df$UTC_TIME = gsub(' ', 0, format(weather_df$UTC_TIME, width=4))
weather_wdf$HrMn = gsub(' ', 0, format(weather_wdf$HrMn, width=4))
weather_kor$datetime_utc = as.POSIXlt(weather_kor$Time, tz='UTC', format='%m/%d/%Y %H:%M')
weather_kor$datetime_utc = weather_kor$datetime_utc + weather_kor$Time_zone*(60*60)
weather_df$datetime = strptime(paste(weather_df$UTC_DATE, weather_df$UTC_TIME), format='%Y%m%d %H%M', tz='UTC')
weather_wdf$datetime = strptime(paste(weather_wdf$Date, weather_wdf$HrMn), format='%Y%m%d %H%M', tz='UTC')

# Add a unique identifier for each row
weather_df$weather_row = 1:nrow(weather_df)
weather_wdf$weather_row = 1:nrow(weather_wdf)
weather_kor$weather_row = 1:nrow(weather_kor)

# A function that finds the closest time
nearest_time = function(x, timelist){
	timediffs = abs(x-timelist)
	best_times = which(timediffs==min(timediffs))
	if(length(best_times)>1) best_times = best_times[1]
	best_times	
}

# Merge 2013-2014 records with Person County Airport data
log_wdf = subset(loggers, SampID <19)
weather_row = c()
for(i in 1:nrow(log_wdf)){
	weather_row = c(weather_row, nearest_time(log_wdf$datetimes[i], weather_wdf$datetime))
	if(i %% 100 == 0) print(i)
}
log_wdf$weather_row = weather_row

# Define columns in weather station data to merge into logger data
keep_cols = c('WIND_Spd','WIND_Spd_Q','DEWPT_Dewpt','DEWPT_Q','SLP_Slp','SLP_Q','PRECIP_Pr','PRECIP_Amt','PRECIP_I','PRECIP_Q')

# Merge and save
log_wdf = merge(log_wdf, weather_wdf[,c('weather_row',keep_cols)], by='weather_row', all.x=T)
write.csv(log_wdf, './Data/Derived Tables/logger&weather_data_wdf.csv', row.names=F)
log_wdf = read.csv('./Data/Derived Tables/logger&weather_data_wdf.csv')

# Merge 2014-2015 records with Duke Forest data
log_df = subset(loggers, SampID > 18)
weather_row = c()
for(i in 1:nrow(log_df)){
	weather_row = c(weather_row, nearest_time(log_df$datetimes[i], weather_df$datetime))
	if(i %% 100 == 0) print(i)
}
log_df$weather_row = weather_row

# Define columns in weather station data to merge into logger data
keep_cols = c('AIR_TEMPERATURE','PRECIPITATION','SOLAR_RADIATION','SR_FLAG','RELATIVE_HUMIDITY','RH_FLAG')

# Merge and save
log_df = merge(log_df, weather_df[,c('weather_row',keep_cols)], by='weather_row', all.x=T)
write.csv(log_df, './Data/Derived Tables/logger&weather_data_df.csv', row.names=F)
log_df = read.csv('./Data/Derived Tables/logger&weather_data_df.csv')

# Merge 2014-2015 records with Brenna's Korstian Division HOBO logger data
log_kor = subset(loggers, SampID > 18)
weather_row = c()
for(i in 1:nrow(log_kor)){
	weather_row = c(weather_row, nearest_time(log_kor$datetimes[i], weather_kor$datetime_utc))
	if(i %% 100 == 0) print(i)
}
log_kor$weather_row = weather_row

# Define columns in weather station data to merge into logger data
keep_cols = c('RH','Temp_C')

# Merge and save
log_kor = merge(log_kor, weather_kor[,c('weather_row',keep_cols)], by='weather_row', all.x=T)
write.csv(log_kor, './Data/Derived Tables/logger&weather_data_kor.csv', row.names=F)
log_kor = read.csv('./Data/Derived Tables/logger&weather_data_kor.csv')
# This doesn't actually work because the last date sampled for the Korstian Dataloggers is 2/14/2015

## Calculate VPD based on weather station humidity and data logger temp

# A function to calculate vapor pressure deficit from temp and rh or dewpt
# We do not use the slight correction for pressure
# based on: Buck 1981, New Equations for Computing Vapor Pressure and Enhancement Factor
# http://journals.ametsoc.org/doi/abs/10.1175/1520-0450%281981%29020%3C1527%3ANEFCVP%3E2.0.CO%3B2
# these constants are recommended for -20 : 50 Celcius
calc_vpd = function(temp, dewpt=NA, rh=NA, rh_temp=NA){
	a = 6.1121
	b = 17.502
	c = 240.97

	# Calculate the vapor pressure of saturated air in the canopy based on data logger temp
	vp_canopysat = a*exp((b*temp)/(c+temp))
	
	# If dewpoint temperature is given, calculate the actual vapor pressure in the air at the weather station
	if(is.na(rh)[1]){
		vp_act = a*exp((b*dewpt)/(c+dewpt))
	}

	# If relative humidity is given, calculate
	if(is.na(dewpt)[1]){
		vp_sat = a*exp((b*rh_temp)/(c+rh_temp))
		vp_act = vp_sat*(rh/100)
	}

	# Vapor pressure deficit is the difference between the saturated vapor pressur in the canopy
	# and the actual vapor pressure at the weather station
	vpd = vp_canopysat - vp_act
	
	vpd
}

# Assign NA to erroneous data indicated by flags
table(log_wdf$DEWPT_Q) # all codes ok
table(log_df$RH_FLAG) # 3 indicates erroneous data
log_df[log_df$RH_FLAG==3, 'RELATIVE_HUMIDITY'] = NA
log_df[log_df$RELATIVE_HUMIDITY==-9999 &!(is.na(log_df$RELATIVE_HUMIDITY)), 'RELATIVE_HUMIDITY'] = NA
log_df[log_df$AIR_TEMPERATURE==-9999,'AIR_TEMPERATURE'] = NA
log_wdf[log_wdf$DEWPT_Dewpt==999.9,'DEWPT_Dewpt'] = NA

# Calculate
log_df$vpd = calc_vpd(temp=log_df$Temp, rh=log_df$RELATIVE_HUMIDITY, rh_temp=log_df$AIR_TEMPERATURE)
log_wdf$vpd = calc_vpd(temp=log_wdf$Temp, dewpt=log_wdf$DEWPT_Dewpt)
log_kor$vpd = calc_vpd(temp=log_kor$Temp, rh=log_kor$RH, rh_temp = log_kor$Temp_C)

# Save data
write.csv(log_df, './Data/Derived Tables/logger&weather_data_df.csv', row.names=F)
write.csv(log_wdf, './Data/Derived Tables/logger&weather_data_wdf.csv', row.names=F)
write.csv(log_kor, './Data/Derived Tables/logger&weather_data_kor.csv', row.names=F)


# Compre vpd calculated using Duke Forest weather station vs. Korstian Data Logger
names(log_kor)[13] = 'vpd_kor'
names(log_df)[17] = 'vpd_df'
log_kor$ID = paste(log_kor$SampID, log_kor$datetimes)
log_df$ID = paste(log_df$SampID, log_df$datetimes)
comp_wd = merge(log_kor[,c('SampID','datetimes','Temp','RH','Temp_C','vpd_kor', 'ID')], log_df[,c('ID','AIR_TEMPERATURE','RELATIVE_HUMIDITY','vpd_df')])

cutoff = which((comp_wd$datetimes > weather_kor[nrow(weather_kor),'datetime_utc']))
comp_wd = comp_wd[-cutoff,]

pdf('./Data/Weather_stations/compare_vpd_korstian-duke_forest.pdf', height=11, width=8)
par(mfrow=c(6,3))
par(mar=c(4,4,1,1))
for(i in 19:72){
	plot(vpd_df~vpd_kor, data=subset(comp_wd, SampID==i))
	abline(0,1, col=2, lwd=2)
	mtext(paste('Sample', i), 3, -1)
}
dev.off()

plot(AIR_TEMPERATURE~Temp_C, data=comp_wd, xlab='Korstian Temp', ylab='Duke Forest Temp')
abline(0,1, col=2, lwd=2)
plot(RELATIVE_HUMIDITY~RH, data=comp_wd, xlab='Korstian RH', ylab='Duke Forest RH')
abline(0,1, col=2, lwd=2)



# Make a dataframe with just environmental variables to be analyzed
keep_cols = c('Index','SampID','Temp','Lux','vpd','Time_zone','Date','Time','datetimes','times')
loggers_env = rbind(log_wdf[,keep_cols], log_df[,keep_cols])
write.csv(loggers_env, './Data/Derived Tables/loggers_env.csv', row.names=F)

# DID NOT QUITE PLOT RIGHT
# SOME TIMES OUT OF ORDER
# BRANCH HEIGHT OUT OF ORDER
pdf('./Analysis/Figures/logger vpd time series.pdf', height=8.5, width=11)
par(mfrow=c(3,1))
par(mar=c(4,6,3,1))
for(i in plotorder){
	use_loggers = subset(loggers_env, SampID==i)
	this_samp = subset(samples, SampID==i)
	plotlabel = paste('T',this_samp$TreeID,'-S',i,', ',this_samp$Pair, ' ',this_samp$BranchPos, ' branch', sep='')
	plot(use_loggers$datetimes, use_loggers$vpd, type='l', ylim=c(-5,100), 
		main=plotlabel, xlab='', ylab='', las=1)
	mtext('Vapor Pressure Deficit',2,4)
	
}
dev.off()

###########################################################
### Calculate summary data for samples
options(stringsAsFactors=F)
loggers_env = read.csv('./Data/Derived Tables/loggers_env.csv')
loggers_env$datetimes = as.POSIXlt(loggers_env$datetimes, tz='Etc/GMT+5')

## Make array of logger data - puts in 0 for missing data

# Factor dates
dates = unique(loggers_env$datetimes)
dates = dates[order(dates)]
loggers_env$datefactor = factor(loggers_env$datetimes, levels=dates)

temp_mat = as.matrix(xtabs(Temp~datefactor+SampID, data=loggers_env))
lux_mat = as.matrix(xtabs(Lux~datefactor+SampID, data=loggers_env))
vpd_mat = as.matrix(xtabs(vpd~datefactor+SampID, data=loggers_env))

# Note that absent datetime levels get a 0 entry.
# Use temp_mat to identify datefactors not recorded for each sample and set to NA.
# Works because loggers never logged temperature exactly = 0.
temp_0 = temp_mat==0
temp_mat[temp_0] = NA
lux_mat[temp_0] = NA
vpd_mat[temp_0] = NA

# Set date window for which data will be analyzed
# Change to NA any values recorded outside of these times

start2013_date = as.POSIXlt('2013-06-19 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT+5')
end2013_date = as.POSIXlt('2014-03-29 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT+5')
start2014_date = as.POSIXlt('2014-06-19 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT+5')
end2014_date = as.POSIXlt('2015-03-29 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT+5')
keep2013 = (dates>=start2013_date)&(dates<=end2013_date)
keep2014 = (dates>=start2014_date)&(dates<=end2014_date)

temp_mat[!keep2013,1:18]<-NA
lux_mat[!keep2013,1:18]<-NA
vpd_mat[!keep2013,1:18]<-NA
temp_mat[!keep2014,19:72]<-NA
lux_mat[!keep2014,19:72]<-NA
vpd_mat[!keep2014,19:72]<-NA

image(temp_mat)
image(lux_mat)
image(vpd_mat)


## Calculate temperature variables June - Oct and Nov - Mar
times = dates$hour + dates$min/60
summer = (dates$mon<=10) & (dates$mon>=5)
winter = (dates$mon>=11) | (dates$mon<=3) # note that Dec is 0
plot(dates$mon, col=c(1,2)[summer+1])
plot(dates$mon, col=c(1,2)[winter+1])


Temp_mean = apply(temp_mat, 2, mean, na.rm=T)
Temp_max = apply(temp_mat, 2, max, na.rm=T) 
Temp_min = apply(temp_mat, 2, min, na.rm=T)
Temp_sd = apply(temp_mat, 2, sd, na.rm=T)
Temp_mean_sum = apply(temp_mat[summer,], 2, mean, na.rm=T)
Temp_max_sum = apply(temp_mat[summer,], 2, max, na.rm=T) 
Temp_min_sum = apply(temp_mat[summer,], 2, min, na.rm=T)
Temp_sd_sum = apply(temp_mat[summer,], 2, sd, na.rm=T)
Temp_mean_win = apply(temp_mat[winter,], 2, mean, na.rm=T)
Temp_max_win = apply(temp_mat[winter,], 2, max, na.rm=T)
Temp_min_win = apply(temp_mat[winter,], 2, min, na.rm=T)
Temp_sd_win = apply(temp_mat[winter,], 2, sd, na.rm=T)

Temp = data.frame(Temp_mean, Temp_max, Temp_min, Temp_sd,
	Temp_mean_sum, Temp_max_sum, Temp_min_sum, Temp_sd_sum, 
	Temp_mean_win, Temp_max_win, Temp_min_win, Temp_sd_win)


# Relationships among temperature variables
chart.Correlation(Temp)


## Calculate daytime light variables June - Oct and Nov - Mar

# Find night times
nightTime=rowSums(lux_mat==0, na.rm=T)>0
table(rowSums(lux_mat==0, na.rm=T)[nightTime])
hist(rowSums(lux_mat==0, na.rm=T)[nightTime])
# Peaks at 18 and 54 indicate that most loggers read 0 light during the same intervals

# When was 0 light recorded between 7am and 5pm?
day_times = times>7 & times < 17 
lux_mat[day_times,][rowSums(lux_mat[day_times,]==0, na.rm=T)>0,] # Dec and Jan

# Omit all instances with 0 light from calculations- assume these occur at night
dayTime = lux_mat != 0

## Calculate light variables for June-Oct and Nov-Mar

# Replace 0 light with NA
lux_mat[lux_mat==0]<- NA

Light_mean = apply(lux_mat, 2, mean, na.rm=T)
Light_max = apply(lux_mat, 2, max, na.rm=T)
Light_sd = apply(lux_mat, 2, sd, na.rm=T)
Light_mean_sum = apply(lux_mat[summer,], 2, mean, na.rm=T)
Light_max_sum = apply(lux_mat[summer,], 2, max, na.rm=T)
Light_sd_sum = apply(lux_mat[summer,], 2, sd, na.rm=T)
Light_mean_win = apply(lux_mat[winter,], 2, mean, na.rm=T)
Light_max_win = apply(lux_mat[winter,], 2, max, na.rm=T)
Light_sd_win = apply(lux_mat[winter,], 2, sd, na.rm=T)

# How frequently are high light events recorded?
highLight = lux_mat>25000
nlightobs = apply(highLight, 2, function(x) sum(!is.na(x)))
nlightobs_sum = apply(highLight[summer,], 2, function(x) sum(!is.na(x)))
nlightobs_win = apply(highLight[winter,], 2, function(x) sum(!is.na(x)))

Light_p90 = apply(lux_mat, 2, function(x) quantile(x, .9, na.rm=T))
Light_high = apply(highLight, 2, sum, na.rm=T)/nlightobs
Light_p90_sum = apply(lux_mat[summer,], 2, function(x) quantile(x, .9, na.rm=T))
Light_high_sum = apply(highLight[summer,], 2, sum, na.rm=T)/nlightobs_sum
Light_p90_win = apply(lux_mat[winter,], 2, function(x) quantile(x, .9, na.rm=T))
Light_high_win = apply(highLight[winter,], 2, sum, na.rm=T)/nlightobs_win

# When were high light levels observed? This gets at duration of time a sample can receive high light
Light_window = sapply(1:72, function(x){

	use_times = times[which(highLight[,x])]

	if(length(use_times)>0){
		range(use_times, na.rm=T)
	} else {
		c(NA,NA)
	}	
})
Light_window_sum = sapply(1:72, function(x){

	use_times = times[which(highLight[summer,x])]

	if(length(use_times)>0){
		range(use_times, na.rm=T)
	} else {
		c(NA,NA)
	}	
})
Light_window_win = sapply(1:72, function(x){

	use_times = times[which(highLight[winter,x])]

	if(length(use_times)>0){
		range(use_times, na.rm=T)
	} else {
		c(NA,NA)
	}	
})

Light_range = apply(Light_window, 2, diff)
Light_range_sum = apply(Light_window_sum, 2, diff)
Light_range_win = apply(Light_window_win, 2, diff)

Light_range[is.na(Light_range)] = 0 
Light_range_sum[is.na(Light_range_sum)] = 0 
Light_range_win[is.na(Light_range_win)] = 0 


# Combine light variables into a data frame
Light = data.frame(Light_mean, Light_max, Light_sd, Light_high, Light_p90, Light_range,
	Light_mean_sum, Light_max_sum, Light_sd_sum, Light_high_sum, Light_p90_sum, Light_range_sum,
	Light_mean_win, Light_max_win, Light_sd_win, Light_high_win, Light_p90_win, Light_range_win)

# Relationships among temperature variables
chart.Correlation(Light)

# Figure: distribution of light across different times of the day
#mycolor = colorRampPalette(mycol[10:1])(100)[cut(dates, 100)]
mycolor = 'black'
plotorder = samples[order(samples$TreeID, samples$Pair, samples$BranchPos,  decreasing=T),'SampID']

# randomize data order
random_order = sample(1:nrow(lux_mat), nrow(lux_mat), replace=F)

pdf('./Analysis/Figures/Diurnal light distribution across loggers.pdf', height=8, width=10.5)
layout(matrix(1:18, nrow=3))

for(i in plotorder){

	this_samp = subset(samples, SampID==i)
	plotlabel = paste('T',this_samp$TreeID,'-S',i,', ',this_samp$Pair, ' ',this_samp$BranchPos, ' branch', sep='')

	par(mar=c(4,3,1,1))
	plot(lux_mat[random_order,i]~times[random_order], xlab='hour', ylab='lumens', ylim=c(0,254000), xlim=c(0,23.5), col=mycolor)
	text(0,254000,plotlabel, adj=c(0,1))
	abline(h=25000, col=2)
	abline(v=range(times[highLight[,i]], na.rm=T), col='blue')
}
dev.off()


## Calculate VPD variables
Vpd_mean = apply(vpd_mat, 2, mean, na.rm=T)
Vpd_max = apply(vpd_mat, 2, max, na.rm=T) 
Vpd_min = apply(vpd_mat, 2, min, na.rm=T)
Vpd_sd = apply(vpd_mat, 2, sd, na.rm=T)
Vpd_mean_sum = apply(vpd_mat[summer,], 2, mean, na.rm=T)
Vpd_max_sum = apply(vpd_mat[summer,], 2, max, na.rm=T) 
Vpd_min_sum = apply(vpd_mat[summer,], 2, min, na.rm=T)
Vpd_sd_sum = apply(vpd_mat[summer,], 2, sd, na.rm=T)
Vpd_mean_win = apply(vpd_mat[winter,], 2, mean, na.rm=T)
Vpd_max_win = apply(vpd_mat[winter,], 2, max, na.rm=T)
Vpd_min_win = apply(vpd_mat[winter,], 2, min, na.rm=T)
Vpd_sd_win = apply(vpd_mat[winter,], 2, sd, na.rm=T)

nvpdobs = apply(vpd_mat, 2, function(x) sum(!is.na(x)))
nvpdobs_sum = apply(vpd_mat[summer,], 2, function(x) sum(!is.na(x)))
nvpdobs_win = apply(vpd_mat[winter,], 2, function(x) sum(!is.na(x)))

Vpd_satfreq = apply(vpd_mat, 2, function(x) sum(x<=0, na.rm=T))/nvpdobs
Vpd_satfreq_sum = apply(vpd_mat[summer,], 2, function(x) sum(x<=0, na.rm=T))/nvpdobs_sum
Vpd_satfreq_win = apply(vpd_mat[winter,], 2, function(x) sum(x<=0, na.rm=T))/nvpdobs_win

vpd_mat_day = vpd_mat
vpd_mat_day[!dayTime] = NA
nvpdobs_day = apply(vpd_mat_day, 2, function(x) sum(!is.na(x)))
nvpdobs_day_sum = apply(vpd_mat_day[summer,], 2, function(x) sum(!is.na(x)))
nvpdobs_day_win = apply(vpd_mat_day[winter,], 2, function(x) sum(!is.na(x)))

Vpd_daysatfreq = apply(vpd_mat_day, 2, function(x) sum(x<=0, na.rm=T))/nvpdobs_day_sum
Vpd_daysatfreq_sum = apply(vpd_mat_day[summer,], 2, function(x) sum(x<=0, na.rm=T))/nvpdobs_day_sum
Vpd_daysatfreq_win = apply(vpd_mat_day[winter,], 2, function(x) sum(x<=0, na.rm=T))/nvpdobs_day_win

Vpd = data.frame(Vpd_mean, Vpd_max, Vpd_min, Vpd_sd, Vpd_satfreq, Vpd_daysatfreq,
	Vpd_mean_sum, Vpd_max_sum, Vpd_min_sum, Vpd_sd_sum, Vpd_satfreq_sum, Vpd_daysatfreq_sum,
	Vpd_mean_win, Vpd_max_win, Vpd_min_win, Vpd_sd_win, Vpd_satfreq_win, Vpd_daysatfreq_win)

chart.Correlation(Vpd)


# Make a data frame
loggerdata = data.frame(SampID=1:72, Temp, Light, Vpd, nlightobs, nvpdobs, nlightobs_sum, nlightobs_win, nvpdobs_sum, nvpdobs_win)
write.csv(loggerdata, './Data/Derived Tables/loggerdata.csv', row.names=F)


############################################################################
### Explore spatial variation in env


## GOAL: Can I predict environment using spatial information?

loggerdata = read.csv('./Data/Derived Tables/loggerdata.csv', row.names=1)

samples = merge(samples, loggerdata) 


### Plot vertical env profiles by sample position

redcols = colorRampPalette(c('#FF8080', '#FF0000'))(12)
bluecols = colorRampPalette(c('#8080FF', '#0000FF'))(12)

env = 'Vpd_mean'
xlab = 'Mean VPD (mb))'
xlim = c(0,10)

samples$BranchPos = factor(samples$BranchPos, levels=c('low','mid','high'))
samples = samples[order(samples$BranchPos),]


pdf(paste('./Analysis/Figures/Height vs.', env, 'profiles.pdf'), height=3.5, width=4.1)
par(mar=c(5,5,.5,.5))
plot(Height~samples[,env], data=samples, type='n', xlab=xlab, 
	ylab='Sample height (m)', las=1, ylim=c(0,25), xlim=xlim)
for(i in 1:12){
topsamps = subset(samples, TreeID==i & Pair=='top')
points(topsamps$Height~topsamps[,env], type='l', col=redcols[i], lwd=2)
sidesamps = subset(samples, TreeID==i & Pair=='side')
points(sidesamps$Height~sidesamps[,env], type='l', col=bluecols[i], lwd=2)

}
legend('bottomright',c('top','side'), col=c('red','blue'), 
	lty=1, lwd=2, bty='n', title='Branch\nPosition')
dev.off()

# By branch angle
colorby = cut(samples$Angle, breaks=seq(0,90, 10), include.lowest=T)

plot(Height~samples[,env], data=samples, pch=16, col=mycol[10:1][colorby])




########### OLD #################

### Map light and temp in the tree
library(sp) # spDistsN1
library(geosphere) # bearing
library(rgl) # plot3D
library(scatterplot3d) # scatterplot3d

# Calculate sample positions in catesian coordinates
posXYZ = sapply(rownames(samples), function(i){

	# Sample positions are given in cylindical coordinates (Height, PosR, PosTheta)
	X = samples[i,'PosR']*sin(samples[i,'PosTheta']) # note this is sin b/c theta is a bearing
	Y = samples[i,'PosR']*cos(samples[i,'PosTheta']) # note this is cos b/c theta is a bearing
	Z = samples[i,'Height']
	
	c(X,Y,Z)
})
rownames(posXYZ) = c('X','Y','Z')
posXYZ = t(posXYZ) # transpose

# Lower samples on branch sides by 10cm so there are no overlapping samples
posXYZ[which(samples$Pair=='side'),'Z'] = posXYZ[which(samples$Pair=='side'),'Z'] -.1

# Alter x and y coordinates to reflect distances among trees
# Use position of tree 1 as the origin
tree_locs = data.frame(LON=c(-79.026831,-79.026681,-79.027053),LAT=c(36.267149, 36.266045, 36.265785))
tree_locs = as.matrix(tree_locs)
tree_dists = spDistsN1(tree_locs, tree_locs[1,], longlat=T)*1000
tree_theta = bearing(tree_locs[1,],tree_locs)/360*2*pi
deltaY = cos(tree_theta)*tree_dists; deltaY[1] = 0
deltaX = sin(tree_theta)*tree_dists; deltaX[1] = 0
posXYZ[,'Y'] = posXYZ[,'Y'] + deltaY[samples$TreeID]
posXYZ[,'X'] = posXYZ[,'X'] + deltaX[samples$TreeID]

# Plot sample locations colored by tree
mycol = c('darkblue','orange','purple')
plot3d(posXYZ[,'X'], posXYZ[,'Y'], posXYZ[,'Z'], xlab='x', ylab='y', zlab='z',
	size=10, col=mycol[samples$TreeID])

# Colored by environmental variables
mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = apply(mycol,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=256))
mycol = mycol[10:1]

myenv = samples$Light_max_sum
myfact = cut(myenv, breaks=10)

plot3d(posXYZ[,'X'], posXYZ[,'Y'], posXYZ[,'Z'], xlab='', ylab='', zlab='Height (m)',
	size=3, col=mycol[myfact], box=F, type='s')

rgl.postscript("light_max_sum 3D.pdf","pdf")

svg('temp_mean_sum color bar.svg', height=5, width=4)
par(mar=c(0,7,4,0))
color.bar(mycol, min=min(myenv), max=max(myenv), nticks=length(mycol)+1, 
	title='Mean Temp.', lab=expression(paste(degree,'C')), sciNote=F, round=1, 
	cex=2, labLine=4)
dev.off()

## Correlations among variables

library('corrplot')

env_cor = cor(loggerdata[,grep('sum', names(loggerdata))], use='pairwise.complete.obs')
rownames(env_cor) = sapply(rownames(env_cor), function(s) substr(s, 1, nchar(s)-4))
colnames(env_cor) = sapply(colnames(env_cor), function(s) substr(s, 1, nchar(s)-4))


svg('summer environment correlations.svg', height=5, width=8)
corrplot(env_cor, method='ellipse', type='upper', diag=F, 
	order='hclust', hclust.method='complete', tl.cex=1.5, 
	tl.col=1, tl.srt=60, cl.cex=1.5, cl.align.text='l', mar=c(0,0,0,0), col=mycol)
dev.off()
# Have to adjust document dimensions in Inkscape






## Linear models
samples$AspectNS = cos((samples$Aspect/360)*2*pi)
samples$AspectEW = sin((samples$Aspect/360)*2*pi)

temp_mod = lm(Temp_mean~Angle+AspectNS+AspectEW+Height+PosR, data=samples)
summary(temp_mod)

light_mod = lm(Light_mean~Angle+AspectNS+AspectEW+Height+PosR, data=samples)
summary(light_mod)

hilight_mod = lm(Light_high~Angle+AspectNS+AspectEW+Height+PosR, data=samples)
summary(hilight_mod)

# variation explained

sapply(c('Height','Angle','AspectNS','AspectEW','PosR'), function(x){
	sapply(c('Temp_mean','Light_mean','Light_high'), function(y){
		summary(lm(samples[,y]~samples[,x]))$r.squared
})
}) -> r2table

sapply(c('Height','Angle','AspectNS','AspectEW','PosR'), function(x){
	sapply(c('Temp_mean','Light_mean','Light_high'), function(y){
		coef(summary(lm(samples[,y]~samples[,x])))[2,4]
})
}) -> Ptable

write.csv(r2table, 'univariate r2 env vs canopy position.csv')



## Env vs height across pair samples

labs = apply(expand.grid(paste('Tree',1:3), c('side','top')), 1, function(x) paste(x[1],x[2]))

env = 'Temp_mean_sum'

pdf('./Analysis/Figures/Height and position vs mean temp.pdf', height=6, width=15)
par(mar=c(6,6,1,1))
layout(t(matrix(1:3)), widths=c(.4,.4,.2))
plot(samples[,env]~as.numeric(BranchPos), data=samples, type='n', axes=F, xlab='Branch Height', ylab='')
points(samples[,env]~ factor(samples$BranchPos, levels=c('low','mid','high')), data=samples, cex=2, cex.lab=2,
	pch=(0:2)[samples$TreeID], lwd=2, col=c('blue','red')[factor(samples$Pair)])
axis(1, at=1:3, labels=c('low','mid','high'), cex.axis=2)
axis(2, las=2, cex.axis=2)
box()
plot(samples[,env] ~ Height, data=samples, cex=2, xlim=c(0,20), xlab='Height (m)', ylab='',
	pch=(0:2)[samples$TreeID], lwd=2, col=c('blue','red')[factor(samples$Pair)], cex.lab=2, cex.axis=2, las=2)
plot.new()
legend('center',labs, title=env, col=rep(c('blue','red'), each=3), pch=0:2, bty='n', pt.lwd=2,cex=2)
dev.off()


# By vertical angle and aspect
library(plotrix)

myvars = c('Light_mean','Light_max','Light_high','Light_p90','Light_range','Temp_min','Temp_mean','Temp_max')

par(mfrow=c(2,4))
for(i in myvars){
	plot(samples$Angle, samples$Aspect, col=mycol[cut(samples[,i], 10)], cex=2, pch=15,
		xlab='', ylab='', xlim=c(0,90), ylim=c(0,360), main=i)
}

ylabs = c('# Events', 'Lux','°C')
names(ylabs) = c('Light_high','Light_mean','Temp_mean')

png('./angle vs light and temp.png', height=500, width=1500)
par(mfrow=c(1,3))
par(mar=c(4,8,1,1))
for(i in c('Light_high','Light_mean','Temp_mean')){
	plot(samples$Angle, samples[,i], bg=mycol[cut(samples[,i], 10)], cex=4, pch=22,
		ylab='', xlab='', xlim=c(0,90), cex.axis=2.5, las=1)
	mtext(ylabs[i],2,6,cex=2) 
	
	usemod = lm(samples[,i]~samples$Angle)
	abline(usemod, lty=c(2,1)[(summary(usemod)$coefficients[2,4]<0.05 )+1], lwd=2)
	mylabel = bquote(italic(r)^2 == .(format(summary(usemod)$r.squared, digits = 2)))
	legend('topleft',legend=mylabel, bty='n', cex=3)

}

dev.off()



png('./height vs light and temp.png', height=500, width=1500)
par(mfrow=c(1,3))
par(mar=c(4,8,1,1))
for(i in c('Light_high','Light_mean','Temp_mean')){
	plot(samples$Height, samples[,i], bg=mycol[cut(samples[,i], 10)], cex=4, pch=22,
		ylab='', xlab='', xlim=c(0,40), cex.axis=2.5, las=1)
	mtext(ylabs[i],2,6,cex=2) 
	
	usemod = lm(samples[,i]~samples$Height)
	abline(usemod, lty=c(2,1)[(summary(usemod)$coefficients[2,4]<0.05 )+1], lwd=2)
	mylabel = bquote(italic(r)^2 == .(format(summary(usemod)$r.squared, digits = 2)))
	legend('topleft',legend=mylabel, bty='n', cex=3)

}

dev.off()


par(mfrow=c(2,4))
for(i in myvars){
	plot(samples$Aspect, samples[,i], col=mycol[cut(samples[,i], 10)], cex=2, pch=15,
		ylab='', xlim=c(0,360), main=i)
}



# By height and distance from trunk

par(mfrow=c(2,4))
for(i in myvars){
	plot(jitter(samples$PosR, amount=.5), samples$Height, col=mycol[cut(samples[,i], 10)], cex=2, pch=15,
		main=i)
}

par(mfrow=c(2,4))
for(i in myvars){
	plot(samples$Height, samples[,i], col=mycol[cut(samples[,i], 10)], cex=2, pch=15,
		ylab='', main=i)

	usemod = lm(samples[,i]~samples$Height)
	abline(usemod, lty=c(2,1)[(summary(usemod)$coefficients[2,4]<0.05 )+1])
}
par(mfrow=c(2,4))
for(i in myvars){
	plot(samples$Angle, samples[,i], col=mycol[cut(samples[,i], 10)], cex=2, pch=15,
		ylab='', main=i)

	usemod = lm(samples[,i]~samples$Angle)
	abline(usemod, lty=c(2,1)[(summary(usemod)$coefficients[2,4]<0.05 )+1])
}



par(mfrow=c(2,4))
for(i in myvars){
	plot(samples$PosR, samples[,i], col=mycol[cut(samples[,i], 10)], cex=2, pch=15,
		ylab='', main=i)
	
	usemod = lm(samples[,i]~samples$PosR)
	abline(usemod, lty=c(2,1)[(summary(usemod)$coefficients[2,4]<0.05 )+1])
}






## How correlated are light and temp?
## Is there a threshold below which they are not correlated?

plot(Temp_mean~Light_mean, data=loggerdata)
plot(Temp_max~Light_mean, data=loggerdata)
plot(Temp_min~Light_mean, data=loggerdata)
plot(Temp_max~Temp_mean, data=loggerdata)

cor(loggerdata[,2:ncol(loggerdata)])




#############
## Calculate how many high light events we will detect if we sample for n days every x minutes.
## Assume a sample can only receive light for 3.5 hrs = 210 minutes each day (the lowest range of all non-trunk samples)

obs_events = 23
obs_days = 24
obs_events_day = 21

exp_light = function(x,n) (210/x)*(obs_events/obs_days/obs_events_day)*n


exp_light(120,180)


###########################################
#### Make plots comparing two samples for presentation

# Make data matrices
temp_mat = as.matrix(xtabs(Temp~as.character(datetimes)+SampID, data=loggers))
lux_mat = as.matrix(xtabs(Lux~as.character(datetimes)+SampID, data=loggers))

# Make vectors of dates
dates = as.POSIXlt(rownames(temp_mat), format='%Y-%m-%d %H:%M:%S')
times = dates$hour + dates$min/60

# Choose samples
S1 = 1
S2 = 2


## Time series
png(paste('compare timeseries light temp S',S1,' and S',S2,'.png', sep=''), height=300, width=900)
par(mar=c(3,5,1,1))
par(mfrow=c(2,1))
plot(dates, temp_mat[,S1], type='l', ylim=temprange, col='red',
	ylab='', xlab='', las=1)
lines(dates, temp_mat[,S2], col='blue')
mtext('Temperature', 3, -1.1, adj=0.01)
mtext(expression(paste(degree,"C")), 2, 3)
plot(dates, lux_mat[,S1], type='l', ylim=lightrange, col='red',
	ylab='', xlab='', las=1)
lines(dates, lux_mat[,S2], col='blue')
mtext('Light Intensity',3, -1.1, adj=0.01)
mtext('Lux', 2, 4)
dev.off()


	use_loggers = subset(loggers, SampID==i)
	this_samp = subset(samples, SampID==i)

plotlabel = paste('T',this_samp$TreeID,'-S',i,', ',this_samp$Pair, ' ',this_samp$BranchPos, ' branch', sep='')
	plot(use_loggers$datetimes, use_loggers$Temp, type='l', ylim=temprange, 
		main=plotlabel, xlab='', ylab='Temperature', las=1)





###########################################
### Functions


color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), 
	title='', lab='', round=NA, sciNote=F, cex=1, labLine=5){
	scale = length(lut)/(max-min)
	tick_labs = ticks
	if(!is.na(round)){ tick_labs = format(round(tick_labs, round), nsmall=round )}
	if(sciNote==T){ tick_labs = format(tick_labs, scientific=T, digits=round) }
		
	plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
	axis(2, ticks, labels=tick_labs, las=1, cex.axis=cex*.8)
	for (i in 1:length(lut)) {
		y = (i-1)/scale + min
		rect(0,y,10,y+1/scale, col=lut[i], border=NA)
	}
	mtext(lab, 2, labLine, cex=cex)
	mtext(title, 3, 0, cex=cex)

}


######################################################
### OLD CODE

## used to read in logger data

logIDs = c(9729228,9729230,9729231,9729232,9729233,9729234,9729235,9729236,
	9729237,9729240,9729241,9729242,9729243,9729244,10126123, 10126130,10126131,10126134)

logfiles = paste('./May-Nov2013/log',logIDs,'.csv',sep='')

newdata = read.csv(logfiles[1], skip=1, row.names=1) # Examine file structure

# Make an array to hold data from all loggers
readings = array(NA, dim=c(length(logIDs), 3303, 3), 
	dimnames=list(LogID=logIDs,MeasID=1:3303,c('Time','Temp','Light')))

for(i in 1:length(logIDs)){
	# Read in file
	newdata = read.csv(logfiles[i], skip=1, row.names=1)

	# Take only data up to 11/19/2013 at 11:30pm
	newdata = newdata[1:3303,1:3]
	
	# Store in array
	readings[i,,] <- as.matrix(newdata[1:dim(readings)[2],])
}

# Check that times match up across loggers: there should be only one unique row:
unique(readings[,,'Time'])

# Read in again as numeric, now that we know 'Time' will work as a rowname
readings = array(NA, dim=c(length(logIDs), 3303, 2), 
	dimnames=list(LoggerID=logIDs, MeasID=readings[1,,'Time'], c('Temp','Light'))
)

for(i in 1:length(logIDs)){
	newdata = read.csv(logfiles[i], skip=1, row.names=1)
	newdata = newdata[,2:3]
	
	readings[i,,] <- as.matrix(newdata[1:dim(readings)[2],])
}

