### This script reads in HOBO data logger CSV files.
### and appends them to the logger table in the Canopy Lichens SQL database

options(stringsAsFactors=F)

#library(RSQLite)

project_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/'
sql_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/Data/SQLite Tables/'
logger_data_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/DataLoggers/'

setwd(project_dir)

# Set up directory and file names for location of raw data locally
# If there is an existing table logger data, then read it in to be appended to
existing_logger_file = T
this_logger_dir = paste(logger_data_dir, 'Nov2014-Apr2015/', sep='')
csv_samples_filename = 'samples.csv'
csv_logger_filename = 'loggers.csv'

# Set up which sample year this data is from and whether it is the initial or second deployment
year = 2014
deploy = 2

# Find loggerfiles in the logger directory
logger_filenames = list.files(this_logger_dir, pattern='*.csv', full.name=T)

# Get logger dataframe from file (if it exists)
if( existing_logger_file){
	loggers = read.csv(paste(sql_dir, csv_logger_filename, sep=''))
} else { loggers = data.frame() } 

# Get samples data frame
samples  = read.csv(paste(sql_dir, csv_samples_filename, sep=''))

# Find start number for unique index
start_index = ifelse(existing_logger_file, max(loggers$Index)+1, 1)

for( f in logger_filenames){
	this_data = read.csv(f, skip=1, row.names=1)

	# Which data logger is this?
	f_grep = regexpr('log[0-9]+', f)
	this_logID = substr(f, f_grep+3, f_grep+attr(f_grep, 'match.length')-1)

	# Which sample does this correspond to?
	if(deploy==1) this_sampID = subset(samples, LoggerID1==this_logID&Year==year)$SampID
	if(deploy==2) this_sampID = subset(samples, LoggerID2==this_logID&Year==year)$SampID

	# Extract data to be added to sql_loggers file
	add_data = data.frame(SampID = rep(this_sampID, nrow(this_data)), this_data[,1:3])
	colnames(add_data) = c('SampID','Date_time','Temp','Lux')

	# Remove rows that do not correspond to environmental measurements (e.g. logging events)
	if(sum(is.na(add_data$Temp))>0) add_data = add_data[-which(is.na(add_data$Temp)),]

	# Put date/time into different format: YYYY-MM-DD HH:MM:SS
	form1 = '%m/%d/%y %I:%M:%S %p'
	form2 = '%m/%d/%Y %H:%M'


	if(!is.na(strptime(add_data$Date_time[1], format=form1))){
		use_form=form1
	} else {
		if(!is.na(strptime(add_data$Date_time[1], format=form2))) use_form=form2
	}

	date_time = strptime(add_data$Date_time, format=use_form)
	add_data$Date = as.character(format(date_time, format='%Y-%m-%d'))
	add_data$Time = as.character(format(date_time, format='%H:%M:%S'))

	# Add unique index
	add_data$Index = start_index:(start_index+nrow(add_data)-1)
	start_index = start_index + nrow(add_data)

	# Append this data to rfull data
	loggers = rbind(loggers, add_data[,c('Index','SampID','Date','Time','Temp','Lux')])
}


# Re-index
loggers$Index = 1:nrow(loggers)

## read back in previously saved full logger data
loggers = read.csv(paste(logger_data_dir, 'all_logger_data.csv', sep=''))


### Remove observations during times when data loggers were not deployed at sample locations

### Note that DST is GMT+4 and EST is GMT+5, contrary to normal convention
good_loggers = data.frame()
for(i in 1:72){

	this_samp = subset(samples, SampID==i)
	these_loggers = subset(loggers, SampID==i)
	
	# Calculate time zone for when logger data was dumped in Novemeber
	# In 2014, some loggers were read out before DST ended and some were read out after DST ended
	wintz = ifelse(this_samp$Year==2014 & as.Date(this_samp$Logger_dump_date) < as.Date('2014-11-02'), 'Etc/GMT+4', 'Etc/GMT+5') 

	# Summer (May-Oct)
	these_dates = strptime(paste(these_loggers$Date, these_loggers$Time), '%Y-%m-%d %H:%M:%S', tz='Etc/GMT+4') # Loggers were initially deployed during DST
	begin = strptime(paste(this_samp$Logger_start_date, this_samp$Logger_start_time), '%Y-%m-%d %H:%M:%S', tz='Etc/GMT+4') # Loggers were started during DST
	end = strptime(paste(this_samp$Logger_dump_date, this_samp$Logger_dump_time), '%Y-%m-%d %H:%M:%S', tz=wintz) 

	keep_loggers_sum = these_loggers[these_dates>begin & these_dates<end,]
		
	# Winter (Nov-Apr)
	these_dates = strptime(paste(these_loggers$Date, these_loggers$Time), '%Y-%m-%d %H:%M:%S', tz=wintz) # Some loggers were re-deployed during DST others were not
	begin = strptime(paste(this_samp$Logger_dump_date, this_samp$Logger_dump_time), '%Y-%m-%d %H:%M:%S', tz=wintz) 
	end = strptime(paste(this_samp$Logger_end_date, this_samp$Logger_end_time), '%Y-%m-%d %H:%M:%S', tz='Etc/GMT+4') # Loggers were read during DST

	keep_loggers_win = these_loggers[these_dates>begin & these_dates<end,]

	keep_loggers = rbind(keep_loggers_sum, keep_loggers_win)
	keep_loggers = keep_loggers[!duplicated(keep_loggers),]

	good_loggers = rbind(good_loggers, keep_loggers)

}

# Add column that indicates time zone for each logger measurement
good_loggers$Time_zone = 4




## Convert times to UTC
times_gmt4 = strptime(paste(good_loggers$Date, good_loggers$Time), '%Y-%m-%d %H:%M:%S', tz='Etc/GMT+4')

for(i in 1:72){
	
	this_samp = subset(samples, SampID==i)

	# Calculate time zone for when logger data was dumped in Novemeber
	# In 2014, some loggers were read out before DST ended and some were read out after DST ended
	wintz = ifelse(this_samp$Year==2014 & as.Date(this_samp$Logger_dump_date) < as.Date('2014-11-02'), 'Etc/GMT+4', 'Etc/GMT+5') 

	# If the time-zone changed to GMT-05 when data was dumped then update time zone of winter measurements to reflect this
	if(wintz=='Etc/GMT+5'){

		# Calculate date when data was retrived from loggers
		dump_date = strptime(paste(this_samp$Logger_dump_date, this_samp$Logger_end_time), '%Y-%m-%d %H:%M:%S', tz='Etc/GMT+5')
		dump_date_gmt4 = strptime(dump_date + 60*60, '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-4')
		
		# Identify events logged after the dump date
		these_events = (good_loggers$SampID==i) & (times_gmt4 > dump_date)
		
		# Change time zone for these events
		good_loggers[these_events,'Time_zone'] = 5
	}
}


# Write to file - be careful not to overwrite data.
write.csv(good_loggers, paste(logger_data_dir, 'loggers_2015-04-30.csv', sep=''), row.names=F)
write.csv(loggers, paste(logger_data_dir, 'all_logger_data.csv', sep=''), row.names=F)
