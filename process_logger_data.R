### This script reads in HOBO data logger CSV files.
### and appends them to the logger table in the Canopy Lichens SQL database

options(stringsAsFactors=F)

# Set up directory and file names
logger_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Canopy Functional Traits/DataLoggers/Nov2013-Apr2014' #'May-Nov2013'
sql_dir = 'Z:/NC Lichens/Canopy Lichens/SQLite'
sql_samples_filename = 'samples.csv'

# If there is an existing file with logger data, then read it in to be appended to
existing_logger_file = F
sql_logger_filename = 'loggers.csv'

# Set up which sample year this data is from
year = 2013

# Find loggerfiles in the logger directory
logger_filenames = list.files(logger_dir, pattern='*.csv', full.name=T)

# Read in existing data from SQL table
sql_samples = read.csv(paste(sql_dir, sql_samples_filename, sep='/'))

if( existing_logger_file){ 
	sql_loggers = read.csv(paste(sql_dir, sql_logger_filename, sep='/'))
} else {sql_loggers = data.frame()}


for( f in logger_filenames){
	this_data = read.csv(f, skip=1, row.names=1)

	# Which data logger is this?
	f_grep = regexpr('log[0-9]+', f)
	this_logID = substr(f, f_grep+3, f_grep+attr(f_grep, 'match.length')-1)

	# Which sample does this correspond to?
	this_sampID = subset(sql_samples, LoggerID==this_logID&Year==year)$SampID
	
	# Extract data to be added to sql_loggers file
	add_data = data.frame(SampID = rep(this_sampID, nrow(this_data)), this_data[,1:3])
	colnames(add_data) = c('SampID','Date_time','Temp','Lux')

	# Remove ending rows that do not correspond to environmental measurements
	if(sum(is.na(add_data$Temp))>0) add_data = add_data[-which(is.na(add_data$Temp)),]

	# Put date/time into different format: YYYY-MM-DD HH:MM:SS
	date_time = strptime(add_data$Date_time, format='%m/%d/%y %I:%M:%S %p')
	add_data$Date = as.character(format(date_time, format='%Y-%m-%d'))
	add_data$Time = as.character(format(date_time, format='%H:%M:%S'))
	
	# Append this data to rfull datat
	sql_loggers = rbind(sql_loggers, add_data[,c('SampID','Date','Time','Temp','Lux')])
}

sql_loggers = data.frame(Index = 1:nrow(sql_loggers), sql_loggers)

# Write to file - be careful not to overwrite data.
write.csv(sql_loggers, paste(sql_dir, sql_logger_filename, sep='/'), row.names=F)

