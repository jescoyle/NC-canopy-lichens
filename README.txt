Code and Data for "Intraspecific variation in epiphyte functional traits reveals limited effects of microclimate on community assembly in temperate deciduous oak canopies"

Generated by Jessica Coyle
Sept. 15, 2015

http://github.com/jescoyle/NC-canopy-lichens

----------------------------------------------------
Description
This repository contains code and data products for analyses performed in "Intraspecific variation in epiphyte functional traits reveals limited effects of microclimate on community assembly in temperate deciduous oak canopies". This paper examined functional trait variation of macrolichen communities in oak canopies in the Piedmont of North Carolina, USA. Code used in analyses are located in R scripts in the Code directory. Code used in analyses that are no longer a part of the project are stored in directories labeled OLD. Most users will not need to access these files. Data are stored in the Data directory.

---------------------------------------------------
Data

Lichens were collected from 72 quadrat samples in 12 trees located at 2 sites. Traits were measured on macrolichens from trees 1-10, only. All lichens were identified in samples from trees 4-12 only (Duke Forest). Data are stored in a SQLite Database which is located in this repository. Individual csv files that comprise this database are also available. Derived data tables generated by SQL queries and R scripts are also included.

See meta_data.txt for a detailed description of fields in each table.

blue2red_10colramp.txt
	A useful color ramp for plotting.

CanopyLichens.sqlite
	An SQLite databased containing all raw data tables and tables derived from them using SQL queries.

sites.csv
	A table with information about the two forest sites sampled.

trees.csv
	A table with information about each of the 12 trees sampled.

samples.csv
	A table with information about the locations and sampling of each of the 72 quadrat samples.

loggerdata.csv
	A table with summary environmental data for each of the 72 samples derived from datalogger records.

derived_traits.csv
	A table with trait values for each macrolichen thallus.

community.csv
	A table with one record for each lichen species found in each sample for samples collected in Duke Forest (trees 4-12). Includes crustose species.

comm_thalli.csv
	A sample X species matrix with the presence (1) or absence (0) of each taxon (columns) in each sample (rows) for lichen thalli used to conduct trait analyses. Only includes records in derived_traits table.

comm_macro.csv
	A sample X species matrix with the presence (1) or absence (0) of each taxon (columns) in each sample (rows) for all macrolichens in samples collected in Duke forest. Combines records in community and derived_traits tables.


---------------------------------------------------
Code

process_logger_data.R
	This script was used to append individual datalogger records into loggers.csv that contained all datalogger records.

format_logger_data.R
	This script matches datalogger records with weather station records to calculate vapor pressure deficit from temperature logs. It then calculates winter, summer, and annual summary environmental variables for each sample.

CFT_functions.R
	This script contains miscellaneous functions used across multiple scripts in the analysis.

examine_community_structure.R
	This script performs ordination-based analyses of community structure across samples and analyzes individual species distributions using logstic regression.

model_func_traits.R
	This script is used to analyze trait data. Trait distributions are examined and outliers removed. Trait correlations and spatial variation in environmental variables are visualized. Taxonomic trait differences are assessed using ANOVA. Individual traits are modeled using linear mixed effects models. Community mean traits are analyzed using method of trait variance decomposition into interspecific and intraspecific components proposed by Leps et al 2011.

/OLD_Bayesian_models/
	This project orgininally fit multi-level trait models using MCMC estimation of Bayesian models. However a comparison of MLE and Baysian methods for detecting the different sources of trait variance (genus vs. sample) suggested that the Baysian methods were not performing well on simulated data. (See script compare_MLE_BUGS.R). Final analyses were performed using ML analysis of mixed effects models, but older Baysian code is included in the directory for reference. Code in this directory is NOT COMPLETE! Take care if re-using.
