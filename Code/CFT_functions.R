## This script holds functions used in the Canopy Functional Traits project



# A function that plots a vertical color ramp on the side of a plot
# cols    : the colors to use
# n       : number of divisions
# barends : location of whole bar c(xleft, ybottom, xright, ytop)
# labels    : vector of labels for bar, assumes 1st and last numbers correspond to 1st and last colors
# title   : title to print above bar
plotColorRamp = function(cols, n, barends, labels=NA, title=NA, mycex=1.5){
	dX = barends[3] - barends[1]
	dY = barends[4] - barends[2]
	dy = dY/n
	
	xpd.old = par('xpd')
	par(xpd=T)

	usecols = colorRampPalette(cols)(n)

	for(i in 1:n){
		rect(barends[1], barends[2]+dy*(i-1), barends[3], barends[2]+dy*i, col=usecols[i], border=NA)
	}

	if(!is.na(labels)){
		dZ = labels[length(labels)]-labels[1]
		dz = dY/dZ
		Yposition = barends[2] + dz*(labels-labels[1])

		text(barends[3]+dX*0.5, Yposition, round(labels,2), pos=4, cex=mycex)
		
		segments(barends[3], Yposition, barends[3]+dX*0.5, Yposition)	
	}
	if(!is.na(title)){
		labels.round = round(labels, 2)
		
		## Determine how many characters away to place title
		digits = max(nchar(round(labels, 2))) # Maximum number of digits in a label
		largest = labels.round[which(nchar(labels.round)==digits)] # Which labels are longest
		no.decimal = sum(largest == floor(largest))>0 # Does one largest label lack a decimal?
			if(!no.decimal) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		no.negative = sum(largest >= 0)>0 # Does one largest label lack a negative sign?
			if(!no.negative) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		
		text(barends[3]+dX*0.5+par('cxy')[1]*mycex*(digits+.5), barends[2]+0.5*dY, labels=title, srt=-90, cex=mycex)
	}
	par(xpd=xpd.old)
}



# A functions that plots margin text in the correct orientation
mtexti <- function(text, side, off = 0.25,
                   srt = if(side == 2) 90  else
                         if(side == 4) 270 else 0, ...) {
    # dimensions of plotting region in user units
    usr <- par('usr')
    # dimensions of plotting region in inches
    pin <- par('pin')
    # user units per inch
    upi <- c(usr[2]-usr[1],
             usr[4]-usr[3]) / pin
    # default x and y positions
    xpos <- (usr[1] + usr[2])/2
    ypos <- (usr[3] + usr[4])/2
    if(1 == side)
        ypos <- usr[3] - upi[2] * off
    if(2 == side)
        xpos <- usr[1] - upi[1] * off
    if(3 == side)
        ypos <- usr[4] + upi[2] * off
    if(4 == side)
        xpos <- usr[2] + upi[1] * off
    text(x=xpos, y=ypos, text, xpd=NA, srt=srt, ...)
}


## A function that calculate unique taxa from a sample when samples may contain genera
# x = a vector of TaxonIDs
# taxa = a dataframe with a column for TaxonID and TaxonConcept ('species','genus','unknown')


calc_unique_taxa = function(x, taxa){
	rownames(taxa) = taxa$TaxonID
	these_taxa = taxa[x,]

	species = subset(these_taxa, TaxonConcept=='species')
	genera = subset(these_taxa, TaxonConcept=='genus')

	sp_genera = unique(species$Genus)
	unique_genera = subset(genera, !(Genus %in% sp_genera))

	new_species = c(species$TaxonID, unique_genera$TaxonID)

	new_species
}

## A function that attempts to assing a species name to thalli only IDed to genus based on an additional list of taxa found in the sample
# x = a vector of TaxonIDs to be checked
# others = a vector of other TaxonIDs known to exist in the sample
# taxa = a dataframe with a column for TaxonID and TaxonConcept ('species','genus','unknown') 

resolve_genera = function(x, others, taxa){
	
	x = calc_unique_taxa(x, taxa)

	rownames(taxa) = taxa$TaxonID
	these_taxa = taxa[x,]

	genera = subset(these_taxa, TaxonConcept=='genus')$Genus

	if(length(genera)>0){
		other_taxa = taxa[others,]
		other_species = subset(other_taxa, TaxonConcept=='species')
		
		# Go through each genus	
		for(g in genera){
			potential_species = subset(other_species, Genus==g)
			
			# Only replace the unidentified thallus if there is only one other potential species it could be	
			if(nrow(potential_species)==1){
				which_x = which(taxa[x,'Genus']==g) # Find the TaxonID to replace
				x[which_x] = potential_species$TaxonID
			}
		}
	} 

	# Return the modified TaxonID list
	x	
}


# Define Hellinger distance function that finds 0 distance between genus sp. and species in the same genus.
# i and j are vectors of species abundances
# taxon_dependency is a binary matrix whose columns indicate the taxonomic hierarchy
#	rownames are the taxa that are ambiguous (i.e. the genera)
calc_Hellinger = function(i, j, taxon_dependency){
	# Determine which taxa are ambiguous
	ambig_taxa = which(colnames(taxon_dependency) %in% rownames(taxon_dependency))
	
	# Determine whether there is a lichen only identified to genus
	if(sum(i[ambig_taxa])>0){
	for(n in names(which(i[ambig_taxa]>0))){
		# If the other site does have a species in that genus
		# Distribute individuals among those species according to their abundance distribution at that site
		dependents = names(which(j[taxon_dependency[n,]]>0))

		if(length(dependents)>0){
			Ni = i[n]
			Nj = sum(j[dependents])
			i[dependents] = Ni*j[dependents]/Nj
			if(n!=dependents) i[n] = 0 # Checks whether dependent taxon is actually the same as the ambiguous taxon
		}
	}}

	if(sum(j[ambig_taxa])>0){
	for(n in names(which(j[ambig_taxa]>0))){
		# If the other site does have a species in that genus
		# Distribute individuals among those species according to their abundance distribution at that site
		dependents = names(which(i[taxon_dependency[n,]]>0))

		if(length(dependents)>0){
			Nj = j[n]
			Ni = sum(i[dependents])
			j[dependents] = Nj*i[dependents]/Ni
			if(n!=dependents) j[n] = 0
		}
	
	}}
	
	# Calculate Hellinger distance
	sqrt(sum((sqrt(i/sum(i)) - sqrt(j/sum(j)))^2))
}


# A function that returns the y-coordinates of rectangles for plotting aov decomposition
# Assumes numbers express relative components
# x is a vector of 3 numbers: (fixed average, intraspecific varaibility, covariation)
aov_heights = function(x){
	y_mat = matrix(NA, nrow=2, ncol=3)
	colnames(y_mat) = names(x)

	if(x[3] >= 0){
		y_mat[,1] = c(x[3] + x[2], sum(x))
		y_mat[,2] = c(x[3], x[3] + x[2])
		y_mat[,3] = c(0, x[3])		
	} else {
		y_mat[,1] = c(x[3] + x[2], sum(x))
		y_mat[,2] = c(0, x[3] + x[2])
		y_mat[,3] = c(x[3], 0)
	}
	y_mat
}