# This program takes as input decimal coordinates of a place of origin
# along with a file from http://download.geonames.org/export/dump/.
# This can be a single country file or several such files concatenated.
# An example of how to run it follows, requiring the AE.txt country file.
# Save this program in the same folder as the geonames file, open it,
# and read it in (press edit > run all). Now, in the console, write
#      IT(25, 55, 24.5, 56, "AE.txt")
# This will produce a route between two locations in the United Arab Emirates
# and the total distance.
# The route is saved in a file called route.txt. The distance is printed to
# the console but can also be made an object in R, called, say, d if
# this is indicated by running the program as d <- migrate(...).
# If a different name for the output file is desired this can be 
# specified after a comma: migrate(..., outputfile="other.txt")

IT <- function(Ala,Alo,Bla,Blo,geonames_file,outputfile="route.txt") {
	library(argosfilter)
	# read the geonames file and reduce the resulting table to those rows
	# for which col. 7 has a "P" and then to the columns for lats and longs
	x <- read.table(file=geonames_file,header=FALSE, sep="\t", quote="", stringsAsFactors=FALSE, encoding="UTF-8", comment.char = "")
	P <- which(x[,7]=="P")
	x <- x[P,c(5,6)]
	# get rid of repeated locations; cp is coordinates pasted together
	cp <- paste(x[,1], "_", x[,2], sep="")
	dup <- which(duplicated(cp))
	if ( length(dup) > 0 ) {
		x <- x[-dup,]
	}
	# the program assumes that the origin coordinates are already in x, so
	# the following subroutine adds them if they are not
	where_Ala <- which(Ala==x[,1])
	where_Alo <- which(Alo==x[,2])
	if ( length(intersect(where_Ala,where_Alo))==0 ) { 
		x <- rbind(c(Ala,Alo),x)
	}
	where_Bla <- which(Bla==x[,1])
	where_Blo <- which(Blo==x[,2])
	if ( length(intersect(where_Bla,where_Blo))==0 ) { 
		x <- rbind(x,c(Bla,Blo))
	}
	# start a matrix with the origin in the first row
	route <- matrix(c(Ala,Alo),ncol=2)
	# make the origin a station
	Sla <- Ala; Slo <- Alo
	# measure the GCD from origin to goal
	orig_dist <- distance(Ala,Bla,Alo,Blo)
	# run the route-finding algorithm
	while ( orig_dist > 0 ) {
		N <- 2
		repeat {
			dist <- as.vector(apply(x, 1, function(z) distance(z[1], Sla, z[2], Slo)))
			nearest <- x[order(dist)[N],]
			orig_dist <- distance(Sla,Bla,Slo,Blo)
			if (orig_dist==0) break
			new_dist <- distance(nearest[1,1],Bla,nearest[1,2],Blo)
			if ( orig_dist - new_dist <= 0 ) {
				N <- N + 1
			} else {
			Sla <- nearest[1,1]; Slo <- nearest[1,2]
			route <- rbind(route, c(Sla,Slo))
			break
			}
		}
	}
	# output the route
	write.table(route, file=outputfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	# sum up the GCDs in the route and round off to integers
	Km <- 0
	for (i in 1:(length(route[,1])-1)) {
		Km <- Km + distance(route[i,1],route[i+1,1],route[i,2],route[i+1,2])
	}
	Km <- round(Km)	
	return(Km)
}
