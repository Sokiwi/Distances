# This requires the following packages to have been installed:
# igraph, deldir, argosfilter, mapproj.
# IT takes as one of its inputs a country file from
# http://download.geonames.org/export/dump/, so files for the countries within
# which the user wishes to get distances should be in the working directory.
# Users can concatenate files, for instance files pertaining to a whole
# continent, but the format should be as in the individual files.
# Read in the function (Edit > Select all).
# Run it by typing IT(Ala,Alo,Bla,Blo,filename), where Ala and Alo are latitude
# and longitude of the origin, Bla and Blo lat and long of the destination,
# and filename is  the desired country file, for instance:
# DDall(22.75, 53.16667, 24.1, 51.4, "AE.txt")
# Coordinates should be in decimal degrees.
# Output is the distance rounded off to an integer.

library(igraph) # for shortest path
library(deldir) # for Delaunay
library(argosfilter) # for computing great circle distances
library(mapproj) # for working with spatial data in general

DDall <- function(Ala,Alo,Bla,Blo,geonames_file) {
	x <- read.table(file=geonames_file,header=FALSE, sep="\t", quote="", stringsAsFactors=FALSE, encoding="UTF-8", comment.char = "")
	P <- which(x[,7]=="P")
	x <- x[P,c(5,6)]
	cp <- paste(x[,1], "_", x[,2], sep="")
	dup <- which(duplicated(cp))
	if ( length(dup) > 0 ) {
		x <- x[-dup,]
	}
	x <- as.matrix(x)
	colnames(x) <- NULL; rownames(x) <- NULL
	add <- matrix(0, nrow=2, ncol=2)
	x <- rbind(add,x)
	# start and goal are deleted from x if they are there
	# and then put in the first and second row of x
	where_Ala <- which(Ala==x[,1])
	where_Alo <- which(Alo==x[,2])
	whereA <- intersect(where_Ala,where_Alo)
	if ( length(whereA) > 0 ) { 
		x <- x[-whereA,]
	}
	where_Bla <- which(Bla==x[,1])
	where_Blo <- which(Blo==x[,2])
	whereB <- intersect(where_Bla,where_Blo)
	if ( length(whereB) > 0 ) { 
		x <- x[-whereB,]
	}
	x[1,1] <- Ala; x[1,2] <- Alo; x[2,1] <- Bla; x[2,2] <- Blo
	xp <- mapproject(x[,2], x[,1], "mercator")
	delaunay <- deldir(xp$x, xp$y)
	e <- delaunay$delsgs
	L <- length(e[,1])*2
	v <- rep(NA,L)
	a <- seq(1, (L-1),2)
	b <- seq(2, L, 2)
	v[a] <- e[,5]; v[b] <- e[,6]
	w <- rep(NA, L/2)
	for (i in 1:(L/2)) {
		w[i] <- distance(x[e[i,5],1], x[e[i,6],1], x[e[i,5],2], x[e[i,6],2])
	}
	g <- make_graph(v, directed = FALSE)
	s <- shortest.paths(g, 1, 2, algorithm="dijkstra", weights=w)
	return(round(as.vector(s)))
}


