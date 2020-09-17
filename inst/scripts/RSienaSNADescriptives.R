##################################################################################
###
### ---- RscriptSNADescriptives.R: a script for the introduction to RSiena -------
###
###                               version: September 8, 2020
##################################################################################
#
# Rscript01DataFormat.R is followed by
# RScriptSNADescriptives.R, code for descriptive analysis of the data, and
# Rscript02SienaVariableFormat.R, which formats data and specifies the model, and
# Rscript03SienaRunModel.R, which runs the model and estimates parameters
# Rscript04SienaBehaviour.R, which illustrates an example of analysing the
# coevolution of networks and behaviour
#
# The entire model fitting is summarised at the end of Rscript03SienaRunModel.R
# (without comments)
#
# This is an R script for getting started with RSiena, written by
# Robin Gauthier, Tom Snijders, Ruth Ripley, Johan Koskinen, and
# Paulina Preciado, with some examples borrowed from Christian Steglich.
# Lines starting with # are not processed by R but treated as comments.
# The script has a lot of explanation of R possibilities that will be
# familiar for readers well acquainted with R, and can be skipped by them.

# A visual inspection of the adjacency matrices can sometimes be useful.
# This will, for example, help in highlighting outliers with respect to
# outdegrees or indegrees, if there are any of such outliers.
# This requires package "sna":

        library( network )
        library( sna )

# For this script, you will need the data read and modified in the script
# Rscript01DataFormat.R. If you have already ran that script, you may
# load the required workspace:
# load("WorkspaceRscript01.RData")

# If not, to make this script self-contained, you may run the commands:

    library(RSiena)
    friend.data.w1 <- s501
    friend.data.w2 <- s502
    friend.data.w3 <- s503
    drink <- s50a
    smoke <- s50s

    net1 <- as.network( friend.data.w1 )
    net2 <- as.network( friend.data.w2 )
    net3 <- as.network( friend.data.w3 )
    plot.sociomatrix( net1,drawlab = F, diaglab = F, xlab = 'friendship t1' )
    plot.sociomatrix( net2,drawlab = F, diaglab = F, xlab = 'friendship t2' )
    plot.sociomatrix( net3,drawlab = F, diaglab = F, xlab = 'friendship t3' )

# The class,
    class( net1 )
# with attributes
    attributes( net1 )
# has special methods associated with it.
# while  plot( friend.data.w1 ) only produces a rather dull plot of
# the first two columns
#     plot( net1, xlab = 'friendship t1' )
# produces a nice sociogram
#
# Some further descriptives you can do for the data are plotting and
# calculating some statistics.

# add the attribute drink to the network object

    net1 %v% "drink" <- drink[ , 1 ]

# color the nodes by drink

    plot( net1, vertex.col = "drink", xlab = 'friendship t1' )


# Now let's color the nodes by drink and scale the vertex by degree of nodes!
#
# First calculate degree:

    deg <- rowSums( as.matrix( net1 ) )# NB:  rowSums() is defined for class matrix

# have a look at the degree distribution

    table( deg, useNA = 'always' )

# Now do the desired plot:

    plot( net1, vertex.col = "drink", vertex.cex = (deg + 1)/1.5 )

# ---- Plot the three waves of data --------------------------------------------

# Add drink to waves 2 and 3
    net2 %v% "drink" <- drink[ , 2 ]
    net3 %v% "drink" <- drink[ , 3 ]
    deg2 <- rowSums( as.matrix( net2 ) )
    deg3 <- rowSums( as.matrix( net3 ) )

# Create a set of panels ( 1 row by 3 columns, or 3 columns by 1 row)

    par( mfrow = c( 1, 3 ) )

# creating three plots after each other will place them in consecutive panels

    plot( net1, vertex.col = "drink", vertex.cex = (deg + 1)/1.5 )
    plot( net2, vertex.col = "drink", vertex.cex = (deg2 + 1)/1.5 )
    plot( net3, vertex.col = "drink", vertex.cex = (deg3 + 1)/1.5 )

# Each time we make a plot the coordinates move - because always
# the starting values are random. We can also save coordinates
# and use them for later plotting:

    par( mfrow = c( 1, 3 ) )
    coordin <-  plot( net1, vertex.col = "drink", vertex.cex = (deg +1 )/1.5 )
    plot( net2, coord = coordin, vertex.col = "drink", vertex.cex = (deg2 + 1)/1.5 )
    plot( net3, coord = coordin, vertex.col = "drink", vertex.cex = (deg3 + 1) /1.5 )

# To get coordinates based on all three waves: coordin <-  plot( net1 + net2 + net3 )
# For more plotting options, try the gplot function in the "sna" library

    ?gplot
    ?gplot.layout


# ---- Basic network statistics ------------------------------------------------

# The package "sna" can be used for a variety of descriptions and analyses.
# The following are examples.
# some important graph level statistics

    gden( net1 ) # density
    grecip( net1 ) # proportion of dyads that are symmetric
    grecip( net1, measure = "dyadic.nonnull" ) # reciprocity, ignoring the null dyads
    gtrans( net1 ) # transitivity

# dyad and triad census

    dyad.census( net1 )
    triad.census( net1 )

# out degree distribution (of course for a symmetric network outdegree=indegree)

    outdegree <- degree( net1, cmode = "outdegree" )
    outdegree #outgoing ties of each note

    hist( outdegree )
    quantile( outdegree )

# measures of connectivity and distance

    dist <- geodist(net1, inf.replace = Inf, count.paths = TRUE)
# calculate the geodesic distance (shortest path length) matrix
    dist$gd
# matrix of geodesic distances
    dist$counts
    table(dist$counts)
# reachability matrix:
    ?reachability
    reach <- reachability( net1 )  # calculate the reachability matrix
    reach


# ---- Network autocorrelation  ------------------------------------------------

# Moran's autocorrelation for outgoing ties:

nacf(net1, drink[, 1], type="moran", neighborhood.type='out')[2]
nacf(net2, drink[, 2], type="moran", neighborhood.type='out')[2]
nacf(net3, drink[, 3], type="moran", neighborhood.type='out')[2]

# Moran's autocorrelation for outgoing and incoming ties:
nacf(net1, drink[, 1], type="moran", neighborhood.type='total')[2]
nacf(net2, drink[, 2], type="moran", neighborhood.type='total')[2]
nacf(net3, drink[, 3], type="moran", neighborhood.type='total')[2]

################################################################################
###
### -- PROCEED TO Rscript02SienaVariableFormat.R FOR PREPARING DATA FOR RSIENA -
###
################################################################################
