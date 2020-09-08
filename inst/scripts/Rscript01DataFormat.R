################################################################################
###                                                                          ###
### ---- Rscript01DataFormat.R: a script for the introduction to RSiena ---- ###
###                                                                          ###
###                version: September 8, 2020                                ###
################################################################################
#
# Rscript01DataFormat.R is followed by
# RScriptSNADescriptives.R, code for descriptive analysis of the data, and
# Rscript02SienaVariableFormat.R, which formats data and specifies the model,
# and Rscript03SienaRunModel.R, which runs the model and estimates parameters
# Rscript04SienaBehaviour.R, which illustrates an example of analysing the
# coevolution of networks and behaviour
#
# The entire model fitting is summarised at the end of RscriptSienaRunModel.R
# (without comments)
#
# This is an R script for getting started with RSiena, written by
# Tom Snijders, with earlier contributions from Robin Gauthier,
# Ruth Ripley, Johan Koskinen, Paulina Preciado, and Zsofia Boda,
# with some examples borrowed from Christian Steglich.
# Lines starting with # are not processed by R but treated as comments.
# The script has a lot of explanation of R possibilities that will be
# familiar for readers well acquainted with R, and can be skipped by them.
#
# There are various really easy online introductions to R. See, for example
#
#    http://www.statmethods.net/
#    http://www.burns-stat.com/pages/Tutor/hints_R_begin.html
#    http://data.princeton.edu/R/gettingStarted.html
#    https://stats.idre.ucla.edu/r/
#
# You can go to any of these sites to learn the basics of R
# or refresh your knowledge.
# There is a lot of documentation available at
#    https://www.r-project.org/other-docs.html
# including some short introductions, handy reference cards,
# and introductions in many languages besides English.
#
#
# Some general points to note:
#
# R is case-sensitive. Be aware of capitalization!
#
# The left-arrow "<-" is very frequently used: it denotes an assignment,
# "a <- b" meaning that object a gets the value b.
# Often b is a complicated expression that has to be evaluated by R,
# and computes a result that then is stored as the object a.
#
# Help within R can be called by typing a question mark and the name of the
# function you need help for. For example
# ?library
# will bring up a file titled "loading/attaching and listing of packages".
# Comments are made at the end of commands after #,
# or in lines starting with # telling R to ignore everything beyond it.
# That is why everything up to now in this file is on lines starting with #.
#
# This session will be using the s50 data which are available in RSiena.
# You can also get them from
#    http://www.stats.ox.ac.uk/~snijders/siena/siena_datasets.htm
#
# Any command in R is a function, and ends by parentheses
# that enclose the arguments of the function,
# or enclose nothing if no argument is needed, such as for the function q()
# In general the command syntax for calling R's functions is function(x) where
# function is an available function and x the name of the object operated on.
#
#################### - CALLING THE DATA AND PRELIMINARY MANIPULATIONS - ########

# The library command loads the packages needed during the session.

    library(RSiena)

# Some additional packages are used by RSiena,
# the so-called required packages; these will be loaded automatically.

# You need to have INSTALLED all of them.

    ?install.packages

# Or click on the tab "Packages", "Install package(s)", then select a
# CRAN mirror close to you (e.g. Bristol if you are in the UK) and finally
# select from the list the package you wish to install.

# Where are you?

    getwd()

# By something like
    # setwd('C:/SienaTest')
# you can set the directory but note the quotes and forward slash.
# It is also possible to set the directory using the menus if you have them.
# On a windows machine, you can predetermine the working directory
# in the <Properties> of the shortcut to R;
# these are accessible by right-clicking the shortcut.

# If you want to set the working directory to for example
    # "C:\Documents and Settings\johan\My Documents\RSiena_course"
# simply copy and paste from windows explorer or type
#    setwd('C:/Documents and Settings/johan/My Documents/RSiena_course')
# or
#    setwd('C:\\Documents and Settings\\johan\\My Documents\\RSiena_course')
#
# but note that "\" has to be changed; both '/' and '\\' work!!!

# What is in this directory?

    list.files()

# What is available in RSiena?

    ?RSiena

# (these are .htlm HELP PAGES)

# At the bottom of this page, when you click on "Index",
# a list of all the available functions is shown in your browser.
# The same list is shown in the graphical user interface for R by requesting

    library(help=RSiena)

# An extensive manual is at the Siena website
# at http://www.stats.ox.ac.uk/~snijders/siena/RSiena_Manual.pdf ;
# it is updated frequently.

# Each data is named (for example below we name it friend.data.w1)
# so that we can call it as an object within R.
# If you read an object straight into R, it will treat it as a
# dataset, or in R terminology a "data frame".
# Here this is not what we want, therefore on reading
# we will immediately convert it to a matrix.
# R will read in many data formats, the command to read them is read.table.
      ?read.table
# In the help page for read.table, look at the section "Value",
# which is there in every help page:
# it indicates the class of the object that is produced by the function.
# For read.table, the value is a data frame; below we see what this is.
#
# If we wished to read a .csv file we would have
# used the read.csv command.
# The pathnames must have forward slashes, or double backslashes.
# If single backslashes are used, one of the error messages will be:
#   1: '\R' is an unrecognized escape in a character string

#-----------------------------------------------------------------------------

# Quick start (Data assignment).
# Please make sure the s50 data set is in your working directory.
# The data set is on the Siena website ("Data sets" tab) and must be
# unzipped in your working directory.

#    friend.data.w1 <- as.matrix(read.table("s50-network1.dat"))
#    friend.data.w2 <- as.matrix(read.table("s50-network2.dat"))
#    friend.data.w3 <- as.matrix(read.table("s50-network3.dat"))
#    drink <- as.matrix(read.table("s50-alcohol.dat"))
#    smoke <- as.matrix(read.table("s50-smoke.dat"))

# To make this script self-contained, you may instead run the commands:
    library(RSiena)
    friend.data.w1 <- s501
    friend.data.w2 <- s502
    friend.data.w3 <- s503
    drink <- s50a
    smoke <- s50s
# because s501 etc. are internal data objects of RSiena.

# Explanation of data structures and formats below

################# - DIFFERENT DATA FORMATS - ###################################
###
### The assignments here involve reading in data to a "data frame"
###    data <- read.table("data.dat")
### reads your text file into a "data frame";
### check the class of the object "data"
###    >  class(data)
###    [1] "data.frame"
### If your data is not a ".dat" file, alternative import methods are
### ----- ".csv" ---------------------------------------------------------------
###    data <- read.table("c:/data.csv", header=TRUE,
###                sep=",", row.names="iden")
### ---- ".csv" ----------------------------------------------------------------
###    data <- read.csv("data.csv",header=T)
### look at ?read.csv and ?read.table to understand these function calls
### ---- ".spss" ---------------------------------------------------------------
###    library(foreign)
###    data <- read.spss("c:/data.spss")
### ---- ".dta" ----------------------------------------------------------------
###    library(foreign)
###    data <- read.dta("c:/data.dta")
###
################################################################################

################# - FROM DATA FRAME TO MATRIX - ################################
###
### ---- Data frame ------------------------------------------------------------
### The data frame is like a spreadsheet of cases by variables,
### where the variables are the columns, and these have the names
###    names( data )
### a spreadsheet view of a data frame is given by the fix() command
###    fix( data )
###
### All changes you make in this spreadsheet will be saved automatically,
### so beware.
###
### As an example create two vectors:
    height <- c( 120, 150, 180, 200, 190, 177, 170, 125, 141, 157 )
    weight <- c( 11, 14, 17, 18, 17, 18, 11, 12, 10, 15 )
### The function c() combines its argument into a vector
### (or into a list, but we are not concerned with that possibility now.).
### These two vectors can be collected as variables in a data frame
    mydata <- data.frame( height, weight )
### and look at the results
    mydata
### The columns of a data frame may be extracted
### using a "$" sign and their names.
### For example:
    names( mydata )
    mydata$height
### Or by "[]" and column number, e.g.
    mydata[1]
    mydata[  , 1]
    mydata[ 1,  ]
### If you wish to see the structure of an object, such as mydata, then request
    str(mydata)
### Objects often have attributes:
    attributes(mydata)

### ---- Matrix ----------------------------------------------------------------
### A "matrix" is a numeric object ordered like a matrix with dimensions
### ( dim() ) given by the number of rows and columns, e.g.
    dim( friend.data.w1 )

### If you wish to play around with a copy of the matrix,
###e.g. having the name "mydata", you can make the copy by the command
    mydata <- friend.data.w1
### The earlier object that we created with the name "mydata" now has been lost.
### Elements if matrices can be accessed by using square brackets.
### For example, the element of "mydata" in row 2 and column 3 is given by
    mydata[ 2, 3 ]
### the first three rows of a matrix called mydata are given by
    mydata[ 1:3, ]
### columns 2, 5, and 6 are given by
    mydata[ , c( 2, 5, 6) ]
### Indeed there are a lot of zeros - networks tend to be sparse.

### ---- Converting data frame to matrix ---------------------------------------
###
### Most classes can be converted to other classes through "as.typeofclass ()",
### e.g., if "mydata" would be an object with a matrix-like structure,
### then it could be converted to the class "matrix" by the command
    mydata <- as.matrix( mydata )

################################################################################


################## - EXAMPLE FOR ARC LIST - ####################################
###
### From the Siena website you can download the data set arclistdata.dat.
### The function read.table can be used with the internet location:
    ArcList <-
      read.table( "https://www.stats.ox.ac.uk/~snijders/siena/arclistdata.dat",
                header=FALSE )
### Note the capitalization.
### Now ArcList is a data.frame
### (we saw this above in the help page for read.table).
### What are its dimensions?
       dim(ArcList)
### You can get an impression of it by looking at the start of the object:
       head(ArcList)
### This is a data file in arclist format, with columns:
### sender id, receiver id, value of tie, wave.
### Add names to the columns:
    names(ArcList) <- c( "sid", "recid", "bff", "wid" )
### and again request
    head(ArcList)
### The bff ("best friend") variable does not have much variability:
    table(ArcList$bff)
### This tells us that non-ties are not included (they would have the value 0),
### and that there are no tie values other than 1.
### It may be nice to order the rows by sender, then by receiver, then by wave.
### The tedious way to do this is
    ArcList <- ArcList[ order( ArcList$sid, ArcList$recid, ArcList$wid), ]
### To understand this, you may look at the help page for function order()
### The rows of Arclist are ordered first by ArcList$sid, then by ArcList$recid,
### and then by ArcList$wid;
### these reordered rows then are put into the object ArcList,
### thus overwriting the earlier contents of this object.
### This way is tedious because it is repeated all the time that the names
### sid, recid, wid refer to the ArcList object.
### The less tedious way uses the function "with".
### The with(a, b) function tells R that b must be calculated,
### while the otherwise unknown names refer to data set a:
    ArcList <- with( ArcList, ArcList[ order( sid, recid, wid), ] )
### An arc list does not give information about the number of nodes,
### because isolates are not recorded.
### The set of non-isolates is
    union(unique(ArcList$sid), unique(ArcList$recid))
### and, given that we have the information that there are 50 nodes
### labeled 1 to 50, the isolates are the following two nodes:
    setdiff(1:50, union(unique(ArcList$sid), unique(ArcList$recid)))
### Now suppose we want to create a separate set of records for each wave.
### Select by wave:

    SAff.1 <- with(ArcList, ArcList[ wid == 1, ] ) #extracts edges in wave 1
    SAff.2 <- with(ArcList, ArcList[ wid == 2, ] ) #extracts edges in wave 2
    SAff.3 <- with(ArcList, ArcList[ wid == 3, ] ) #extracts edges in wave 3

### This can be arranged more efficiently as
    SAff <- lapply(1:3, function(m){ with(ArcList, ArcList[ wid == m, ] ) } )
### with the correspondence that SAff.1 is SAff[[1]], etc.

### For transforming an adjacency matrix, e.g., friend.data.w1, into an arclist,
### create indicator matrix of the non-zero entries:
    ones <- !friend.data.w1 %in% 0
### create empty edge list of desired length
    edges <- matrix(0, sum(ones), 3)
### fill the columns of the edge list
    edges[, 1] <- row(friend.data.w1)[ones]
    edges[, 2] <- col(friend.data.w1)[ones]
    edges[, 3] <- friend.data.w1[ones]
### if desired, order edge list by senders and then receivers
    edges <- edges[order(edges[, 1], edges[, 2]), ]
### What did we make:
    edges

### For transforming an arclist into a matrix,
### first remove the fourth column indicating the wave, so that we are left
### with sender, receiver and value of the tie,
### and transform it into matrix format (at first it is a data.frame)
    SAff.1.copy <- SAff.1[, 1:3]
    SAff.1.copy <- as.matrix(SAff.1.copy)
### create empty adjacency matrix
    adj <- matrix(0, 50, 50)
### put edge values in desired places
    adj[SAff.1.copy[,1:2]] <- SAff.1.copy[, 3]
### Now adj is the adjacency matrix.
### Also see Section 4.1 of the Siena manual
### and the help page for sienaDependent.

################################################################################
################ - READING IN PAJEK DATA - #####################################
###
### Skip this section if handling Pajek data is of no concern to you.
### If you have data in Pajek format you can use the package "network" in order
### to convert it to a network object. This example is from ?read.paj
###   require( network )
###
###   par( mfrow = c( 2, 2 ) )
###
###   test.net.1 <-
###    read.paj("http://vlado.fmf.uni-lj.si/pub/networks/data/GD/gd98/A98.net" )
###   plot( test.net.1,main = test.net.1$gal$title )
###   test.net.2 <-
###    read.paj("http://vlado.fmf.uni-lj.si/pub/networks/data/mix/USAir97.net" )
###   plot( test.net.2,main = test.net.2$gal$title )
###
################################################################################

# Before we work with the data, we want to be sure it is correct. A simple way
# to check that our data is a matrix is the command class()

    class( friend.data.w1 )

# to list the properties of an object attributes( friend.data.w1 )
# (different classes have different attributes)

# To check that all the data has been read in, we can use the dim() command.
# The adjacency matrix should have the same dimensions as the original data
# (here, 50 by 50).

    dim(friend.data.w1)
    dim(drink)

# To check the values are correct, including missing values, we can use
# the following commands to tabulate the variables.

    table( friend.data.w1, useNA = 'always' )
    table( friend.data.w2, useNA = 'always' )
    table( friend.data.w3, useNA = 'always' )
    table( drink, useNA = 'always' )
    table( smoke, useNA = 'always' )

# NA is the R code for missing data (Not Available).
# This data set happens to have no missings
# (see the data description on the Siena website).
# If there are any missings,
# it is necessary to tell R about the missing data codes.
# Let us do as if the missing codes for the friendship network were 6 and 9.
# This leads to the following commands.
# (For new R users: the c() function used here as "c(6,9)" constructs
#  a vector [c for column] consisting of the numbers 6 and 9.
#  This function is used a lot in basic R.)

    friend.data.w1[ friend.data.w1 %in% c(6,9) ] <- NA
    friend.data.w2[ friend.data.w2 %in% c(6,9) ] <- NA
    friend.data.w3[ friend.data.w3 %in% c(6,9) ] <- NA

# Commands for descriptive analysis are in the script: RscriptSNADescriptives.R

############## - SELECTING SUBSETS OF DATA - ###################################

# To select a subset of the data based on an actor variable, say,
# those who have the value 2 or 3 on drinking at time 1
# (the possibilities are endless, but hopefully this will serve as a pattern)

     use <- drink[, 1] %in% c(2, 3)

# This creates a logical vector which is TRUE for the cases where the condition
# is satisfied. To view or check, display the vectors next to each other:

    cbind(drink[ , 1 ], use)
    data.frame(drink[ , 1 ], use)

# and the number of selected cases is displayed by

     sum( use )

# or

    table( use )

# Given this selection, submatrices can be formed in case the analyses
# are to be done for this subset only:

    friend1.data.w1 <- friend.data.w1[ use, use ]
    friend1.data.w2 <- friend.data.w2[ use, use ]
    drink1 <- drink[ use, ]

# A useful option in R that allows you to save your workspace:
save.image("WorkspaceRscript01.RData")
# Later you can load this in a new session by
# load("WorkspaceRscript01.RData")

# If next time you would like to continue from here, you will not need to open
# and run this script again, since you will be able to load this current state
# of your workspace. But packages will have to be attached again.

################################################################################
###
### ---- PROCEED TO RscriptSNADescriptives.R FOR DESCRIPTIVE ANALYSIS ----------
###
################################################################################
