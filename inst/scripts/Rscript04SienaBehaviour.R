################################################################################
###
### ---- Rscript04SienaBehaviour.R: a script for the introduction to RSiena ----
###
###                         version September 8, 2020
################################################################################
#
# The introductory script is divided into the following script files:
# Rscript01DataFormat.R, followed by
# RScriptSNADescriptives.R, code for descriptive analysis of the data, and
# Rscript02SienaVariableFormat.R, that formats data and specifies the model, and
# Rscript03SienaRunModel.R, that runs the model and estimates parameters
# Rscript04SienaBehaviour.R, that illustrates an example of analysing the
# coevolution of networks and behaviour
# Written with contributions by Robin Gauthier, Tom Snijders, Ruth Ripley,
# and Johan Koskinen.
#

# Here is a short script for analysing the co-evolution of the
# friendship network and drinking behaviour for the 50 girls in the
# Teenage Friends and Lifestyle Study data
# (http://www.stats.ox.ac.uk/~snijders/siena/s50_data.zip), described in
# http://www.stats.ox.ac.uk/~snijders/siena/s50_data.htm

# Read in the adjacency matrices, covariates and dependent behavioural variable
# assuming data are in current working directory

#		friend.data.w1 <- as.matrix(read.table("s50-network1.dat")) # network
#		friend.data.w2 <- as.matrix(read.table("s50-network2.dat"))
#		friend.data.w3 <- as.matrix(read.table("s50-network3.dat"))
#		drink <- as.matrix(read.table("s50-alcohol.dat")) # behaviour
#		smoke <- as.matrix(read.table("s50-smoke.dat")) # covariate

# If you wish to make it easier, you can use this data set as included
# in the package - but the above is included to show you
# how to use data from files.
# To use the internal data set:

        friend.data.w1 <- s501
        friend.data.w2 <- s502
        friend.data.w3 <- s503
        drink <- s50a
        smoke <- s50s

# At this point it is a good idea to use the sna package to plot the networks
# and the behavioural variable. Descriptive measures of the similarity of
# "friends" with respect to behaviour (like Moran's I) are given by the function
# nacf() in the sna package. See the script RscriptSNADescriptives.R.

# Tell RSiena that the adjacency matrices are network data and in what order
# they should be treated

	friendship <- sienaDependent( array( c( friend.data.w1, friend.data.w2,
							friend.data.w3 ),
							dim = c( 50, 50, 3 ) ) )# create dependent variable

# Tell RSiena that the variable "drink" should be treated
# as a dependent variable

		drinkingbeh <- sienaDependent( drink, type = "behavior" )
		smoke1 <- coCovar( smoke[ , 1 ] )

# Define the data set and obtain the basic effects object
		myCoEvolutionData <- sienaDataCreate( friendship, smoke1, drinkingbeh )
		myCoEvolutionEff <- getEffects( myCoEvolutionData )

# Run reports to check that data is properly formated and
# to get some basic descriptives

        print01Report( myCoEvolutionData, modelname = 's50_3_CoEvinit' )

# Define the effects to include in the coevolution model
# Start with some structural effects (use the shortnames that you find in
# effectsDocumentation(myCoEvolutionEff) )

		myCoEvolutionEff <- includeEffects( myCoEvolutionEff, transTrip, cycle3)

# Include a homophily effect for the constant covariate smoking

		myCoEvolutionEff <- includeEffects( myCoEvolutionEff, simX,
											interaction1 = "smoke1" )

# If we want to parse out whether there is a selection or influence (or both)
# effect for drinking behaviour,
# we need to also include sender, receiver and homophily effects
# of drinking for friendship formation:

		myCoEvolutionEff <- includeEffects(myCoEvolutionEff, egoX, altX, simX,
										   interaction1 = "drinkingbeh" )

# For the influence part, i.e. the effect of the network on behaviour,
# we specify the following effects:
# indegree, outdegree and assimilation effects for drinking

		myCoEvolutionEff <- includeEffects( myCoEvolutionEff,
								name = "drinkingbeh",
								avAlt,indeg, outdeg,
								interaction1 = "friendship" )

# Check what effects you have decided to include:

		myCoEvolutionEff

#
# Now we have to define the algorithm settings.
# The defaults are adequate. You only have to specify the filename
# that will receive the results in text format.

		myCoEvAlgorithm <- sienaAlgorithmCreate( projname = 's50CoEv_3' )

# Finally, estimate the model; the whole command is put in parentheses
# to have the results printed directly to the screen.

        (ans <- siena07( myCoEvAlgorithm, data = myCoEvolutionData,
                        effects = myCoEvolutionEff ))

# THE RESULTS

# Note that the "convergence t-ratio" is the t-ratio for convergence checking,
# not the t statistic for testing the significance of this effect.
# (See Section 6.1.2 of the manual.)
# For good convergence, the t-ratios for convergence
# all should be less than .1 in absolute value,
# and the overall maximum convergence ratio should be less than 0.25.
# If this is not yet the case, you should try again, starting from the
# last estimate as the "previous answer":

        (ans1 <- siena07( myCoEvAlgorithm, data = myCoEvolutionData,
                        effects = myCoEvolutionEff, prevAns = ans ))

# which can be repeated if necessary.

# For this small data set, the model for behavior dynamics is over-specified,
# leading to some very large standard errors.
# For this data set it is better to drop the degree effects on behaviour,
# because the data does not contain enough information to estimate them.

        myCoEvolutionEff2 <- includeEffects( myCoEvolutionEff,
								name = "drinkingbeh", indeg, outdeg,
								interaction1 = "friendship", include = FALSE)

       (ans2 <- siena07( myCoEvAlgorithm, data = myCoEvolutionData,
                       effects = myCoEvolutionEff2 ))

###############################################################################
##                              Some other effects                           ##
###############################################################################
#
# The set of available effects for this data set can be obtained by requesting
        effectsDocumentation(  myCoEvolutionEff )
# See the manual, Chapter 12, for the meaning of these effects.
# To study the direct effect of the actor covariate smoking on the dependent
# variable drinking, use the effFrom effect:

        myCoEvolutionEff3 <- includeEffects( myCoEvolutionEff2,
                                name = "drinkingbeh", effFrom,
                                interaction1 = "smoke1")

# Since we already have a good result for a simpler model,
# it is helpful to start estimating from these estimates
# as starting values:

       (ans3 <- siena07( myCoEvAlgorithm, data = myCoEvolutionData,
                       effects = myCoEvolutionEff3, prevAns = ans2 ))
# In my case, convergence was not good enough:
       (ans3 <- siena07( myCoEvAlgorithm, data = myCoEvolutionData,
                       effects = myCoEvolutionEff3, prevAns = ans3 ))
# You can get a nicer presentation of the results in a file
# in your working directory in LaTeX by
        siena.table(ans3)
# and in html (can be imported into MS-Word) by
        siena.table(ans3, type="html")


