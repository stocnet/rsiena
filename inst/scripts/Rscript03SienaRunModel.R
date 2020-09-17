################################################################################
###                                                                          ###
###                                                                          ###
### --- Rscript03SienaRunModel.R: a script for the introduction to RSiena -- ###
###                                                                          ###
###                         version September 8, 2020                        ###
################################################################################
#
# The introductory script is divided into the following script files:
# Rscript01DataFormat.R, followed by
# RScriptSNADescriptives.R, code for descriptive analysis of the data, and
# Rscript02SienaVariableFormat.R, which formats data and specifies the model,
# Rscript03SienaRunModel.R, which runs the model and estimates parameters, and
# Rscript04SienaBehaviour.R, which illustrates an example of analysing the
# coevolution of networks and behaviour.
# Written by Tom Snijders, with earlier contributions from Robin Gauthier,
# Ruth Ripley, Johan Koskinen, and Paulina Preciado.
#
# This script, Rscript03SienaRunModel.R, runs the estimation in RSiena for the
# model set up and defined in the script Rscript02SienaVariableFormat.R.
#
# A quick version of the model fitting without comments is given at the end
# of this script

########################### ESTIMATION OF PARAMETERS ###########################

## For this script, you will need the data read and modified in the script
# Rscript02SienaVariableFormat.R. If you have already ran that script, you may
# load the required workspace:

	load("WorkspaceRscript02.RData")

# If not, to make this script self-contained, you may run the commands:

	library(RSiena)
	friend.data.w1 <- s501
	friend.data.w2 <- s502
	friend.data.w3 <- s503
	drink <- s50a
	smoke <- s50s
	friendship <- sienaDependent(
                     array( c( friend.data.w1, friend.data.w2, friend.data.w3 ),
                     dim = c( 50, 50, 3 ) ) )
	smoke1 <- coCovar( smoke[ , 1 ] )
	alcohol <- varCovar( drink )
	mydata <- sienaDataCreate( friendship, smoke1, alcohol )

# and request
	mydata
# to see what you have produced.

# Parameters of the model are estimated by the function siena07.
# This requires the data specification; the effects specification;
# and a number of parameters, or settings, for the estimation algorithm.
# The latter are contained in an object created by the function
# sienaAlgorithmCreate. You can look at the help provided by
# ?sienaAlgorithmCreate
# to find out about options that you may use here;
# for beginning users, only the two options mentioned below are relevant.
#
# Output will be written to a file with name projname.txt, where projname is
# whatever name is given; the default (used if no name is given) is Siena.
# This file will be written to your current directory.
# New estimation runs will append to it.
# A new call to print01Report will overwrite it!

	myalgorithm <- sienaAlgorithmCreate(projname = 's50_3')

# Let us first redefine the model, to obtain a simpler specification
# that will serve as an illustration here.

		myeff <- getEffects( mydata )
		myeff <- includeEffects( myeff, transTrip, cycle3 )
		myeff <- includeEffects( myeff, egoX, altX, egoXaltX,
								 interaction1 = "alcohol" )
		myeff <- includeEffects( myeff, simX, interaction1 = "smoke1" )
		myeff

# The function siena07 actually fits the specified model to the data
# If you wish the pretty picture of Siena on the screen as information
# about the progress of the algorithm, type

		ans <- siena07( myalgorithm, data = mydata, effects = myeff)

# (ans for "answer").
# If however you do not want the pretty picture, or if this leads to
# difficulties (which may happen e.g. on a Mac), then type

#		ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)

# and intermediate information will be written to the console.

# Function siena07 produces a so-called sienaFit object, here called ans;
# and it fills in a few things in the sienaEffects object myeff,
# if this is the first use of myeff in a siena07 call.
# By using various different effects objects, i.e., with different names,
# you can switch between specifications.

# The batch = FALSE parameters will give a graphical user interface being opened
# which reports on the progress of the estimation algorithm;

# verbose = TRUE leads to extensive diagnostic information being sent
# to the console during the estimation, and results after the estimation
# (these results are also copied to the output file projname.txt, see above);
# while batch=TRUE gives only a limited amount of printout sent to the console
# during the estimation (which is seen when clicking in the console,
# or more immediately if the Buffered Output is deselected in the Misc menu)
# which monitors the progress of the estimation algorithm in a different way.

# The call of siena07 leads to output in the file s50_3.txt
# (or more generally projname.txt,
# where projname is the name given in sienaAlgorithmCreate)
# and to the creation of the object which here is called ans (for "answer").

# To use multiple processors, e.g., if you wish to use 2 processes, request

#		ans <- siena07( myalgorithm, data = mydata, effects = myeff,
#					  nbrNodes = 2, useCluster = TRUE)

# Adjust the nbrNodes to the number that you wish to use.
# If you wish to work on with other programs while running siena07,
# it is advisable to use one node less than the number of available processors.
# If you wish to use other machines as well,
# see the more detailed instructions below.
# You will then need to use the clusterString argument as well.
#
# For more advanced use, it can be helpful to have access to the networks
# simulated in the so-called third phase of the estimation algorithm.
# These networks can be used, e.g., for checking goodness of fit.
# This can be achieved by using the parameter returnDeps=TRUE.
# The fitted object ans will then have a component named "sims"
# which contains a list (each iteration) of lists (each data object)
# of lists (each dependent network or behavior variable) of edgelists for
# networks or vectors for behavior variables.
# See the manual for further explanation.


################### LOOKING AT THE RESULTS ################################

# The file "s50_3.txt" will contain the results of the estimation.
# It is contained in the current directory ("getwd()").
# This file can be read by any text editor.
# A summary of the results is obtained on the screen by

		ans

# and a larger summary by

		summary(ans)

# Depending on the random seed and the model specification,
# the results could be something like the following.

# Estimates, standard errors and convergence t-ratios
#
#										Estimate   Standard	  Convergence
#													 Error		t-ratio
#
#Rate parameters:
#  0.1      Rate parameter period 1      6.6141  ( 1.1703   )
#  0.2      Rate parameter period 2      5.1578  ( 0.8523   )
#
#Other parameters:
#  1.  eval outdegree (density)         -2.7149  ( 0.1236   )   0.0618
#  2.  eval reciprocity                  2.4114  ( 0.2188   )   0.0412
#  3.  eval transitive triplets          0.6381  ( 0.1491   )   0.0673
#  4.  eval 3-cycles                    -0.0564  ( 0.3012   )   0.0842
#  5.  eval smoke1 similarity            0.2631  ( 0.2078   )   0.0894
#  6.  eval alcohol alter               -0.0217  ( 0.0688   )   0.0399
#  7.  eval alcohol ego                  0.0387  ( 0.0834   )   0.0319
#  8.  eval alcohol ego x alcohol alter  0.1270  ( 0.0512   )   0.0412

# Overall maximum convergence ratio:    0.2287
# Total of 2244 iteration steps.

# The results can also be viewed externally in the output file s50_3.txt
# It is advisable that you have a look at all three reports and
# understand how information is organized in each of them.

# To understand the table above, note that the "convergence t-ratio"
# is the t-ratio for convergence checking,
# not the t statistic for testing the significance of this effect!
# See Section 6.2 of the manual to understand this better.
# In the external output file, these are called
# "t-ratios for deviations from targets".
# The rule of thumb is that all t-ratios for convergence
# should ideally be less than 0.1 in absolute value,
# and the "Overall maximum convergence ratio" should be less than 0.25;
# this signifies good convergence of the algorithm.
# In the example here, this is the case.
# If this would not be the case, the best thing to do would be
# to continue the estimation, using the estimates produced here,
# and contained in ans, as the new initial values.
# This is explained below.
# Because the estimation algorithm is based on random simulations of the
# network evolution, there always will be small differences between
# different runs of the algorithm.
# To obtain "publication grade" estimates, where this variability is minimized,
# choose the value of parameter n3 in sienaAlgorithmCreate()
# ("Number of iterations in phase 3") larger than the default value of 1000;
# e.g., n3=5000.

# With function siena07 we made ans as the object containing
# all the results of the estimation. For example,

		ans$theta

# contains the vector of parameter estimates,

    ans$se

# contains the standard errors, while

		ans$covtheta

# contains the covariance matrix of the estimates.
# There are several "methods" available for viewing the object
# containing the results of the estimation.
# Above we already mentioned the commands
#		ans
# and
#		summary( ans )
# The command

	siena.table( ans )

# will produce in your working directory a table formatted
# for inclusion in a LaTeX document.
# The command

	siena.table( ans, type="html" )

# produces a table formatted in html,
# which can be included, e.g., in a Word document.
# See

		?print.sienaFit

# for further information, e.g., about the use of the xtable package
# for RSiena; if you use xtable, see the set of vignettes for xtable at
# http://cran.r-project.org/web/packages/xtable ,
# which gives more options.


############## MORE ON INITIALIZING PARAMETERS FOR ESTIMATION ########

# If the estimation algorithm has not produced good estimates
# (it 'has not converged well'),
# as will be indicated by some of the t-ratios for convergence being larger
# than 0.1 or the overall maximum convergence ratio being more than 0.25,
# (the precise values of these thresholds may be taken with a gain of salt)
# the best thing to do is continuing the estimation,
# using the estimates produced here,
# and contained in ans, as the new initial values.
# This is done by the option prevAns ('previous ans') as in

		ans <- siena07(myalgorithm, data=mydata, effects=myeff, prevAns=ans)

# the parameter estimates in ans then are extracted and
# used in the new estimation;
# moreover, Phase 1 will be omitted from the algorithm,
# as derivatives and covariance matrix are used from the previous run.
# This should be used only if the model specification in myeff
# has not changed, and if the provisional parameter estimates obtained
# in ans are reasonable; if they are not reasonable,
# make a fresh estimation without the prevAns parameter.

# In Section 6.3.1 of the manual you can read more about the initial
# values used for the estimation algorithm; but this rarely is of any concern.

# Sections 6.3-6.5 of the manual give further help.


################################################################################
###
### ---- Testing effects -------------------------------------------------------
###
################################################################################
#
# Three types of tests are available in SIENA.

# 1. t-type tests of single parameters can be carried out by dividing
# the parameter estimate by its standard error.
# Under the null hypothesis that the parameter is 0, these tests have
# approximately a standard normal distribution.
# 2. Score-type tests of single and multiple parameters are described
# in the manual.
# Parameters can be restricted by putting TRUE in the
# include, fix and test columns of the effects object.
# For example, to request a score test for the indegree popularity effect,
# the commands can be as follows.

		myeff <- setEffect(myeff, inPopSqrt, fix=TRUE, test=TRUE,
										  initialValue=0.0)
		ans <- siena07(myalgorithm, data=mydata, effects=myeff)

# After such an operation, again request

		summary(ans)

# to see the results, including those for the score test.
# You can also simply request

    score.Test(ans)

# see ?score.Test for further explanation.

# 3. Wald tests of single and multiple parameters can be obtained by means
# of the functions Wald.RSiena and Multipar.RSiena;
# see the help pages for these functions and also see the manual.


################################################################################
###
### ---- Time test -------------------------------------------------------------
###
################################################################################
#
# An application of the score test is given for the special case of parameter
# heterogeneity by Lospinoso et al. (2010) and implemented in RSiena.
# To apply the test to the results obtained above, request, e.g.,
		tt2 <- sienaTimeTest(ans)
		tt2
# If you wish more information, see
    summary(tt2)

################################################################################
###
### ---- Summary of model fitted -----------------------------------------------
###
################################################################################

# define original data
	friend.data.w1 <- s501
	friend.data.w2 <- s502
	friend.data.w3 <- s503
	drink <- s50a
	smoke <- s50s
# define RSiena data structures
	  friendship <- sienaDependent( array( c( friend.data.w1,
											  friend.data.w2, friend.data.w3 ),
									dim = c( 50, 50, 3 ) ) )
	  smoke1 <- coCovar( smoke[ , 1 ] )
	  alcohol <- varCovar( drink )

	  mydata <- sienaDataCreate( friendship, smoke1, alcohol )
# create effects structure
	myeff <- getEffects( mydata )
# get initial description
	  print01Report( mydata, modelname = 's50_3_init' )
# specify model
	  myeff <- includeEffects( myeff, transTrip, cycle3 )
	  myeff <- includeEffects( myeff, egoX, altX,
							   egoXaltX, interaction1 = "alcohol" )
	  myeff <- includeEffects( myeff, simX, interaction1 = "smoke1" )
# estimate
	  myAlgorithm <- sienaAlgorithmCreate( projname = 's50_3' )
	(ans <- siena07( myAlgorithm, data = mydata, effects = myeff))
# (the outer parentheses lead to printing the obtained result on the screen)
# if necessary, estimate further
    (ans <- siena07( myAlgorithm,
	                    data = mydata, effects = myeff, prevAns=ans))

################################################################################
###
### -- PROCEED TO Rscript04SienaBehaviour.R FOR
###									MODELING NETWORKS AND BEHAVIOUR BY RSIENA --
###
################################################################################
