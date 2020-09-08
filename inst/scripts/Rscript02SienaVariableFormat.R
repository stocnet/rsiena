###############################################################################
###                                                                         ###
### - Rscript02SienaVariableFormat.R: a script for the introduction to RSiena -
###                                                                         ###
###                               version September 8, 2020                 ###
###############################################################################
#
# The introductory script is divided into the following script files:
# Rscript01DataFormat.R, followed by
# RScriptSNADescriptives.R, code for descriptive analysis of the data,
# Rscript02SienaVariableFormat.R, which formats data and specifies the model,
# Rscript03SienaRunModel.R, which runs the model and estimates parameters
# Rscript04SienaBehaviour.R, which illustrates an example of analysing the
# coevolution of networks and behaviour
# Written by Tom Snijders, with earlier contributions from Robin Gauthier,
# Ruth Ripley, Johan Koskinen, and Paulina Preciado.
#
# This script, Rscript02SienaVariableFormat.R, sets up variables for analysis.
# The manipulations in this script requires that you have gone through the
# first part, "CALLING THE DATA AND PRELIMINARY MANIPULATIONS",
# of the script "Rscript01DataFormat.R" beforehand.

#### FORMATTING DATA ACCORDING TO THEIR ROLES AS VARIABLES IN A SIENA MODEL ####

# For this script, you will need the data read and modified in the script
# Rscript01DataFormat.R. If you have already ran that script, you may
# load the required workspace:
	load("WorkspaceRscript01.RData")
# If not, to make this script self-contained, you may run the commands:

	library(RSiena)
	friend.data.w1 <- s501
	friend.data.w2 <- s502
	friend.data.w3 <- s503
	drink <- s50a
	smoke <- s50s

# A number of objects need to be created in R, as preparations to letting
# siena07 execute the estimation. This will be indicated by
# A: dependent variables;
# B: explanatory variables;
# C: combination of dependent and explanatory variables;
# D: model specification.

# ---- A. ----------------------------------------------------------------------
# First we have to create objects for the dependent variables.

# sienaDependent creates a sienaDependent object, here a network,
# from a matrix or array or list of sparse matrix of triples.
# This object will have the role of a dependent variable in the analysis.
# The name of this network object (here: friendship) will be used
# in the output file.

 friendship <- sienaDependent(
                     array( c( friend.data.w1, friend.data.w2, friend.data.w3 ),
                     dim = c( 50, 50, 3 ) ) )

# The integers in the dim() here refer to the number of nodes (senders,
# receivers) and the number of waves.
# This object is an array of dimension 50 x 50 x 3, representing
# three adjacency matrices, with a number of attributes.
# You can ask for a brief decription of this object simply by typing

    print(friendship)

# for which the shorthand is

    friendship

# Note that this is an object of class

        class(friendship)

# with specific attributes and methods associated with it.
# You can get the detailed information by requesting

        dim( friendship )
        attributes( friendship )

# If you only are interested in the value of one particular attribute,
# you can request this by, e.g.,

    attr( friendship, "type")

# An extensive description of the friendship data is obtained by typing

   str(friendship)

# The function sienaDependent can also be used to create a behavior variable
# object with the extra argument type = "behavior".
# (Non-mentioned attributes get the default value, and in this case
# oneMode is the default; see below.)

# The 'drink' data (created in RscriptDataFormat.R ) is made available as a
# dependent behavior variable by the function

        drinkingbeh <- sienaDependent( drink, type = "behavior" )

# the class, class(drinkingbeh), is still sienaDependent.

# (Note: only use the variable in ONE role in a given model:
#  behavior variable or changing covariate!)

# The options available for defining a sienaDependent object
# are displayed by typing

        ?sienaDependent

# This shows that next to one-mode (unipartite)
# and behavior dependent variables,
# also two-mode (bipartite) dependent variables are possible.
# You can infer that oneMode is the default type from the fact
# that it is mentioned first.

# To create bipartite network objects you need two node sets
# and must create the node sets too.
# See

    ?sienaNodeSet

# where the Examples section shows an example of the syntax.


# ---- B. ----------------------------------------------------------------------
# Second we construct objects for the explanatory (independent) variables.
# From the help request
#       ?sienaDataCreate
# we see that these can be of five kinds:
	# coCovar            Constant actor covariates
	# varCovar           Time-varying actor covariates
	# coDyadCovar        Constant dyadic covariates
	# varDyadCovar       Time-varying dyadic covariates
	# compositionChange  Composition change indicators

# You can get help about this by the following requests:
#       ?coCovar
#       ?varCovar
#       ?coDyadCovar
#       ?varDyadCovar
#       ?sienaCompositionChange

# The variables available for this data set all are changing actor covariates.
# For illustrative purposes, we use smoking as observed at the first wave
# as a constant covariate:

        smoke1 <- coCovar( smoke[ , 1 ] )

# This selects the first column of smoke,
# which contains the first wave observations,
# and makes it available as a constant covariate.
# This is the pattern for for evey covariate file, e.g.
#       Attr1 <- coCovar( Covariate1 )
# where Covariate1 is a matrix with dim(Covariate1) equal to n x 1
# Note, if Covariates is a matrix with dim(Covariates) equal to n x p
# you can create constant covariates through
#       Attr1 <- coCovar(Covariates[,1])
#       ...
#       Attrp <- coCovar(Covariates[,p])

# We use the drinking data as a changing covariate.
# The function varCovar creates a changing covariate object from a matrix;
# the name comes from 'varying covariate'.

        alcohol <- varCovar( drink )

# You need at least three waves in the data set to define a varying covariate
# by the function varCovar; the reason is that the previous wave is used
# as a predictor of the next wave.

# The command

        attributes( alcohol )

# will tell you the information that RSiena now has added to the drink data.

# ---- C. ----------------------------------------------------------------------
# We now combine the dependent and independent variables.
# The function sienaDataCreate creates a Siena data object from input networks,
# covariates and composition change objects;
# the objects that earlier were created by sienaDependent will have the role
# of dependent variables, and similarly the other roles are predetermined
# by creation by the functions coCovar, varCovar,
# coDyadCovar, varDyadCovar, and sienaCompositionChange.

        mydata <- sienaDataCreate( friendship, smoke1, alcohol )

# You may check the result by requesting

		mydata

# You should now understand how this differs from the result of
#       mybehdata <- sienaDataCreate( friendship, smoke1, drinkingbeh )

# If you would like to use different names, you could request this as follows:
#        mydata <- sienaDataCreate( nominations = friendship, smoke1,
#                                   drinking = alcohol )

# Another type of dependent network is a two-mode network,
# which here is also called a bipartite network.
# To get advice for a first step in the use of two-mode networks,
# see the script RscriptSienaBipartite.R on the Siena scripts page.

# This finishes the data specification. Now we have to specify the model.
#

# ---- D. ----------------------------------------------------------------------
################################################################################
###
### ---- DEFINING EFFECTS ------------------------------------------------------
###
################################################################################
# The data set as combined in mydata implies a certain set of effects
# that can be included in the specification of the model.
# The function getEffects creates a dataframe of effects with a number of extra
# properties for use in RSiena:

        myeff <- getEffects( mydata )

# mydata is needed as an argument as the effects depend on the number
# and types of covariates and dependent variables.
# Before we explain the object myeff and how we shall be going to use it,
# we first produce a data description which is available now:

        print01Report( mydata, modelname = 's50_3_init' )

# This writes a basic report of the data to the file
# s50_3_init.txt in the current working directory. Locate and open it!
# Inspecting this is important because it serves as a check and also contains
# a number of basic descriptives.
# In this description you can see that the third wave for alcohol is not used.
# This is because changing covariates are assumed to be constant from one wave
# until immediately before the next wave,
# so that the values for the last wave are ignored.

# Let us now consider the myeff object, which is used to specify the model.
# It is of the class "sienaEffects", and contains the model specification.
# You can inspect the current model specification by simply requesting

       myeff

# For starting, the model specification is just a very limited default
# (including rates of change, outdegree and reciprocity only).
# To make a meaningful analysis, you will need to add to it.

# The rows of myeff correspond to the effects.
# By requesting

        names( myeff )

# you see the type of information that is stored about the effects,
# i.e., the columns (characteristics) defined for the effects.
# If desired, more information about these variables can be obtained
# from the help files:
#      ?getEffects
# where the characteristics are described.

# Some often used variables are effectName, shortName, type, and parameter.
# The set of available effects and their most used columns
# can be inspected as follows:

    effectsDocumentation(myeff)

# This creates an html file with the list of all effects available in myeff;
# the effects are all defined in Section 12 of the RSiena manual,
# but only those that are meaningful for dataset mydata are used in myeff.
# The "include" column defines whether effects are included in the model.

        myeff$include

# Here the TRUE values correspond to the default model specification which,
# however, is not meant as a serious model, being too limited.
# There are 3 main ways to operate on myeff.
# 1. Using RSiena functions "includeEffects", "setEffects", etc;
# 2. Changing myeff in spreadsheet form by the function fix();
# 3. Changing myeff directly by operating on its elements.
# The first way is most in line with the design philosophy of R,
# and allows you to save scripts that also can be used when
# there will be new versions of RSiena.
# Therefore, we suggest that for starting you only study
# option 1, ' Adding/removing effects using includeEffects'.
# The other two options are treated only for special
# and more difficult occasions; and for loooking 'behind the screens'.

# For identifying your effects you need the "shortName"s,
# which can be read in the manual
# (section "Mathematical definition of effects"),
# or obtained from the "effectsDocumentation()" function mentioned above.


# ---- 1. Adding/removing effects using includeEffects -------------------------
# The best way of specifying the model is by the includeEffects function.
# This function uses short names instead of full names.
# The short names are given by the effectsDocumentation() function
# mentioned above, and also are listed in the descriptions given in
# Section 12 of the manual.

# For illustration, let us start from scratch with a new sienaEffects object,
# and add the transitive triples and 3-cycles effects

        myeff <- getEffects( mydata )
        myeff <- includeEffects( myeff, transTrip, cycle3 )

# To see the current model specification,

        myeff

# Note that we can set several effects in one go!
# To remove an effect, e.g., the 3-cycle effects

        myeff <- includeEffects( myeff, cycle3, include=FALSE )

# Check again which effects now are included in the model

        myeff

# ---- Adding/removing covariate related effects -------------------------------
# The short names do not differentiate between the covariates:
# e.g., the effects 'alcohol ego' and 'smoke1 ego' both have shortName 'egoX',
# and the command

        myeff <- includeEffects( myeff, egoX )

# results in a message that does not (like the earlier one)
# confirm the newly included effect.
# The covariates are indicated by the variable "interaction1"
# in the sienaEffects object, listed as "inter1" in the result of
# effectsDocumentation(), and this has to be mentioned to include these effects:

        myeff <- includeEffects( myeff, egoX, altX, egoXaltX,
                                 interaction1 = "alcohol" )
        myeff <- includeEffects( myeff, simX, interaction1 = "smoke1" )

# We check the results again:

        myeff

# By looking at the help offered by

        ?includeEffects

# you can see how to include endowment and creation effects.
# Effects that depend on other variables, such as egoX, altX, etc. above,
# need the specification of these variables to define them.
# This is done by the interaction1 parameter
# when only one variable name is needed,
# and by interaction2 if there is a second variable involved,
# such as AltsAvAlt (see the manual).
# Although the names of these parameters are interaction1 and interaction2,
# this does not refer to the common concept of interaction
# in statistical modelling!

# ---- Creating interaction effects --------------------------------------------
# As a special topic, let us show how interaction effects are created.

# A convenient method to include an interaction is offered by the
# includeInteraction function.
# This can be used to interact two or three effects
# (if the interactions are allowed, which depends on their interactionType;
# see the manual for this).
# The interaction between smoke1 ego and reciprocity, for instance,
# can be defined by the command

        myeff <- includeInteraction( myeff, egoX, recip,
                                    interaction1 = c("smoke1","") )
		myeff

# This shows the interaction as an "unspecified interaction effect";
# but when printing results of the estimation the names of the
# interacting effects will be mentioned.
# E.g., an interaction between smoke1 ego and alcohol ego is defined by

#        myeff <- includeInteraction( myeff, egoX, egoX,
#                                    interaction1 = c( "smoke1", "alcohol" ) )

# Note that the keyword 'interaction1' used by RSiena is used for identifying
# the covariate for which the ego effect is selected, and does not
# refer to the interaction effect itself.
# If at least one of the interacting effects requires the interaction1 parameter
# for it specification, then this parameter is also required for the
# includeInteraction function.
# Then the two or three interaction1 parameters must be combined using c();
# the same goes for interaction2, if that also is necessary for the definition.
# As shown above for the recip effect, the interaction1 or interaction2 parameter
# is "" (i.e., an empty string) for effects where it is not needed.

# A second special topic is how to access other characteristics of effects.
# This can be done by the setEffect function.
# E.g., the dense triads effects
# counts the number of triplets with at least xx ties,
# where xx is the parameter of the effect, which can be 5 or 6
# (note that 6 is the maximum number of ties in a triplet).
# The default is 5. This is changed to 6 by the command

        myeff <- setEffect(myeff, denseTriads, parameter = 6)
		myeff

# The 'parameter' keyword refers to the effect parameter, described in
# Section 12 of the manual.

# ---- 2. Adding/removing effects using fix() ----------------------------------
# fix calls a data editor internal to R, so we can manually edit the effects.

#        fix( myeff )

# How to use fix() is presented merely for getting to know what myeff is.
# In practical analysis it is more convenient
# to use routine "includeEffects" instead, as explained above.
# fix() may not be usable if you do not have tcl/tk available!
# Note that the top of the dataframe shows the names of the columns:
# name, effectName, etc.
# You can edit the "include" column by changing the TRUE and FALSE values
# as required; when the editor is closed, the new values are stored.
# When you make an error, however, the effects object may be corrupted.
# Therefore, this way of adding and removing effects is more risky.

# ---- 3. Adding/removing effects by direct manipulation of myeff --------------
# Alternatively we can edit the dataframe directly by using R functions.
# You are advised to skip this part ("3.") at first and second reading,
# and read it only if you wish to get more understanding
# of the internal structure of the effects object.
# The commands below are used to set "include" to TRUE or FALSE,
# as an alternative to using the data editor.
# The "include" column with values TRUE or FALSE will always be
# located at the 9th column,
# but transitive triplets will not always be at the 16th row as this depends
# on the number of periods and variables; further, the list of available effects
# changes over different versions of RSiena.
# Some examples are the following
# (preceded by # because not proposed to be applied).

        #myeff[16,9] <- TRUE   #transitive triples
        #myeff[34,9] <- TRUE   #3 cycles
        #myeff[37,9] <- TRUE   #transitive ties
        #myeff[77,9] <- TRUE   #indegree popularity (sqrt)
        #myeff[86,9] <- TRUE   #outdegree popularity (sqrt)
        #myeff[93,9] <- TRUE   #indegree based activity (sqrt)
        #myeff[102,9] <- TRUE  #outdegree based activity (sqrt)
        #myeff[155,9] <- TRUE   #indegree-indegree assortativity
        #myeff[319,9] <- TRUE   #alcohol alter
        #myeff[331,9] <- TRUE   #alcohol ego
        #myeff[370,9] <- TRUE   #alcohol similarity
        #myeff[439,9] <- TRUE   #alcohol ego x alcohol alter

# But in other choices of data, the effect numbers will change.
# This is a reason why this is not a convenient method.

save.image("WorkspaceRscript02.RData")

################################################################################
###
### -- PROCEED TO Rscript03SienaRunModel.R FOR PARAMETER ESTIMATION BY RSIENA --
###
################################################################################
