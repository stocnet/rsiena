#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: transformScript.r
# *
# * Description: This module contains the code for the meta analysis of a
# * collection of Siena fits.
# *****************************************************************************/

##@transformScript transform R scripts to the new names
transformScript <- function(theoldscript, fileName="newScript.R",
					linelength=80, initialblanks=6,	logfile="log.txt",
					keepOldAlgorithms=FALSE)
{
# This function transforms theoldscript,
# supposed to be a script utilizing RSiena,
# to utilize the new function names as of RSiena version 1.6.2.
# Written by Tom Snijders, version February 24, 2026.

# Operation of TransformScript; see structure after "proper"
# 1. For functions of which nothing changes except the name,
#    all old function names are replaced by the new names.
#    These are in list 'namePairs' of name pairs.
# 2. For all other function names in 'furtherNames',
#    if it is "sienaAlgorithmCreate":
#             this command is executed, to have the algorithm object available.
#    for all others in 'furtherNames':
#             the function name is prepended by the single character 'z',
#             and this modified function, which is defined below, is executed.
#             This creates the new function call.
#    if it is "siena07" or "sienaBayes":
#             "zsiena07" creates not only the new call of "siena",
#             and "zsienaBayes" creates not only the new call of "multi_siena",
#             but also, if necessary, calls of algorithm creation functions.
#             For this purpose, it takes the algorithm object,
#             which was created by processing 
#             a preceding call of "sienaAlgorithmCreate",
#             and compares its values to the default values;
#             if any are different, then the appropriate is constructed
#             among set_model_saom, set_algorithm_saom, and set_output_saom,
#             or set_algorithm_mcm in the case of "sienaBayes",
#             and added to the new script 
#             before the call of "siena" or "multi_siena". 
#             This is done again for each new call of "siena07" or "sienaBayes",
#             which may result in superfluous calls of algorithm constructions.
#    For these function calls (i.e., in 'furtherNames'):
#             * Everything that comes after a hash (#) sign is disregarded;
#               if the function call itself comes after a hash sign,
#               it is not processed and the line is left unchanged.
#             * If the new function call is longer than linelength characters,
#               it is broken up in multiple lines,
#               after the position of commas,
#               and each line except the first starts with initialblanks blanks.
#             * For creation of a time-constant dyadic covariate ("coDyadCovar")
#               the keyword 'oneMode' is always used;
#               the user should check whether,
#               instead, it should be 'bipartite'.
# 3. Finally, for the remaining occurrences of names in 'furtherNames'
#             which will occur in lines after a "#" or without a "<-",
#             or for "Wald.RSiena", "testSame.RSiena", "score.Test" without a "("", 
#             the old function name is replaced by the new name.
# A log file is created indicating changes made
# and their line numbers in theoldscript.


################################################################################
### Data for alias substitution    #############################################
################################################################################

namePairs <- list(c("sienaDependent", "as_dependent_rsiena"),
   c("coCovar", "as_covariate_rsiena"),
   c("varDyadCovar", "as_covariate_rsiena"),
   c("sienaCompositionChange", "as_composition_rsiena"),
   c("sienaCompositionChangeFromFile", "as_composition_file_rsiena"),
   c("sienaNodeSet", "as_nodeset_rsiena"),
   c("sienaDataCreate", "make_data_rsiena"),
   c("sienaModelCreate", "sienaAlgorithmCreate"),
   c("sienaGroupCreate", "make_group_rsiena"),
   c("getEffects", "make_specification"),
   c("sienaDataConstraint", "make_constraint"),
   c("siena.table", "write_result"),
   c("sienaTimeTest", "test_time"),
   c("sienaGOF", "test_gof"),
   c("siena08", "meta_siena"),
   c("descriptives.sienaGOF", "descriptives"),
   c("selectionTable", "interpret_selection"),
   c("influenceTable", "interpret_influence"),
   c("sienaRIDynamics", "interpret_size_dynamics"),
   c("print01Report", "write_report"),
   c("modelname", "outputName"),
   c("updateTheta", "update_theta"),
   c("updateSpecification", "update_specification"),
   c("extract.posteriorMeans", "extract_posteriorMeans"),
	c("glueBayes", "multi_glue"),
	c("extract.posteriorMeans", "extract_posteriorMeans"),
	c("simpleBayesTest", "test_parameter"),
	c("multipleBayesTest", "test_parameter"),
	c("testedPar", "tested"),
	c("tested0", "tested_null")
	)

################################################################################
### Names of functions for which the change    #################################
### involves more than name replacement:       #################################
################################################################################

furtherNames <- list("sienaAlgorithmCreate", "siena07",
					"varCovar", "coDyadCovar", "includeEffects", "setEffect",
					"includeInteraction", "sienaRI", "Multipar.RSiena",
					"Wald.RSiena", "testSame.RSiena", "score.Test", 
					"sienaBayes", "simulateData")
found <- rep(0, length(namePairs) + length(furtherNames))
names(found) <- c(sapply(namePairs, function(x){x[[1]]}), furtherNames)

# With their new names and some others:
furtherNamePairs <- list(
		c("siena07", "siena"),
		c("varCovar", "as_covariate_rsiena"),
		c("coDyadCovar", "as_covariate_rsiena"),
		c("includeEffects", "set_effect"),
		c("setEffect", "set_effect"),
		c("includeInteraction", "set_interaction"),
		c("sienaRI", "interpret_size"),
		c("Multipar.RSiena", "test_parameter"),
		c("Wald.RSiena", "test_parameter"),
		c("testSame.RSiena", "test_parameter"),
		c("score.Test", "test_parameter"),
	c("sienaBayes", "multi_siena"),
	c("simulateData", "simulate"),
    c("extract.sienaBayes", "extract_sienaBayesFit"),
# perhaps the next ones are superfluous:
	c("interaction1", "covar1"),
	c("interaction2", "covar2"),
	c("interaction3", "covar3"))
furtherfound <- rep(0, length(furtherNamePairs))
names(furtherfound) <- c(sapply(furtherNamePairs, function(x){x[[1]]}))


################################################################################
### Helper functions    ########################################################
################################################################################

Reduce2 <- function(x, i=NULL){
# searches in string x for an occurrence of namePairs[[h]][[1]]
# and replaces it by namePairs[[h]][[2]]
# reports this as line i in file logfile and in array found
	out <- x
	for (h in 1:length(namePairs))
	{
		find <- grep(namePairs[[h]][[1]], out, useBytes=TRUE)
		if  (any(find))
		{
			write(paste("line ", i,'; found and replaced ', namePairs[[h]][[1]]),
							file=logfile, append=TRUE)
			out <- sub(namePairs[[h]][[1]], namePairs[[h]][[2]], out)
			found[h] <<- found[h]+1
		}
	}
	out
}


Reduce22 <- function(x, i=NULL){
# searches in string x for an occurrence of furtherNamePairs[[h]][[1]]
# and replaces it by furtherNamePairs[[h]][[2]]
# reports this as line i in file logfile and in array found
	out <- x
	for (h in 1:length(furtherNamePairs))
	{
		find <- grep(furtherNamePairs[[h]][1], out, useBytes=TRUE)
		if  (any(find))
		{
			write(paste("line ", i,'; found and replaced ', 
			            furtherNamePairs[[h]][[1]]), file=logfile, append=TRUE)
			out <- sub(furtherNamePairs[[h]][1], furtherNamePairs[[h]][2], out)
			furtherfound[h] <<- furtherfound[h]+1
		}
	}
	out
}

Reduce3 <- function(x, sep=", "){
# pastes a list of strings to one string
	out <- NULL
	if (length(x) == 0){return(NULL)}
	for (h in which(sapply(x, function(a){!is.null(a)})))
	{
		if (is.null(out))
		{
			out <- x[[h]]
		}
		else
		{
			out <- paste(out,x[[h]], sep=sep)
		}
	}
	return(out)
}

arr <- function(x,xx){
	paste(", ", x , "=", xx, sep="")
}

adtext <- function(s,t){
	return(paste(", ", s, "='", as.character(t), "'", sep=""))
}

pastent <- function(a,b){
	if (b)
	{
		return(paste(a, "=", b, sep=""))
	}
	else
	{
		return(NULL)
	}
}

noccur <- function(x, theline)
# number of occurrences of x in theline
{
	gr <- gregexpr(x, theline, fixed=TRUE, useBytes=TRUE)[[1]]
	ngr <- ifelse(gr >= 1, length(gr), max(gr,0))
	ngr
}

balanceParentheses <- function(x)
{
	nstarts <- noccur("(", x)[[1]]
	nends <- noccur(")", x)[[1]]
	while (nends != nstarts)
	{
		if (nends > nstarts)
		{
			x <- paste("(", x, sep="")
		}
		else
		{
			x <- paste(x, ")", sep="")
		}
		nstarts <- noccur("(", x)[[1]]
		nends <- noccur(")", x)[[1]]
	}
	x
}

textBeforeHash <- function(st, dropBlanks=FALSE)
# The text in st before "#" and perhaps without blanks
 {
	position <- regexpr("#", st, useBytes=TRUE)
	if (position >0 )
	{
		stnew <- substring(st, 1, position-1)
	}
	else
	{
		stnew <- st
	}
	if (dropBlanks)
	{
		stnew <- gsub(" ", "", stnew, fixed=TRUE, useBytes=TRUE)
	}
	stnew
}

shortentext <- function(stext, linewidth=linelength,
							startingblanks=initialblanks){
# transforms string stext to vector of strings of length <= linewidth,
# with breaks occurring after commas (if any),
# prepends startingblanks blanks at the start of follow-up lines
	if (linewidth <= startingblanks)
	{
		stop("startingblanks should be less than linewidth")
	}
	blanks <- Reduce(paste, rep(" ", startingblanks))
	if (is.null(stext)) return(NULL)
	stext_remain <- stext
	shorttext <- list(stext_remain)
	commas <- gregexpr(",", stext_remain, fixed=TRUE, useBytes=TRUE)[[1]]
	if (min(commas) > 0)
	{
		h <- 1
		textlength <- nchar(stext_remain)
		comma_pos <- ifelse((any(commas <= linewidth)),  max(commas[commas<=linewidth]),
					 min(commas[commas > linewidth]))
		while ((textlength > linewidth) && (comma_pos < textlength)){
			breakpos <- ifelse(textlength <= linewidth, textlength, comma_pos)
			shorttext[[h]] <- substr(stext_remain, 1, breakpos)
			shorttext[[h+1]] <- paste(blanks,
							substr(stext_remain, breakpos+1, textlength ))
			stext_remain <- shorttext[[h+1]]
			textlength <- nchar(stext_remain)
			commas <- gregexpr(",", stext_remain, fixed=TRUE, useBytes=TRUE)[[1]]
			if (min(commas) < 0)
			{
				comma_pos <- textlength + 1
			}
			else
			{
				comma_pos <- ifelse((any(commas <= linewidth)),  max(commas[commas<=linewidth]),
					 min(commas[commas > linewidth]))
			}
			h <- h+1
		}
	}
	unlist(shorttext)
}

findTheLines <- function(i, script){
# is the call distributed over several lines?
# Then these lines are concatenated and processed together.
# This considers only text before "#".
	linei <- textBeforeHash(script[i])
	nstarts <- noccur("(", linei)[[1]]
	nends <- noccur(")", linei)[[1]]
	newline <- linei
	ii <- i
	while (nends < nstarts)
	{
		ii <- ii+1
		if (ii > length(script))
		{
			stop(paste("No matching parentheses in what follows after line ", ii, sep=""))
		}
		lineii <- textBeforeHash(script[ii])
		newline <- paste(newline, lineii, collapse="")
		nstarts <- nstarts + noccur("(", lineii)[[1]]
		nends <- nends +  noccur(")", lineii)[[1]]
	}
	list(newline, ii)
}


################################################################################
### Function to process each function defined below that starts with z #########
################################################################################


transformCommand <- function(funName, oldscript, i, logfile)
{
# funName is a character string,
# oldscript is a string vector, newscript is a list
# transformCommand produces a string;
# but if funName=="siena07" or "sienaBayes", a vector of strings
# 	require(stringr) # for str_squish; not used now
	newline <- NULL
	oldline <- textBeforeHash(oldscript[i], dropBlanks=TRUE)
	if (funName %in% c("Multipar.RSiena", "Wald.RSiena", "testSame.RSiena", "score.Test"))
	{ # These may be used without an assignment to a name
		funNameWith <- paste(funName, "(", sep="")
	}
	else
	{
		funNameWith <- paste("<-",funName, "(", sep="")
	}
	if (any(grep(funNameWith, oldline, fixed=TRUE, useBytes=TRUE)))
	{
		found[funName] <<- found[funName]+1
# Define new function name, starting with z:
		newFunName <- paste("z", funName, sep="")
	# is the call distributed over several lines?
		TheLines <- findTheLines(i, oldscript)
		newline <-  Reduce(paste, TheLines[[1]])
#		newline <-  str_squish(Reduce(paste, TheLines[[1]]))
		ii <- TheLines[[2]]
		if (ii > length(oldscript))
		{
			newline <- paste("### The folliwng command is incomplete: ", newline, sep="")
		}
		else
		{
			theparts <- strsplit(newline, "<-", fixed=TRUE)
			if (is.list(theparts))
			{
				theparts <- theparts[[1]]
			}
			if (length(theparts)==1) # No occurrence of "<-" in newline
			{ 
				assignment <- FALSE
				newcommand <- newline
			}
			else
			{
				assignment <- TRUE
				newcommand <- theparts[2]
			}
			newcommand <- sub(funName, newFunName, newcommand)
			newcommand1 <- balanceParentheses(newcommand)
			if ((funName=="siena07") || (funName=="sienaBayes"))
# zsiena07 and zsienaBayes can produce a string vector with more than one element
			{
# Suppose that siena07 or sienaBayes is used in a function with a "..." argument:
				newcommand1 <- sub(", ...", "", newcommand1, fixed=TRUE)
				newcommand1 <- sub(",...", "", newcommand1, fixed=TRUE)
				thelines <- eval(parse(text=newcommand1))
				thelines$siena_text <- balanceParentheses(paste(theparts[1], "<-",
							thelines$siena_text, collapse=""))
				thelines$siena_text <- shortentext(thelines$siena_text)
# the other components of thelines were shortened in zsiena07 and zsienaBayes
				newline <- unlist(thelines)
			}
			else
			{
				newcommand <- unlist(eval(parse(text=newcommand1)))
				if (assignment)
				{
					newline <- paste(theparts[1], "<- ", newcommand, sep="")
				}
				else
				{
					newline <- newcommand
				}
				newline1 <- balanceParentheses (newline)
				newline <- shortentext(newline1)
			}
			write(paste("\nline ", i, "to", ii, '; found and transformed', funName),
							file=logfile, append=TRUE)
			i <- ii
			if (funName=="simulateData")
			{
				write("NOTE SCRIPT For simulate, nsim=100 is arbitrary; please modify." ,
								file=logfile, append=TRUE)
				cat("NOTE: For simulate (replacing simulateData), nsim=100 is arbitrary; please modify.\n" )
			}
		}
		write("lines \r", file=logfile, append=TRUE)
	}
	list(newline=newline, i=i)
}

################################################################################
### Functions to process the algorithm given in siena07 or sienaBayes   ########
################################################################################

modelalgorithm <- function(algo){
#  model algorithm
	alg0 <- list()
	if (!is.null(algo$modelType))
	{
		textmt <- paste(names(algo$modelType),
			algo$modelType, sep = "=", collapse = ",")
		alg0[[1]] <- paste("modelType=c(", textmt, ")", sep="")
	}
	if (!is.null(algo$behModelType))
	{
		textmt <- paste(names(algo$behModelType),
			algo$behModelType, sep = "=", collapse = ",")
		alg0[[2]] <- paste("behModelType=c(", textmt, ")", sep="")
	}
	if (is.null(algo$MaxDegree))
	{
		alg0[[3]] <- NULL
	}
	else
	{
		textmd <- paste(names(algo$MaxDegree),
				algo$MaxDegree, sep = "=", collapse = ",")
		alg0[[3]] <- paste("MaxDegree=c(", textmd, ")", sep="")
	}

	if (!is.null(algo$UniversalOffset))
	{
		textoff <- paste(names(algo$UniversalOffset),
			algo$UniversalOffset, sep = "=", collapse = ",")
		alg0 <- c(alg0, paste("OffSet=c(", textoff, ")", sep=""))
	}
	mod_content <- Reduce3(alg0)
	if (is.null(mod_content))
	{
		alg_model_text <- NULL
	}
	else
	{
		alg_model_text <- paste("alg_model <- set_model_saom(",
									mod_content, ")", sep="")
		alg_model_text <- shortentext(alg_model_text)
	}
	alg_model_text
}

algalgorithm <- function(algo){
	ndtext <- function(a) # not default
	{
		a_value <- with(algo, eval(parse(text=a)))
		defalg <- sienaAlgorithmCreate(silent=TRUE)
		def_value <- with(defalg, eval(parse(text=a)))
		if ((as.character(a_value) != as.character(def_value)) &
		      (!(is.null(a_value) & (def_value == "NULL")))) # comparison with default value
		{
			if (is.character(a_value))
			{
#				a_value <- paste("'", a_value, "'", sep="")
			}
			return(paste(a,"=",a_value, sep=""))
		}
		else
		{
			return(NULL)
		}
	}
# Begin algorithm algorithm
	alg0 <- list()
  # Sets SAOM estimation procedure:
	alg0 <- c(alg0, ndtext("maxlike"))
	alg0 <- c(alg0, ndtext("gmm"))
	if (!is.na(algo$cconditional))
	{
		alg0 <- c(alg0, paste("cond=", algo$cconditional, sep=""))
	}
	alg0 <- c(alg0, ndtext("condvarno"))
	alg0 <- c(alg0, ndtext("condname"))
	alg0 <- c(alg0, ndtext("simOnly"))
	if (!is.null(algo$randomSeed))
	{
		if (algo$randomSeed != "NULL")
		{
			alg0 <- c(alg0, paste("seed=",algo$randomSeed, "", sep=""))
		}
	}
  # Specification of Robbins Monro algorithm:
	alg0 <- c(alg0, ndtext("n3"))
	alg0 <- c(alg0, ndtext("nsub"))
	alg0 <- c(alg0, ndtext("n2start"))
#	if (!is.null(algo$n2start)) #  ndtext("n2start")
#	{
#		alg0 <- c(alg0, paste("n2start=", algo$n2start), sep="")
#	}
	alg0 <- c(alg0, ndtext("firstg"))
	alg0 <- c(alg0, ndtext("reduceg"))
	alg0 <- c(alg0, ndtext("truncation"))
	alg0 <- c(alg0, ndtext("doubleAveraging"))
	alg0 <- c(alg0, ndtext("diagonalize"))
	alg0 <- c(alg0, ndtext("standardizeVar"))
	alg0 <- c(alg0, ndtext("dolby"))
	alg0 <- c(alg0, ndtext("useStdInits"))
	alg0 <- c(alg0, pastent("findiff", algo$FinDiff.method))
  # Specification of likelihood algorithm:
	alg0 <- c(alg0, ndtext("mult"))
#	if ((!is.null(algo$mult)))
#	{
#		alg0 <- c(alg0, paste("mult=", algo$mult, sep=""))
#	}
	alg0 <- c(alg0, ndtext("prML"))
	alg0 <- c(alg0, ndtext("maximumPermutationLength"))
	alg0 <- c(alg0, ndtext("minimumPermutationLength"))
	alg0 <- c(alg0, ndtext("initialPermutationLength"))
	alg0 <- c(alg0, ndtext("localML"))
#browser()
	alg0
}

################################################################################
### Functions for handling function names in furtherNames   ####################
################################################################################


zsiena07 <- function(alg, data=NULL, effects=NULL, batch=FALSE,
		verbose=FALSE, silent=FALSE,
        useCluster=FALSE, nbrNodes=2,
        thetaValues = NULL,
        returnThetas = FALSE,
        thetaBound = 50,
        targets = NULL,
        initC=TRUE, clusterString=rep("localhost", nbrNodes), tt=NULL,
        parallelTesting=FALSE, clusterIter=!alg$maxlike,
        clusterType=c("PSOCK", "FORK"), cl=NULL,
		prevAns=NULL, returnDeps=FALSE,
		returnChains=FALSE, returnDataFrame=FALSE,
		returnChangeContributions=FALSE, ...)
{
# default algorithm:
	defalg <- sienaAlgorithmCreate(silent=TRUE)
	algname <- deparse(substitute(alg))	
	available <- suppressWarnings(!inherits(try(alg, silent=TRUE),
									"try-error"))
	if (available)
	{
		available <- inherits(alg, "sienaAlgorithm")
	}
	if (!available) 
	{ # not available from a previous run of zsienaAlgorithmCreate in this script
		alg <- sienaAlgorithmCreate(silent=TRUE)
		cat("The algorithm for siena07 is not available in this script.\n")
		cat("The name given in the script is ",algname ,".\n")
		cat("If", algname, "is not a default algorithm object,\n ")
		cat("you either have to supply it in the old script, or modify the new script.\n")
		write(paste("\nAlgorithm ", algname, "of siena07 not found.", sep=" "),
						file=logfile, append=TRUE)
		commentline <- paste("### Algorithm ", algname,
							"is not available in this script. ###", sep="")
		alg_model_text <- NULL
		alg_out_text <- NULL
		alg_alg_text <- NULL
	}
	else # alg does exist
	{
# Define output algorithm
		commentline <- NULL
		alg_out_text <- NULL
		algo <- alg # The algorithm given to siena07
		alg0 <- list()
		alg0[[1]] <- pastent("lessMem", (!algo$sf2.byIteration))
		alg0[[2]] <- pastent("returnThetas", returnThetas)
		alg0[[3]] <- pastent("returnChains", returnChains)
		alg0[[4]] <- pastent("returnDataFrame", returnDataFrame)
		alg0[[5]] <- pastent("returnChangeContributions", returnChangeContributions)
		out_content <- Reduce3(alg0)
		if (is.null(out_content))
		{
			alg_out_text <- NULL
		}
		else
		{
			alg_out_text <- paste("alg_out <- set_output_saom(", out_content, ")", sep="")
			alg_out_text <- shortentext(alg_out_text)
		}
# Define model algorithm
		alg_model_text <- modelalgorithm(algo)
# Define algorithm algorithm
		alg0 <- algalgorithm(algo)
		if (!is.null(targets))
		{
			alg0 <- c(alg0, paste("targets=",
						deparse1(substitute(targets)), "", sep=""))
		}
		if (!is.null(thetaValues))
		{
			alg0 <- c(alg0, paste("thetaValues=",
						deparse1(substitute(thetaValues)), "", sep=""))
		}
		alg_content <- Reduce3(alg0)
		if (is.null(alg_content))
		{
			alg_alg_text <- NULL
		}
		else
		{
			alg_alg_text <- paste("alg_alg <- set_algorithm_saom(",
									alg_content, ")", sep="")
			alg_alg_text <- shortentext(alg_alg_text)
		}
	}
## Text of call of siena:
	coretext <- paste("siena(data=", deparse1(substitute(data)),
						", effects=", deparse1(substitute(effects)), sep="")
	text0 <- list(coretext)
	if (hasArg("thetaBound"))
	{
		text0 <- c(text0, arr("thetaBound", substitute(thetaBound)))
	}
	if (hasArg("returnDeps"))
	{
		text0 <- c(text0, arr("returnDeps", substitute(returnDeps)))
	}
	if (hasArg("batch"))
	{
		text0 <- c(text0, arr("batch", substitute(batch)))
	}
	if (hasArg("silent"))
	{
		text0 <- c(text0, arr("silent" ,substitute(silent)))
	}
	if (hasArg("verbose"))
	{
		text0 <- c(text0, arr("verbose", substitute(verbose)))
	}
	if (hasArg("initC"))
	{
		text0 <- c(text0, arr("initC", substitute(initC)))
	}
	if (useCluster & (nbrNodes >= 2))
	{
		text0 <- c(text0, arr("nbrNodes", substitute(nbrNodes)))
	}
	if (hasArg("clusterString"))
	{
		text0 <- c(text0, arr("clusterString", substitute(clusterString)))
	}
	if (hasArg("clusterType"))
	{
		text0 <- c(text0, adtext("clusterType", clusterType))
	}
	if (hasArg("cl"))
	{
		text0 <- c(text0, arr(cl, substitute()))
	}
	if (hasArg(prevAns))
	{
		text0 <- c(text0, paste(", prevAns=",
					deparse1(substitute(prevAns)), sep=""))
	}
	text0 <- Reduce3(text0, sep="")
	textc <- list()
	anyalg <- FALSE
	if (!is.null(alg_model_text))
	{
		textc <- c(textc, "control_model=alg_model")
		anyalg <- TRUE
	}
	if (!is.null(alg_alg_text))
	{
		textc <- c(textc, "control_algo=alg_alg")
		anyalg <- TRUE
	}
	if (!is.null(alg_out_text))
	{
		textc <- c(textc, "control_out=alg_out")
		anyalg <- TRUE
	}
	if (anyalg)
	{
		siena_text <- paste(text0, ", ", Reduce3(textc), ")", sep="")
	}
	else
	{
		siena_text <- paste(text0, ")", sep="")
	}
# siena_text will be shortened in transformCommand,
# because something will be prepended
	return(list(commentline=commentline, alg_out_text=alg_out_text,
			alg_model_text=alg_model_text, alg_alg_text=alg_alg_text,
			siena_text=siena_text))
}


zvarCovar<- function(val, centered=TRUE, nodeSet="Actors", warn=TRUE,
						imputationValues=NULL)
{
	text0 <- paste("as_covariate_rsiena(", deparse1(substitute(val)), sep="")
	text0 <- paste(text0, "type='monadic'", sep=", ")
	if (!centered)
	{
		text0 <- paste(text0, "centered=FALSE", sep=", ")
	}
	if (hasArg(nodeSet))
	{
		text0 <- paste(text0, ", nodeSet=",
					deparse1(substitute(nodeSet)), sep="")
	}
	if (!warn)
	{
		text0 <- paste(text0, "warn=FALSE", sep=", ")
	}
	if (hasArg(imputationValues))
	{
		text0 <- paste(text0, ", imputationValues=",
					deparse1(substitute(imputationValues)), sep="")
	}
	cov_text <- paste(text0, ")", sep="")
	cov_text
}

zcoDyadCovar<- function(val, centered=TRUE, nodeSets=c("Actors","Actors"),
					   warn=TRUE, sparse=inherits(val,"TsparseMatrix"),
					   type=c("oneMode", "bipartite"))
{
	text0 <- paste("as_covariate_rsienas(", deparse1(substitute(val)), sep="")
	text0 <- paste(text0, "type='oneMode'", sep=", ")
	if (!centered)
	{
		text0 <- paste(text0, "centered=FALSE", sep=", ")
	}
	if (hasArg(nodeSets))
	{
		text0 <- paste(text0, ", nodeSet=",
					deparse1(substitute(nodeSets)), sep="")
	}
	if (!warn)
	{
		text0 <- paste(text0, "warn=FALSE", sep=", ")
	}
	cov_text <- paste(text0, ")", sep="")
	cov_text
}

zincludeEffects <- function(myeff, ..., include=TRUE, name=myeff$name[1],
						   type="eval", interaction1="", interaction2="",
						   fix=FALSE, test=FALSE,
						   character=FALSE, verbose=TRUE)
{
	text0 <- paste("set_effect(", deparse1(substitute(myeff)), sep="")
	dots <- substitute(list(...))
	if (length(dots)==2)
	{
		dots <- as.character(dots)[2]
	}
	text0 <- toString(c(text0, dots))
	if (hasArg(type))
	{
		text0 <- paste(text0, ", type=",
					deparse1(substitute(type)), sep="")
	}
	if (hasArg(name))
	{
		text0 <- paste(text0, ", depvar=",
					deparse1(substitute(name)), sep="")
	}
	if (hasArg(interaction1))
	{
		text0 <- paste(text0, ", covar1=",
					deparse1(substitute(interaction1)), sep="")
	}
	if (hasArg(interaction2))
	{
		text0 <- paste(text0, ", covar2=",
					deparse1(substitute(interaction2)), sep="")
	}
	if (!include)
	{
		text0 <- paste(text0, "include=FALSE", sep=", ")
	}
	if (fix)
	{
		text0 <- paste(text0, "fix=TRUE", sep=", ")
	}
	if (test)
	{
		text0 <- paste(text0, "test=TRUE", sep=", ")
	}
	if (!verbose)
	{
		text0 <- paste(text0, "verbose=FALSE", sep=", ")
	}
	setef_text <- paste(text0, ")", sep="")
	setef_text
}

zsetEffect <- function(x, shortName, parameter = NULL,
    fix = FALSE, test = FALSE, random=FALSE, initialValue = 0,
    timeDummy = ",", include = TRUE,
    name = x$name[1], type = "eval", interaction1 = "",
    interaction2 = "", effect1=0, effect2=0, effect3=0,
    period=1, group=1, character=FALSE, verbose = TRUE)
{
	text0 <- paste("set_effect(", deparse1(substitute(x)), sep="")
	text0 <- paste(text0,
					deparse1(substitute(shortName)), sep=", ")
	if (hasArg(type))
	{
		text0 <- paste(text0, ", type=",
					deparse1(substitute(type)), sep="")
	}
	if (hasArg(name))
	{
		text0 <- paste(text0, ", depvar=",
					deparse1(substitute(name)), sep="")
	}
	if (hasArg(interaction1))
	{
		text0 <- paste(text0, ", covar1=",
					deparse1(substitute(interaction1)), sep="")
	}
	if (hasArg(interaction2))
	{
		text0 <- paste(text0, ", covar2=",
					deparse1(substitute(interaction2)), sep="")
	}
	if (hasArg(effect1))
	{
		text0 <- paste(text0, ", effect1=", effect1, sep="")
	}
	if (hasArg(effect2))
	{
		text0 <- paste(text0, ", effect2=", effect2, sep="")
	}
	if (hasArg(effect3))
	{
		text0 <- paste(text0, ", effect3=", effect3, sep="")
	}
	if (hasArg(period))
	{
		text0 <- paste(text0, ", period=", period, sep="")
	}
	if (hasArg(group))
	{
		text0 <- paste(text0, ", group=", group, sep="")
	}
	if (hasArg(parameter))
	{
		text0 <- paste(text0, ", parameter=", parameter, sep="")
	}
	if (hasArg(initialValue))
	{
		text0 <- paste(text0, ", initialValue=",
						deparse1(substitute(initialValue)), sep="")
	}
	if (hasArg(random))
	{
		text0 <- paste(text0, ", random=", random, sep="")
	}
	if (hasArg(timeDummy))
	{
		text0 <- paste(text0, ", timeDummy=",  deparse1(substitute(timeDummy)), sep="")
	}
	if (!include)
	{
		text0 <- paste(text0, "include=FALSE", sep=", ")
	}
	if (fix)
	{
		text0 <- paste(text0, "fix=TRUE", sep=", ")
	}
	if (test)
	{
		text0 <- paste(text0, "test=TRUE", sep=", ")
	}
	if (!verbose)
	{
		text0 <- paste(text0, "verbose=FALSE", sep=", ")
	}
	sete_text <- paste(text0, ")", sep="")
	sete_text
}

zincludeInteraction <- function(x, ..., include = TRUE, name = x$name[1],
    type = "eval", interaction1 = rep("", 3), interaction2 = rep("", 3),
    fix=FALSE, test=FALSE, random=FALSE,
    initialValue=0,
    character = FALSE, verbose = TRUE)
{
	text0 <- paste("set_interaction(", deparse1(substitute(x)), sep="")
	dots <- substitute(list(...))
	text0 <- toString(c(text0, dots))
	if (hasArg(type))
	{
		text0 <- paste(text0, ", type=",
					deparse1(substitute(type)), sep="")
	}
	if (hasArg(name))
	{
		text0 <- paste(text0, ", depvar=",
					deparse1(substitute(name)), sep="")
	}
	if (hasArg(interaction1))
	{
		text0 <- paste(text0, ", covar1=",
					deparse1(substitute(interaction1)), sep="")
	}
	if (hasArg(interaction2))
	{
		text0 <- paste(text0, ", covar2=",
					deparse1(substitute(interaction2)), sep="")
	}
	if (hasArg(initialValue))
	{
		text0 <- paste(text0, ", initialValue=",
					deparse1(substitute(initialValue)), sep="")
	}
	if (hasArg(random))
	{
		text0 <- paste(text0, ", random=", random, sep="")
	}
	if (!include)
	{
		text0 <- paste(text0, "include=FALSE", sep=", ")
	}
	if (fix)
	{
		text0 <- paste(text0, "fix=TRUE", sep=", ")
	}
	if (test)
	{
		text0 <- paste(text0, "test=TRUE", sep=", ")
	}
	if (!verbose)
	{
		text0 <- paste(text0, "verbose=FALSE", sep=", ")
	}
	seti_text <- paste(text0, ")", sep="")
	seti_text
}

zsienaRI <- function(data, ans=NULL, theta=NULL, effects=NULL,
	getChangeStats=FALSE)
{
	if (!is.null(effects))
	{
		text0 <- paste("interpret_size(",
						deparse1(substitute(effects)), sep="")
		if (hasArg(theta))
		{
			text0 <- paste(text0, ", theta=",
					deparse1(substitute(theta)), sep="")
		}
	}
	else if (hasArg(ans))
# or replace by (deparse1(substitute(ans)) != "NULL") ?
	{
		text0 <- paste("interpret_size(",
						deparse1(substitute(ans)), sep="")
	}
	else
	{
		stop("sienaRI called without effects and ans")
	}
	text0 <- paste(text0, ", data=",
					deparse1(substitute(data)), sep="")
	if (getChangeStats)
	{
		text0 <- paste(text0, "getChangeStats=TRUE", sep=", ")
	}
	setin_text <- paste(text0, ")", sep="")
	setin_text
}


zWald.RSiena <- function(A, x)
{
	paste("test_parameter(", deparse1(substitute(x)), ", tested=", 
				deparse1(substitute(A)), ")", sep="")
}


zMultipar.RSiena <- function(x, tested)
{	
	text0 <- paste("test_parameter(", deparse1(substitute(x)), ", tested=", 
				deparse1(substitute(tested)), ")", sep="")
	text0
}

zscore.Test <- function(x, tested=NULL)
{
	if (is.null(tested))
	{
	text0 <- paste("test_parameter(", deparse1(substitute(x)), ", method='score')", 
						sep="")
	}
	else
	{
	text0 <- paste("test_parameter(", deparse1(substitute(x)), ", method='score', tested=", 
				deparse1(substitute(tested)), ")", sep="")
	}
	text0
}

ztestSame.RSiena <- function(x, e1, e2)
{
	paste("test_parameter(", deparse1(substitute(x)), ", method='same', tested=", 
				deparse1(substitute(e1)), ", tested2=", 
				deparse1(substitute(e2)), ")", sep="")
}

zset_algorithm_mcmc <- function(nprewarm=50, nwarm=50, nmain=250, nrunMHBatches=20, nImproveMH=100,
    targetMHProb=0.25,
    nSampVarying=1, nSampConst=1, nSampRates=0,
    priorPrecFactor=1,
    initgainGlobal=0.1, initgainGroupwise = 0.001, delta=1e-10,
    initML=FALSE,
    incidentalBasicRates=FALSE,
	initfgain=0.2, gamma=0.05)
# This is set_algorithm_mcmc of multiSiena version 1.2.36,
# without class and version statements.
# Sets details of mcmc algorithm specification
{
	algo <- list()
	algo$nprewarm <- nprewarm
	algo$nwarm <- nwarm
	algo$nmain <- nmain
	algo$nrunMHBatches <- nrunMHBatches
	algo$nImproveMH <- nImproveMH
	algo$targetMHProb <- targetMHProb
	algo$nSampVarying <- nSampVarying
	algo$nSampConst <- nSampConst
	algo$nSampRates <- nSampRates
	algo$priorPrecFactor <- priorPrecFactor
	algo$initgainGlobal <- initgainGlobal
	algo$initgainGroupwise <- initgainGroupwise
	algo$initML <- initML
	algo$incidentalBasicRates <- incidentalBasicRates
	algo$initfgain <- initfgain
	algo$gamma <- gamma
	algo$delta <- delta
	algo
}


zsienaBayes <- function(data, effects, algo, saveFreq=100,
				initgainGlobal=NULL, initgainGroupwise =NULL, initfgain=NULL,
				gamma=NULL, initML=NULL,
				priorMeanEta=NULL, priorSigEta=NULL,
				priorMu=NULL, priorSigma=NULL, priorDf=NULL, priorKappa=0,
				priorRatesFromData=2,
				priorPrecFactor=NULL,
				frequentist=FALSE, incidentalBasicRates=NULL,
				reductionFactor=0.5,
				delta=NULL,
				nprewarm=NULL, nwarm=NULL, nmain=NULL, nrunMHBatches=NULL,
				nSampVarying=NULL, nSampConst=NULL, nSampRates=NULL,
				nImproveMH=NULL, targetMHProb=NULL,
				lengthPhase1=round(nmain/5), lengthPhase3=round(nmain/5),
				storeScores = FALSE,
				storeAll = FALSE, prevAns=NULL, usePrevOnly=TRUE,
				prevBayes = NULL, newProposalFromPrev=(prevBayes$nwarm >= 1),
				silentstart=TRUE,
				nbrNodes=1, clusterType=c("PSOCK", "FORK"),
				getDocumentation=FALSE)
{
# Construct mcmc algorithm
# In the call of zsienaBayes, I changed all entries for set_algorithm_mcmc to NULL
# alg mcmc text
	setmcmc <- zset_algorithm_mcmc()
	alg_text <- ""
	i <- 0
# Go through the arguments of sienaBayes that were transferred to set_algorithm_mcmc:
	for (nam in names(setmcmc))
	{
		if (!is.null(eval(parse(text=nam))))
		{
			i <- i+1
			if (i <= 1)
			{
			alg_text <- paste(nam, "=", eval(parse(text=nam)), sep=" ")
			}
			else
			{
			alg_text <- paste(alg_text, ", ", nam, " = ", eval(parse(text=nam)), sep="")
			}
		}
	}
	if (alg_text == "")
	{
		alg_text <- NULL
	}
	else
	{
		alg_text <- paste("alg_mcmc <- set_algorithm_mcmc(",
									alg_text, ")", sep="")
		alg_text <- shortentext(alg_text)
	}

	algname <- deparse(substitute(algo))
	if (!exists(algname)) # Is it available?
	{ # not available from a previous run of zsienaAlgorithmCreate in this script
		alg <- sienaAlgorithmCreate(silent=TRUE)
		cat("The algorithm for sienaBayes is not available in this script.\n")
		cat("The name given in the script is ",algname ,".\n")
		cat("If", algname, "is not a default algorithm object,\n ")
		cat("you either have to supply it in the old script, or modify the new script.\n")
		write(paste("\nAlgorithm ", algname, "of sienaBayes not found.", sep=" "),
						file=logfile, append=TRUE)
		commentline <- paste("### Algorithm ", algname,
							"is not available in this script. ###", sep="")
		alg_model_text <- NULL
		alg_out_text <- NULL
		alg_alg_text <- NULL
	}
	else
	{
		commentline <- NULL
# Define model algorithm
		alg_model_text <- modelalgorithm(algo)
# Define algorithm algorithm
		alg_content <- Reduce3(algalgorithm(algo))
		if (is.null(alg_content))
		{
			alg_alg_text <- NULL
		}
		else
		{
			alg_alg_text <- paste("alg_alg <- set_algorithm_saom(",
									alg_content, ")", sep="")
			alg_alg_text <- shortentext(alg_alg_text)
		}
	}
## Define text of call of multi_siena:
	coretext <- paste("multi_siena(data=", deparse1(substitute(data)),
						", effects=", deparse1(substitute(effects)), sep="")
	text0 <- list()
	if (hasArg("priorMu"))
	{
         text0 <- c(text0, arr("priorMu" ,substitute(priorMu)))
	}

	if (hasArg("priorSigma"))
	{
        text0 <- c(text0, arr("priorSigma" ,substitute(priorSigma)))
	}
	if (hasArg("priorDf"))
	{
		text0 <- c(text0, arr("priorDf" ,substitute(priorDf)))
	}
	if (hasArg("priorKappa"))
	{
		text0 <- c(text0, arr("priorKappa" ,substitute(priorKappa)))
	}
	if (hasArg("priorMeanEta"))
	{
        text0 <- c(text0, arr("priorMeanEta" ,substitute(priorMeanEta)))
	}
	if (hasArg("priorSigEta"))
	{
        text0 <- c(text0, arr("priorSigEta" ,substitute(priorSigEta)))
	}
	if (hasArg("priorRatesFromData"))
	{
		text0 <- c(text0, arr("priorRatesFromData", substitute(priorRatesFromData)))
	}
	if (hasArg("reductionFactor"))
	{
        text0 <- c(text0, arr("reductionFactor" ,substitute(reductionFactor)))
	}
	if (hasArg("saveFreq"))
	{
		text0 <- c(text0, arr("saveFreq" , substitute(saveFreq)))
	}
	if (hasArg("storeScores"))
	{
		text0 <- c(text0, arr("storeScores", substitute(storeScores)))
	}
	if (hasArg(prevAns))
	{
		text0 <- c(text0, paste(", prevAns=",
					deparse1(substitute(prevAns)), sep=""))
	}
	if (hasArg("usePrevOnly"))
	{
		text0 <- c(text0, arr("usePrevOnly", usePrevOnly))
	}
	if (hasArg("newProposalFromPrev"))
	{
		text0 <- c(text0, arr("newProposalFromPrev", substitute(newProposalFromPrev)))
	}
	if (hasArg("silentstart"))
	{
		text0 <- c(text0, arr("silentstart" , substitute(silentstart)))
	}
	if (hasArg("nbrNodes"))
	{
		text0 <- c(text0, arr("nbrNodes", substitute(nbrNodes)))
	}
	if (hasArg("clusterType"))
	{
		text0 <- c(text0, adtext("clusterType", clusterType))
	}
	text0 <- Reduce3(text0, sep="")

	textc <- list()
	anyalg <- FALSE
	if (!is.null(alg_model_text))
	{
		textc <- c(textc, "control_model=alg_model")
		anyalg <- TRUE
	}
	if (!is.null(alg_alg_text))
	{
		textc <- c(textc, "control_algo=alg_alg")
		anyalg <- TRUE
	}
	if (!is.null(alg_text))
	{
		textc <- c(textc, "control_mcmc=alg_mcmc")
		anyalg <- TRUE
	}
	if (anyalg)
	{
		text1 <- paste(coretext, Reduce3(textc), sep=", ")
		siena_text <- paste(text1, text0, ")", sep="")
	}
	else
	{
		siena_text <- paste(coretext, text0, ")", sep="")
	}
# siena_text will be shortened in transformCommand,
# because something will be prepended
	return(list(commentline=commentline, alg_text=alg_text,
			alg_model_text=alg_model_text, alg_alg_text=alg_alg_text,
			siena_text=siena_text))
}

zsimulateData <- function(g, x, xd, thin=1, nstart=x$nwarm,
                            seed=NULL, nbrNodes=1)
{
	text0 <- paste("simulate(", substitute(x), sep="")
#	nm <- sum(!is.na(x$ThinParameters[,1,1])) - x$nwarm
# but x may be unavailable
#   nsim <- nm %/% thin # integer division
	nsim <- 100  # arbitrary
	text0 <- paste(text0, ", nsim=", nsim, sep="")
	if (hasArg(seed))
	{
		text0 <- paste(text0, ", seed=", seed, sep="")
	}
	text0 <- paste(text0, ", data=",
					deparse1(substitute(xd)), sep="")
	if (hasArg(nbrNodes))
	{
		text0 <- paste(text0, ", nbrNodes=", nbrNodes, sep="")
	}
	paste(text0, ")", sep="")
}


################################################################################
### Function transformScript proper   ##########################################
################################################################################

# To enforce "silent=TRUE" in calls of sienaAlgorithmCreate:
FromTransformScript <- TRUE

cat("Processing ", deparse1(substitute(theoldscript)), ".\n")
# write beginning of log file:
write(c(paste("This is file  ", deparse1(substitute(logfile)), sep=""),
		"log file for function TransformScript ",
		paste("transforming '", deparse1(substitute(theoldscript)),
					"' to the new names of version 1.6.", sep=""),
		date(), ""),
		file = logfile)

thenewscript <- list()

write("PASS 1", append=TRUE, file=logfile)
# Replace function names by new names ('alias')
# for those functions which further are unaltered:
i <- 0
repeat{
	i <- i+1
	if (i > length(theoldscript)) break
	thenewscript[[i]] <- Reduce2(theoldscript[i], i )
}

thenewscript <- unlist(thenewscript)
# This makes it into a vector of strings, each line a string.
# This is necessary because commands might be spread over several lines,
# and findTheLines, used below, uses string vectors.


write(c("", "PASS 2"), append=TRUE, file=logfile)
# Process functions in furtherNames.

thenewerscript <- list()
i <- 0
write("lines ", file=logfile, append=TRUE)
repeat{
	i <- i+1
	if (i > length(thenewscript)) break
	write(paste(i,'\r'), file=logfile, append=TRUE)
	if ((i %% 20) == 0)
	{
		write('\n', file=logfile, append=TRUE)
	}
# "\r" suppresses the linefeed
	sumfound <- sum(found)
	theline <-  textBeforeHash(thenewscript[i])
	theShortline <-  textBeforeHash(theline, dropBlanks=TRUE)
	if (any(grep("<-sienaAlgorithmCreate(", theShortline, fixed=TRUE, useBytes=TRUE)))
	{
		found["sienaAlgorithmCreate"] <- found["sienaAlgorithmCreate"]+1
		TheLines <- findTheLines(i, thenewscript)
	# the command sienaAlgorithmCreate is executed
	# to have the algorithm object available for future calls of siena07
		eval(parse(text=TheLines[[1]]))
		if (keepOldAlgorithms)
		{
			thenewerscript <- c(thenewerscript,
			"### The following call of sienaAlgorithmCreate is redundant")
			for (j in (i:TheLines[[2]]))
			{
				thenewerscript <- c(thenewerscript, thenewscript[j])
			}
		}
		i <- TheLines[[2]]
		write(paste('\nline ', i,'; found and executed', "sienaAlgorithmCreate"),
							file=logfile, append=TRUE)
		write("lines ", file=logfile, append=TRUE)
	}
	theTransform <- transformCommand("siena07", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("sienaBayes", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("simulateData", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("varCovar", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("coDyadCovar", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("includeEffects", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("setEffect", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("includeInteraction", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("sienaRI", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("Multipar.RSiena", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline)) 
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("Wald.RSiena", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("testSame.RSiena", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	theTransform <- transformCommand("score.Test", thenewscript, i, logfile)
	i <- theTransform$i
	if (!is.null(theTransform$newline))
	{
		thenewerscript <- c(thenewerscript, theTransform$newline)
	}
	if (sumfound == sum(found))
	{
		thenewerscript <- c(thenewerscript, thenewscript[i])
	}
}


write(c("", "", "PASS 3"), append=TRUE, file=logfile)
# Replace remaining occurrences of function names in furtherNamePairs
# in lines containing no "<-" (presumably after "#")"<-",
# or for "Wald.RSiena", "testSame.RSiena", "score.Test" without a "("", 
# by the new names:

i <- 1
while((i <= length(thenewerscript))){
	thenewerscript[i] <- Reduce22(thenewerscript[i], i )
	i <- i+1
}

thenewerscript <- c("#  Script for RSiena",
		paste("#", "  transforming '", deparse1(substitute(theoldscript)),
					"' to the new names of version 1.6.", sep=""),
		paste("#  log file : ", deparse1(substitute(logfile)), sep=""),
		thenewerscript)

if (sum(found) + sum(furtherfound) > 0)
{
	write(c("", "SUMMARY", ""), append=TRUE, file=logfile)
}

if (sum(found) > 0)
{
write("Total found and processed:", file=logfile, append=TRUE)
	write(paste(names(found)[found>0],  found[found>0],
					sep = ": ", collapse = "\n"), file=logfile, append=TRUE)
}

if (sum(furtherfound) > 0)
{
	write(c("", "Total further found and replaced:"), file=logfile, append=TRUE)
	write(paste(names(furtherfound)[furtherfound>0],
				furtherfound[furtherfound>0],
					sep = ": ", collapse = "\n"), file=logfile, append=TRUE)
}

script <- unlist(thenewerscript)
write(script,  file=fileName)
if (logfile=="")
{
	logfile <- "console"
}
if (fileName=="")
{
	fileName <- "console"
}
cat("New script written to",  fileName, "; log file is",
		logfile, ". \n")
invisible(script)
}
