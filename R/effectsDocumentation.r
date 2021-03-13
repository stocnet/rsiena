#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: effectsDocumentation.r
# *
# * Description: This module contains a function for documenting the shortNames
# * and other fields of an effects object.
# *****************************************************************************/

##@effectsDocumentation Documentation
effectsDocumentation <- function(effects= NULL, type="html",
	display=(type=="html"),
	filename=ifelse(is.null(effects), "effects", deparse(substitute(effects))))
{
	## require(xtable)
	if (is.null(effects)) {
		x <- RSiena::allEffects[, c("effectGroup", "effectName", "shortName",
			"endowment", "interaction1", "interaction2", "parm", "interactionType")]
	} else {
		if (!inherits(effects, "sienaEffects")) {
			stop(paste(deparse(substitute(effects)), "is not a sienaEffects object."))
		}
		x <- as.data.frame(effects[, c("name", "effectName", "shortName", "type",
			"interaction1", "interaction2", "parm", "interactionType")])
	}
#  if (any(x$type=="gmm"))
#  {
#    x <- x[-which(x$type=="gmm"),]
#  }
	storage.mode(x$parm) <- "integer"
	names(x)[4] <- ifelse(is.null(effects), "endow?", "type")
	names(x)[5] <- "inter1"
	names(x)[6] <- "inter2"
	names(x)[8] <- "interactionType"
	if (is.null(effects))
	{
		x$row <- as.integer(row.names(x))
	}
	else
	{
		x$row <- 1:dim(x)[1]
	}
	x <- x[, c(9, 1:8)]

	myorder <- c("nonSymmetricRate",
				"covarNonSymmetricRate",

				"symmetricRate",
				"covarSymmetricRate",

				"bipartiteRate",
				"covarBipartiteRate",

				"behaviorRate",
				"behaviorOneModeRate",
				"behaviorSymmetricRate",
				"covarBehaviorOneModeRate",
				"covarBehaviorSymmetricRate",
				"behaviorBipartiteRate",
				"covarBehaviorRate",
				"continuousRate",

				"nonSymmetricObjective",
				"dyadObjective",
				"covarNonSymmetricObjective",
				"unspecifiedNetInteraction",
				"nonSymmetricNonSymmetricObjective",
				"nonSymmetricSymmetricObjective",
				"nonSymmetricSymmetricSObjective",
				"nonSymmetricBipartiteObjective",
				"covarNetNetObjective",
				"tripleNetworkObjective",
				"dyadANetNetObjective",
				"settingsObjective",

				"symmetricObjective",
				"covarSymmetricObjective",

				"bipartiteObjective",
				"dyadBipartiteObjective",
				"dyadSecondBipartiteObjective",
				"covarBipartiteObjective",
				"unspecifiedNetInteraction",
				"bipartiteNonSymmetricObjective",
				"bipartiteSymmetricObjective",
				"bipartiteBipartiteObjective",

				"unspecifiedNetInteraction",

				"behaviorObjective",
				"behaviorOneModeObjective",
				"behaviorSymmetricObjective",
				"behaviorBipartiteObjective",
				"behaviorOneOneModeObjective",
				"behaviorSymSymObjective",
				"behaviorOneModeSymObjective",
				"behaviorBipBipObjective",
				"covarBehaviorObjective",
				"covarBehaviorNetObjective",
				"dyadBehaviorNetObjective",
				"covarABehaviorBipartiteObjective",
				"covarBBehaviorBipartiteObjective",
				"unspecifiedBehaviorInteraction",
				"continuousFeedback",
				"continuousWiener",
				"continuousIntercept",
				"continuousOneModeObjective",
				"unspecifiedContinuousInteraction")

	mytab <- table(RSiena::allEffects[,1])

	addtorowPos <- cumsum(c(0, mytab[myorder]))[1:length(myorder)]
	addtorowText <- names(mytab[myorder])
	x[is.na(x)] <- "FALSE" ## endow? field

	if (type == "latex") {
		addtorowText <- paste(" \\hline \\multicolumn{4}{l}{", addtorowText,
			"}\\\\ \\hline")
		addtorowText[1] <- paste(" \\endhead ", addtorowText[1], collapse = "")
	} else {
		x[x == ""] <- "<br>"
		addtorowText <- paste(" <TR> <TD colspan=\"8\" >", addtorowText,
			"</TD> </TR>")
	}
	add.to.row <- NULL
	x$interactionType[x$type=='gmm'] <- ''
	if (is.null(effects))
	{
		add.to.row$pos <- lapply(addtorowPos, function(x)x)
		add.to.row$command <- as.vector(sapply(addtorowText, function(x)x))
		order2 <- match(myorder, x[, 2])
		order3 <- as.vector(mytab[myorder])
		order4 <- unlist(apply(cbind(order2, order3), 1,
				function(x) x[1]:(max(x[1] + x[2] - 1, 1))))
		y <- x[order4, -2]
	}
	else
	{
		y <- rbind(x[x$type!='gmm',], x[x$type=='gmm',])
	}
	row.names(y) <- 1:nrow(y)

	if (type =="latex")
	{
		filename2 <- paste(filename,".tex", sep="", collapse="")
		includefile <- paste(filename,".include.tex", sep="", collapse="")
		includepart <- paste(filename,".include", sep="", collapse="")
		print(xtable::xtable(y), add.to.row=add.to.row, file=includefile,
			tabular.environment="longtable", hline.after=c(-1),
			floating=FALSE, include.rownames=FALSE)

		cat(file=filename2, "\\documentclass[12pt,a4paper]{article}\n",
			"\\usepackage[pdftex,dvipsnames]{color}\n",
			"\\usepackage[landscape]{geometry}\n",
			"\\usepackage{longtable}\n",
			"\\pagestyle{empty}\n",
			"\\textheight=7.25in\n",
			"\\topmargin=-1in\n",
			"\\evensidemargin=-0in\n",
			"\\oddsidemargin=-0in\n",
			"\\marginparwidth=-0in\n",
			"\\begin{document}\n",
			"\\include{", includepart, "}\n",
			"\\end{document}\n", sep=""
			)
	}
	else
	{
		filename <- paste(filename,".html", sep="", collapse="")
		print(xtable::xtable(y), add.to.row=add.to.row,
			file=filename,
			type="html", hline.after=c(-1),
			sanitize.text.function=function(x){x},
			floating=FALSE, include.rownames=FALSE,
			html.table.attributes="border=1 cellpadding=3 cellspacing=0")
		if (display) {
			browseURL(paste("file://", getwd(), "/", filename, collapse = "",
				sep = ""))
		}
	}
}
