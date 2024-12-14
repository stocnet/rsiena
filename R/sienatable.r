## /***************************************************************************
##  * SIENA: Simulation Investigation for Empirical Network Analysis
##  *
##  * Web: https://www.stats.ox.ac.uk/~snijders/siena
##  *
##  * File: sienatable.r
##  *
##  * Description: This file contains the code to save a latex or html table of
##  * estimates for a sienaFit or sienaBayesFit object
##  * Written by Charlotte Greenan; modifications by Tom Snijders.
##  *
##  ***************************************************************************/


##@shortBayesResult abbreviated sienaBayesFit results
shortBayesResults <- function(x, nfirst=NULL){
	if (!inherits(x, "sienaBayesFit"))
	{
		stop('x must be a sienaBayesFit object')
	}
	if (is.null(nfirst))
	{
		nfirst <- x$nwarm+1
	}
	df1 <- sienaFitThetaTable(x, fromBayes=TRUE, nfirst=nfirst)[[1]][,
		c("text", "value", "se", "cFrom", "cTo", "postSd", "cSdFrom", "cSdTo" )]
	df1$postSd[is.na(df1$cSdFrom)] <- NA
	df1$postSd <- as.numeric(df1$postSd)
	df1$cSdFrom <- as.numeric(df1$cSdFrom)
	df1$cSdTo <- as.numeric(df1$cSdTo)
	df2 <- as.data.frame(x$requestedEffects[,c("name","shortName", "interaction1", "interaction2",
		"type", "randomEffects", "fix", "parm", "period", "effect1", "effect2", "effect3", "group")])
	df2$period <- as.numeric(df2$period)
	replace1 <- function(x){ifelse(x=="text", "effectName", x)}
	replace2 <- function(x){ifelse(x=="value", "postMeanGlobal", x)}
	replace3 <- function(x){ifelse(x=="se", "postSdGlobal", x)}
	replace4 <- function(x){ifelse(x=="postSd", "postSdBetween", x)}
	dfs <- cbind(df2, df1)
	dfr <- dfs
	names(dfr) <- replace1(replace2(replace3(replace4(names(dfs)))))
	dfr
}

##@siena.table siena07 Saves latex or html table of estimates
## for a sienaFit or sienaBayesFit object
siena.table <- function(x, type='tex',
	file=paste(deparse(substitute(x)),'.',type,sep=""),
	vertLine=TRUE, tstatPrint=FALSE,
	sig=FALSE, d=3, nfirst=NULL)
{
	## fromObjectToText transforms names from the effects object
	## to names readable in LateX or text.

    fromObjectToText <- function(a, type='tex'){
		b <- as.character(a)
		if (type == 'tex')
		{
			b <- gsub('->', '$\\rightarrow$', fixed=TRUE, b)
			b <- gsub('<-', '$\\leftarrow$', fixed=TRUE, b)
			b <- gsub('<>', '$\\leftrightarrow$', fixed=TRUE, b)
			b <- gsub('=>', '$\\Rightarrow$', fixed=TRUE, b)
			b <- gsub('>=', '$\\geq$', fixed=TRUE, b)
			b <- gsub('<=', '$\\leq$', fixed=TRUE, b)
			b <- gsub('sqrt', '$\\sqrt{}$', fixed=TRUE, b)
			b <- gsub('&', '.', fixed=TRUE, b)
		}
# Note: R changes \\ to \ but still displays \\ in printing the string.
		b <- gsub('^(1/1)', '', fixed=TRUE, b)
		b <- gsub('^(1/2)', '(sqrt)', fixed=TRUE, b)
		b <- gsub('^', '', fixed=TRUE, b)
		b <- gsub('_', '-', fixed=TRUE, b)
		b <- gsub('#', '.', fixed=TRUE, b)
		b
	}

	objectName <- deparse(substitute(x))
	fileName <- file
	fromBayes <- FALSE
	xkind.string <- "sienaFit"
	tstat <- tstatPrint
	if (!inherits(x, "sienaFit"))
	{
		if (inherits(x, "sienaBayesFit"))
		{
			fromBayes <- TRUE
			x$gmm <- FALSE # for compatibility
			sig <- FALSE
			tstat <- TRUE # for sienaBayesFit, the place for t-stats is used for between-groups sd.
			if (is.null(nfirst))
			{
				nfirst <- x$nwarm + 1
				message("Note: the print function for sienaBayesFit objects can also use a parameter nfirst,")
				message("    indicating the first run from which convergence is assumed.")
				message("    The default value used now is nfirst =",
					x$nwarm + 1, ".")
			}
			xkind.string <- "sienaBayesFit"
		}
		else
		{
			stop('x must be a sienaFit or sienaBayesFit object')
		}
	}
	if (!(type %in% c('tex','html')))
	{
		stop('type should be either "tex" or "html"')
	}
	if (!gmm(x))
	{
	  effects <- x$requestedEffects
	  p <- x$pp
	}
	else
	{
	  effects <- x$requestedEffects[!x$gmmEffects,]
	  p <- x$pp - sum(x$gmmEffects)
	}
	condrates <- 0
	if (fromBayes)
	{
		nwaves <- attr(x$initialResults$f,'observations')
	}
	else
	{
	nwaves <- dim(x$targets2)[2]
	}
	if (x$cconditional)
	{
		condrates <- length(x$rate)
	}

	if (fromBayes)
	{
		xx <- shortBayesResults(x, nfirst=nfirst)
		theta <- xx$postMeanGlobal
		ses <- xx$postSdGlobal
		sd.between <- xx$postSdBetween
		se.string <- '(psd)'
		sdb.string <- 'psd-between'
	}
	else
	{
	  if (!x$gmm)
	  {
	    theta <- x$theta
#	    theta[diag(x$covtheta) < 0.0 | x$fixed] <- NA
# This is misleading, theta could be nonzero for fixed
	    ses <- sqrt(diag(x$covtheta))
	    ses[x$fixed] <- NA
	    max.t1 <- max(abs(x$tstat[!x$fixed]))
	    dd <- 2
	    max.t <- round(max.t1, digits = dd)
	  }
	  else
	  {
	    fixednogmm <- x$fixed[!x$gmmEffects]
	    theta <- x$theta
#	    theta[diag(x$covtheta) < 0.0 | fixednogmm] <- NA
# This is misleading, theta could be nonzero for fixed
	    theta <- theta[!x$gmmEffects]
	    ses <- sqrt(diag(x$covtheta))
	    ses[fixednogmm] <- NA
	    max.t1 <- max(abs(x$tstat[!fixednogmm]))
	    dd <- 2
	    max.t <- round(max.t1, digits = dd)
	  }
	if (max.t < max.t1)
	{
		max.t <- max.t + 10^{-dd} #needs to be rounded up
	}
	maxlincomb.t1 <- x$tconv.max
	maxlincomb.t <- round(maxlincomb.t1, digits = dd)
	if (maxlincomb.t < maxlincomb.t1)
	{
		maxlincomb.t <- maxlincomb.t + 10^{-dd} #needs to be rounded up
	}
		se.string <- '(s.e.)'
	}
	if (length(x$condvarno) == 0)
	{
		condvarno <- 0
	}
	else
	{
		condvarno <- x$condvarno
	}

	effects$effectName <- fromObjectToText(effects$effectName, type=type)
	extraType <- effects$type
	extraType <- vapply(extraType, function(tp){switch(tp,
			'creation'='creation', 'endow'='maintenance', '')},
			FUN.VALUE='a', USE.NAMES =FALSE)
	effects$effectName <- paste(effects$effectName, extraType)
	max.eff.width <- max(nchar(effects$effectName))
	effects$effectName <- format(effects$effectName,width=max.eff.width)

	max.width <- function(theta)
	{
		max(nchar(as.character(round(abs(theta[!is.na(theta)]))))) +
			2*(min(theta[!is.na(theta)])<0)
	}

	max.ses.width <- max.width(ses)
	max.theta.width <- max.width(theta)
	max.tstat.width <- ifelse(fromBayes, max.width(sd.between), max.width(theta/ses))

	## signif converts t values into daggers and asterisks

	signif <- function(a)
	{
		s <- format("",width=17)

		if (!is.na(a))
		{
			a <- abs(a)
			signif1 <- qnorm(1-0.5*c(0.001,0.01,0.05,0.1))

			if (type=="html")
			{
				signif2 <- c("&#134","*","**","***")
			}
			else
			{
				signif2 <- c(format("$^\\dagger$",width=18),
					format("$^\\ast$",width=18),
					format("$^{\\ast\\ast}$",width=19),
					"$^{\\ast\\ast\\ast}$")
			}
			s2 <- signif2[sum(a>signif1)]

			if (length(signif2[sum(a>signif1)])>0)
			{
				s <- s2
			}
		}
		s
	}

	## mystr rounds a number and splits into its integer and fractional parts

	mystr <- function(theta,int.width=0)
	{
		if (!is.na(theta))
		{
			tc <- as.character(round(theta,d))
			tcsplit <- unlist(strsplit(tc,split=".",fixed=T))

			if (length(tcsplit)==1)
			{
				tcsplit <- c(tcsplit,paste(rep("0",d),sep="",collapse=""))
			}

			if (nchar(tcsplit[2])<d)
			{
				tcsplit[2] <- paste(tcsplit[2],paste(rep("0",d-nchar(tcsplit[2])),
						sep="",collapse=""),sep="")
			}

			if (theta < 0)
			{
				if (tcsplit[1] == "0")
				{
					if (type=='tex')
					{
						tcsplit[1] <- "--0"
					}
					else
					{
						tcsplit[1] <- "&nbsp;&ndash;0"
					}
				}
				else
				{
					if (type=='tex')
					{
						tcsplit[1] <- paste(c("-",tcsplit[1]),sep="",collapse="")
					}
					else
					{
						tcsplit[1] <- sub("-","&nbsp;&ndash;",tcsplit[1])
					}
				}
			}
		}
		else
		{
			if (fromBayes)
			{
				tcsplit <- c("","")
			}
			else
			{
			tcsplit <- c("N","A.")
		}
		}

		if (int.width>0)
		{
			tcsplit[1] <- format(tcsplit[1],width=int.width,justify="right")
		}
		tcsplit
	}

	## mydf creates a data.frame; these will be bound together to form the table

	mydf <- function(pp)
	{
		if (type == "html")
		{
		data.frame(first=rep("", pp),
			effect=rep("", pp),
			amp1 = rep("", pp),
			par1 = rep("", pp),
			amp2 = rep("", pp),
			par2 = rep("", pp),
			amp3s = rep("", pp),
			signif = rep("",pp),
			amp3lpar = rep("", pp),
			se1 = rep("", pp),
			amp4 = rep("", pp),
			se2 = rep("", pp),
			rpar = rep("", pp),
			amp5 = rep("", pp),
			tstat1 = rep("",pp),
			amp6 = rep("",pp),
			tstat2 = rep("",pp),
			ent = rep("",pp),
			stringsAsFactors =FALSE)
		}
		else
		{
		data.frame(first=rep("", pp),
			effect=rep("", pp),
			amp1 = rep("", pp),
			par1 = rep("", pp),
			amp2 = rep("", pp),
			par2 = rep("", pp),
			signif = rep("",pp),
			amp3lpar = rep("", pp),
			se1 = rep("", pp),
			amp4 = rep("", pp),
			se2 = rep("", pp),
			rpar = rep("", pp),
			amp5 = rep("", pp),
			tstat1 = rep("",pp),
			amp6 = rep("",pp),
			tstat2 = rep("",pp),
			ent = rep("",pp),
			stringsAsFactors =FALSE)
	}
	}

	## mydf2 creates a data.frame with latex symbols required
	## on each line of the main body of the table

	mydf2 <- function(pp)
	{
		if (type == "html")
		{
			if (tstat)
			{
				amp5 <- rep("</TD><TD align=\"right\">",pp)
				amp6 <- rep(".",pp)
			}
			else
			{
				amp5 <- rep("", pp)
				amp6 <- rep("", pp)
			}

			data.frame(first = rep("<TR><TD>", pp),
				effect = rep("", pp),
				amp1 = rep("</TD><TD align=\"right\">", pp),
				par1 = rep("", pp),
				amp2 = rep(".", pp),
				par2 = rep("", pp),
				amp3s = rep("</TD><TD align=\"left\">", pp),
				signif = rep("",pp),
				amp3lpar = rep("</TD><TD align=\"right\">&nbsp;&nbsp;&nbsp;(", pp),
				se1 = rep("", pp),
				amp4 = rep(".", pp),
				se2 = rep("", pp),
				rpar = rep(")", pp),
				amp5 = amp5,
				tstat1 = rep("",pp),
				amp6 = amp6,
				tstat2 = rep("",pp),
				ent = rep("</TD></TR>",pp),
				stringsAsFactors =FALSE)
		}
		else
		{
			if (tstat == TRUE)
			{
				amp <- rep(" & ",pp)
			}
			else
			{
				amp <- rep("", pp)
			}

			data.frame(first=rep("", pp),
				effect=rep("", pp),
				amp1 = rep(" & ", pp),
				par1 = rep("", pp),
				amp2 = rep(" & ", pp),
				par2 = rep("", pp),
				signif = rep("",pp),
				amp3lpar = rep(" & (", pp),
				se1 = rep("", pp),
				amp4 = rep(" & ", pp),
				se2 = rep("", pp),
				rpar = rep(")", pp),
				amp5 = amp,
				tstat1 = rep("",pp),
				amp6 = amp,
				tstat2 = rep("",pp),
				ent = rep("\\\\",pp),
				stringsAsFactors =FALSE)
		}
	}

	## tableSection creates lines of latex which can be appended to the table;
	## eg. headings, subtitles

	tableSection <- function(latex)
	{
		r <- mydf(length(latex))
		r$effect <- latex
		r
	}

	## lines of latex to include

	if (type == "html")
	{
		if (tstat)
		{
			if (fromBayes)
			{
				start.tstat <- "<TD>betw. sd.</TD>"
			}
			else
		{
			start.tstat <- "<TD>t stat.</TD>"
		}
		}
		else
		{
			start.tstat <- ""
		}

		startTable <- tableSection(c("<TABLE border=1  rules=none frame = hsides>",
				paste("<TR><TD>Effect</TD><TD>&nbsp;&nbsp;&nbsp;par.</TD><TD></TD>
					<TD>&nbsp;&nbsp;&nbsp;",se.string,"</TD>",start.tstat,"</TR>")))
					midTable <- tableSection(c("",""))
					indentTable <- tableSection("")
		ruleTable <- tableSection("<TR> <TD colspan='9'><HR/></TD> </TR>")

		footnoteStart <- c("</TABLE>","<TABLE border=1  rules=none frame = below>")
		if (fromBayes)
		{
			footnote <- c(paste("<TR> <TD colspan=9 align=left> par=posterior mean;
							psd = posterior standard deviation;
							</TD> </TR> "),
						paste("<TR> <TD colspan=9 align=left>
							betw. sd = posterior between-groups stand. deviation.
							</TD> </TR> <TR> </TR>"), "</TABLE>")
		}
		else
		{
					footnote <- c(paste(" <TR> <TD colspan=9 align=left>
							all convergence t ratios < ",
							max.t,".</TD> </TR> <TR> </TR>",
							" <TR> <TD colspan=9 align=left>
							Overall maximum convergence ratio ",
							maxlincomb.t,".</TD> </TR> <TR> </TR>",
							sep="",collapse=""),"</TABLE>")
		}
		if (sig)
					{
			footnote <- c(footnoteStart, "<TR> <TD colspan=4 align=left> &#134 p < 0.1;
							* p < 0.05; ** p < 0.01; *** p < 0.001;  </TD> </TR> <TR> </TR> " ,footnote)
					}
		else
		{
			footnote <- c(footnoteStart, footnote)
		}

	}
	else
	{
		if (vertLine)
		{
			linesep="|"
		}
		else
		{
			linesep=""
		}
		if (tstat)
		{
			start.tstat <- "r@{.}l"
			if (fromBayes)
			{
				start.tstat2 <- paste("&\\multicolumn{2}{c", linesep, "}{betw. sd}")
		}
		else
		{
				start.tstat2 <- paste("&\\multicolumn{2}{c", linesep, "}{$t$ stat.}")
		}
		}
		else
		{
			start.tstat <- ""
			start.tstat2 <- ""
		}
		if (is.null(x$startingDate))
		{
			startdate <- NULL
		}
		else
		{
			startdate <- paste("%Estimation date",x$startingDate)
		}
		startTable <- tableSection(c(paste("% Table based on", xkind.string, "object",
					objectName, ',', date()),
					startdate,
				paste("\\begin{tabular}{l",
					linesep,
					"r@{.}l r@{.}l",start.tstat,
					linesep,"}"),
				"\\hline",
				"\\rule{0pt}{2ex}\\relax",
				paste("Effect &\\multicolumn{2}{c}{par.}&\\multicolumn{2}{c",
					"}{",se.string,"}",
					start.tstat2,"\\\\[0.5ex]"),
				"\\hline"))
		midTable <- tableSection(c("\\hline",
				"\\rule{0pt}{2ex}\\relax"))
		indentTable <- tableSection("\\rule{0pt}{2ex}\\relax")
		ruleTable <- tableSection("\\hline")
		if (fromBayes)
		{
			footnote <- c(paste("\\multicolumn{7}{l}\n   ",
				"{\\footnotesize{par = posterior mean; psd = posterior standard deviation;}}\\\\\n",
				"\\multicolumn{7}{l}\n   ",
				"{\\footnotesize{betw. sd = posterior between-groups stand. deviation.}}\\\\\n",
				sep="",collapse=""),
			"\\end{tabular}")
		}
		else
		{
		footnote <- c(paste("\\multicolumn{5}{l}\n   ",
				"{\\footnotesize{convergence $t$ ratios all $<$ ", max.t,
				".}}\\\\\n",
				"\\multicolumn{5}{l}",
				"{\\footnotesize{Overall maximum convergence ratio ",
				maxlincomb.t,".}}",
				sep="",collapse=""),
			"\\end{tabular}")

		if (sig == TRUE)
		{
			footnote <- c("\\multicolumn{5}{l}{\\footnotesize{$^\\dagger$ $p$ $<$ 0.1;
				$^\\ast$ $p$ $<$ 0.05; $^{\\ast\\ast}$ $p$ $<$ 0.01;
				$^{\\ast\\ast\\ast}$ $p$ $<$ 0.001;}}\\\\" ,footnote)
		}
	}
	}

	endTable <- tableSection(footnote)

	subtitleLatex <- function(subtitle)
	{
		if (type == "html")
		{
			tableSection(paste("<TR> <TD colspan=4 align=left>",subtitle,
					"  </TD> </TR> <TR> </TR>"))
		}
		else
		{
			tableSection(c(paste("\\multicolumn{",4+2*tstat,"}{l}{\\emph{",subtitle,"}}&\\\\",sep=""),
					"\\hline"))
		}
	}

	## main body of table for each dependent variable

	mainLatex <- function(rows,sections)
	{
		nn <- length(rows)

		if (condvarno == sections)
		{
			nrate <- 0
			m <- 2
		}
		else
		{
			nrate <- sum(effects[rows,]$type == 'rate')
			m <- 0

			if (nn>nrate)
			{
				m <- 2
			}
		}

		mainTable <- rbind(mydf2(nrate),mydf(m),mydf2(nn-nrate))
		mid <- nrate+(1:2)

		if (m == 2)
		{
			mainTable[mid,] <- midTable
		}
		mainTable$effect[-mid] <- effects[rows,]$effectName
		mainTable$par1[-mid] <- sapply(theta[rows],mystr,max.theta.width)[1,]
		mainTable$par2[-mid] <- sapply(theta[rows],mystr)[2,]
		mainTable$se1[-mid] <- sapply(ses[rows],mystr,max.ses.width)[1,]
		mainTable$se2[-mid] <- sapply(ses[rows],mystr)[2,]

		basicRates <- (1:nwaves)
		if (!x$gmm)
		{
		  fixed.2 <- (1:nn)[x$fixed[rows]]
		}
		else
		{
		  fixed.2 <- (1:nn)[fixednogmm[rows]]
		}

		if (condvarno == sections)
		{
			basicRates <- nn+m+1 #so none are removed
		}

		remove <- unique(c(basicRates,fixed.2))

		if (tstat)
		{
			if (fromBayes)
			{
				remove <- is.na(sd.between)
				mainTable$tstat1[-mid][-remove] <-
					sapply(sd.between[rows],mystr,max.tstat.width)[1,][-remove]
				mainTable$tstat2[-mid][-remove] <-
					sapply(sd.between[rows],mystr)[2,][-remove]

				if (any(remove))
				{
					if (type=='tex')
					{
# This had to be dropped, gave a warning; or that is perhaps not serious?
#						mainTable$tstat1[-mid][remove] <- "\\omit"
#						mainTable$tstat2[-mid][remove] <- "-"
					}
				}
			}
			else
		{
			mainTable$tstat1[-mid][-remove] <-
				sapply(theta[rows]/ses[rows],mystr,max.tstat.width)[1,][-remove]
			mainTable$tstat2[-mid][-remove] <-
				sapply(theta[rows]/ses[rows],mystr)[2,][-remove]

			if (min(remove)<nn+m+1)
			{
				if (type=='tex')
				{
					mainTable$tstat1[-mid][remove] <- "\\omit"
					mainTable$tstat2[-mid][remove] <- "-"
				}
			}
		}

		}

		if (sig)
		{
			mainTable$signif[-mid][-remove] <- sapply(theta[rows]/ses[rows],signif)[-remove]
		}

		if (condvarno == sections)
		{
			rateTable <- mydf2(condrates)
			rateTable$effect <- paste('Rate',1:condrates)
			max.rate.width <- max.width(x$rate)
			max.vrate.width <- max.width(x$vrate)
			rateTable$par1 <- sapply(x$rate,mystr,max.rate.width)[1,]
			rateTable$par2 <- sapply(x$rate,mystr)[2,]
			rateTable$se1 <- sapply(x$vrate,mystr,max.vrate.width)[1,]
			rateTable$se2 <- sapply(x$vrate,mystr)[2,]

			if (tstat)
			{
				if (type=='tex')
				{
					rateTable$tstat1 <- "\\omit"
					rateTable$tstat2 <- "-"
				}
			}

			rbind(indentTable,rateTable,mainTable,ruleTable)
		}
		else
		{
			rbind(indentTable,mainTable,ruleTable)
		}
	}

	## Builds the table

	sections <- 0
	table <- startTable

	nBehavs <- sum(x$f$types == "behavior")
	nNetworks <- length(x$f$types) - nBehavs

	if (nBehavs > 0 && nNetworks > 0)
	{
		table <- rbind(table,subtitleLatex("Network Dynamics"))
	}

	if (nNetworks > 0)
	{
		netEffects <- effects[effects$netType != "behavior",]
		netNames <- unique(netEffects$name)
		sections <- 0
		for (i in 1:nNetworks)
		{
			sections <- sections+1
			thisNetTable <- mainLatex((1:p)[effects$name == netNames[i]],sections)
			table <- rbind(table,thisNetTable)
		}
	}

	if (nBehavs > 0 && nNetworks > 0)
	{
		table <- rbind(table,subtitleLatex("Behaviour Dynamics"))
	}

	if (nBehavs > 0)
	{
		behEffects <- effects[effects$netType == 'behavior',]
		behNames <- unique(behEffects$name)

		for (i in 1:nBehavs)
		{
			sections <- sections+1
			thisBehTable <- mainLatex((1:p)[effects$name == behNames[i]],sections)
			table <- rbind(table,thisBehTable)
		}
	}

	table <- rbind(table,endTable)

	##Saves the table to a file
	write.table(table,file=fileName,row.names=F,col.names=F,sep="", quote=FALSE)
	cat('Results for', objectName, 'written to', fileName,'.\n')
}
