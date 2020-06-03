#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: RSienaRDocumentation.r
# *
# * Description: This module contains the code for documenting the
# * RSiena R source.
# *****************************************************************************/
##
##@getInternals Documentation
getInternals <- function()
{
    fnlist <- read.csv("RSienafnlist.csv", as.is=TRUE)
    mylist <- ls(parent.frame())
    ##  print(mylist)
    ## require(codetools)
	## fnlist has a row number at the front: column 3 contains the function names
    mylist <- mylist[mylist %in% fnlist[, 3]]
    mytt <- lapply(mylist, function(x)
       {
           x <- get(x, envir=parent.frame(3))
           if (is.function(x))
           {
               tt <- codetools::findGlobals(x, merge=FALSE)[[1]]
               tt2 <- codetools::findLocals(body(x))
               tt <- c(tt, tt2)
               tt[tt %in% fnlist[, 3]]
           }
           else
           {
                 NULL
           }
       }
           )
    names(mytt) <- mylist
    mytt
}
##@getRSienaDocumentation Documentation
getRSienaRDocumentation <- function(Rdir)
{
	## require(xtable)
	## require(codetools)

    thisdir <- getwd()
    ## temporarily move directory
	on.exit(setwd(thisdir))
    setwd(Rdir)

    ## extract comment lines
    if (.Platform$OS.type == "windows")
    {
        shell('grep "##@" *.r *.R > comments.lis')
    }
    else
    {
        system('grep "##@" *.r *.R > comments.lis')
    }
    ## read them in
    comms <- readLines('comments.lis')
    ## remove the file
    file.remove("comments.lis")
    ## remove the shell line
    comms <- comms[!grepl("comments.lis", comms)]
	## get rid of tabs
	comms <- gsub("\t", "", comms)
    ## split off
    mystr <- paste("##", "@", sep="")
    comms1 <- strsplit(comms, mystr)
    ## join up rest
    ## comms2 <- do.call(rbind, comms1)
    ## turn into dataframe
    comms3 <- sapply(comms1, function(x)
                 {
                     tmp <- strsplit(x[2], " ")[[1]]
                     if (tmp[2] == "internal")
                     {
                         c(x[1], tmp[1], tmp[2], paste('internal to', tmp[3],
                                                       collapse=" "))
                     }
                     else
                     {
                         c(x[1], tmp[1], tmp[2], paste(tmp[-c(1, 2)],
                                                       collapse=" "))
                     }
                 }
                     )
    comms3  <-  t(comms3)

    ## get the calls (global)
    codet <- lapply(comms3[,2], function(x)
                {
                    x <- try(getFromNamespace(x, pkgname), silent=TRUE)
                    if (is.function(x))
                    {
                        tmp1 <- codetools::findGlobals(x, merge=FALSE)[[1]]
                        tmp2 <- codetools::findLocals(body(x))
                        tmp <- c(tmp1, tmp2)
                    }
                    else
                        tmp <- NULL
                    unique(as.vector(tmp[tmp %in% comms3[,2]]))
                }
                    )
    names(codet) <- comms3[, 2]

    ## now the internal ones
    ## find the list of files from comms3
#	browser()
	internals <- comms3[grepl("internal to", comms3[, 4]), ]
    ttmp <- unique(comms3[grepl("internal to", comms3[, 4]), 4])
    ttmp1 <- sub("internal to ", "", ttmp)
	ttmp3 <- lapply(ttmp1, function(x)
				{
					tmp <- x
					while(!is.na(sub <- match(tmp[1], internals[, 2])))
					{
						tmp1 <- internals[sub, 4]
						tmp1 <- sub("internal to ", "", tmp1)
						tmp <- c(tmp1, tmp)
					}
					tmp
				})
    #ttmp2 <- comms3[match(ttmp1, comms3[, 2]), 1]
    #ttmp2 <- sub(":", "", ttmp2)
    ## write out the fnlist in the Rdir
    write.csv(data.frame(comms3), "RSienafnlist.csv")
    ## get the list of internals
	tt <- lapply(ttmp3, function(x)
			 {
				 yy <- getFromNamespace(x[1], pkgname)
                 targs <- formals(yy)
				 targs[1:length(targs)] <- 1
				 if (length(x) > 1)
				 {
					 targs['getDocumentation'] <- x[-1]
				 }
				 else
				 {
					 targs['getDocumentation'] <- TRUE
				 }
                 do.call(yy, targs)
			 })
    #tt <- lapply(1:length(ttmp3), function(x, y)
    #         {
    #             yy <- y[x]
    #             yy <- getFromNamespace(yy, pkgname)
    #             targs <- formals(yy)
    #             n <- length(targs)
    #             myargs <- targs
    #             for (i in 1:n)
#				 {
#                     myargs[[i]] <- 1#
				 #}
                 #myargs['getDocumentation'] <- TRUE
                 #do.call(yy, myargs)
             #}, y=ttmp)
    names(tt) <- ttmp
    ## remove the file
    file.remove("RSienafnlist.csv")
    ## reformat this
    ttt <- lapply(1:length(tt), function(x,y)
              {
                  yy <- y[[x]]
                  n <- length(y[[x]])
                  bb <- names(yy)
                  t1<- lapply(1:n, function(x,  b, a)
                          {
                              y <- a[[x]]
                              bb <- b[[x]]
                              n <- length(y)
                              if ( n > 0)
                                  cbind( rep(bb, n), y)
                              else
                                  c( bb, " ")
                          },  a=yy, b=bb)
                  do.call(rbind,t1)
              }, y=tt
                  )

    tttt <- as.data.frame(do.call(rbind,ttt))
    names(tttt) <- c('Function', 'Calls')

    ## create an object that will tabify to the right output
    tmp2 <- codet

    tmp4 <- lapply(1 : length(tmp2), function(x, y, z, a)
               {
                   n <- length(y[[x]])
                   if (n > 0)
                   {
                       cbind( rep(a[x, 1], n), rep(z[x], n),  y[[x]],
                             rep(a[x, 3], n), rep(a[x, 4], n))
                   }
                   else
                   {
                       cbind(a[x, 1], z[x], " ", a[x, 3], a[x, 4])
                   }

               }, y=tmp2, z=names(tmp2), a=comms3)

    tmp5 <- do.call(rbind, tmp4)
    tmp5 <- as.data.frame(tmp5, stringsAsFactors=FALSE)
    names(tmp5) <- c('Source File', 'Function', 'Calls', 'Type', 'Notes')

    ## now merge in the internals
    tmp5bit <- tmp5[tmp5$Function %in%tttt$Function,]
    tmerge <- merge(tmp5bit, tttt, by="Function")
    tmerge <- tmerge[, -3]
    tmerge <- tmerge[, c(2, 1, 5, 3, 4)]
    names(tmerge)[3] <- "Calls"
    tmp5new <- rbind(tmp5[!tmp5$Function %in% tttt$Function,], tmerge)
    tmp55 <- split(tmp5new, tmp5new$Function)

    ## same for called by
    tmp6 <- lapply(1 : length(tmp2), function(x, y, z)
               {
                   n <- length(y[[x]])
                   if (n > 0)
                   {
                       cbind( rep(z[x], n),  y[[x]] )
                   }
                   else
                   {
                       cbind(z[x], " ")
                   }

               }, y=tmp2, z=names(tmp2))

    tmp7 <- do.call(rbind, tmp6)

    tmp7 <- as.data.frame(tmp7, stringsAsFactors=FALSE)
    names(tmp7) <- c( 'Called from', 'Function')

    tmp7 <- tmp7[order(tmp7[,2],tmp7[,1]), ]

    tttt7 <- tttt
    names(tttt7) <- c("Called from", "Function")
    tttt7 <- tttt7[order(tttt7[,2],tttt7[,1]), ]

   ## tmp7bit <- tmp7[tmp7$Function %in% tttt7$Function, ]

    tmp7new <- merge(tmp7, tttt7, by=c("Function", "Called from"), all=TRUE)

    tmp7new <- tmp7new[order(tmp7new[,1], tmp7new[,2]),]

    tmp77new <- split(tmp7new, tmp7new$Function)

    tmp77new <- tmp77new[-1]

    ## create desired output format
    tmp11 <- lapply(1:length(names(tmp55)), function(x,y,z)
                {
                    thisone <- names(tmp55)[x]
                    yy <- y[[thisone]]
                    zz <- z[[thisone]]
                    d <- max(nrow(yy), nrow(zz))
                    fn <- yy$Function[1]
                    src<- yy$`Source File`[1]
                    type<- yy$Type[1]
                    notes<- yy$Notes[1]
                    if (!is.null(zz))
                    {
                        called <- c(zz[,2], rep(' ', d-nrow(zz)))
                    }
                    else
                    {
                        called <- rep(' ', d)
                    }
                    tmp <- data.frame(src=rep(src,d),
                                      fun=rep(fn, d),
                                      type=rep(type, d),
                                      notes=rep(notes, d),
                                      calls=c(yy[,3], rep(' ', d-nrow(yy))),
                                      called=called, stringsAsFactors=FALSE)
                    tmp
                }, y=tmp55, z=tmp77new)
    ## join into a data frame
    tmp12 <- do.call(rbind, tmp11)
    names(tmp12)[2] <- "Function"

    tmp12 <- tmp12[order(tmp12[, "type"], as.numeric(row.names(tmp12))), ]
    tmp12 <- tmp12[, c(3, 2, 5, 6, 4, 1)]
    ff <- xtable::xtable(tmp12)
    ## go back to start directory
    setwd(thisdir)
    print(ff, tabular.environment="longtable",
          file="RSienaRDocumentation.tex", floating=FALSE)

    write.csv(tmp12, "RSienaRDocumentation.csv")
}
