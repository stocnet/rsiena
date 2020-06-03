maxlikefn<- function(z, x, INIT=FALSE, TERM=FALSE, data, effects=NULL,
	nstart=1000,
	pinsdel=0.6, pperm=0.3, prelins=0.1, multfactor=2.0,
	promul=0.1, promul0=0.5, pdiaginsdel=0.1,
	fromFiniteDiff=FALSE, noSamples=1, sampInterval=50, int=1)
{
	mlInit <- function(z, x, data, effects)
	{
		f <- NULL
		if (!inherits(data, "siena"))
		{
			stop("not valid siena data object")
		}
		if (is.null(effects))
		{
			effects <- getEffects(data)
		}
		if(!inherits(effects, "sienaEffects"))
		{
			stop("effects is not a legitimate Siena effects object")
		}
		effects <- effects[effects$include, ]
		z$theta <- effects$initialValue
		z$fixed <- effects$fix
		z$test <- effects$test
		z$pp <- length(z$test)
		z$posj <- rep(FALSE, z$pp)
		z$targets <- rep(0, z$pp)
		##  effectsNames<- getEffectNames(effects)
		z$posj[grep("basic", effects$effectName)] <- TRUE
		z$posj[grep("constant", effects$effectName)] <- TRUE
		z$BasicRateFunction <- z$posj
		observations <- data$observations
		mats <- vector("list", observations)
		f$mynets <- vector("list", observations)
		types <- sapply(data$depvars, function(x) attr(x, "type"))
		netsubs <- which(types=="oneMode")
		netsub <- min(netsubs) ### only one for now
		actsubs <- which(types=='behavior')
		for (i in 1:observations)
		{
			mats[[i]] <- data$depvars[[netsub]][, , i]
			f$mynets[[i]] <- mats[[i]]
			if (i==1)
			{
				f$mynets[[i]][is.na(mats[[i]])] <- 0
			}
			else ##carry missing forward!
			{
				f$mynets[[i]][is.na(mats[[i]])] <-
					f$mynets[[i - 1]][is.na(mats[[i]])]
			}
			f$mynets[[i]][mats[[i]]==10] <- 0
			f$mynets[[i]][mats[[i]]==11] <- 1
		}
		f$mystructs <- vector('list',observations)
		for (i in 1:observations)
		{
			f$mystructs[[i]] <- mats[[i]]
			f$mystructs[[i]][, ] <- 0
			f$mystructs[[i]][mats[[i]]==11] <- 1
			f$mystructs[[i]][mats[[i]]==10] <- 1
		}
		f$mats <- mats
		for (i in 1:observations)
		{
			f$mats[[i]][mats[[i]]==11] <- 1
			f$mats[[i]][mats[[i]]==10] <- 0
		}
		if (length(actsubs) > 0)
		{
			acts <- matrix(data$depvars[[actsubs[1]]],
				ncol=observations)
			f$acts <- acts
			f$myacts <- acts
			f$myacts[is.na(acts)] <- 0
			f$meanact <- round(mean(acts, na.rm=TRUE))
		}
		f$observations <- observations
		## browser()

		if (any(z$targets!=0))
		{
			Report(c("Targets should be zero for maximum likelihood:',
					'they have been zeroed\n"))
					z$targets <- rep(0, z$pp)
		}
		mat1 <- data$depvars[[netsub]][, , 1]
		mat2 <- data$depvars[[netsub]][, , 2]
		## f$mat1<- mat1
		## f$mat2<- mat2
		startmat <- mat1
		startmat[is.na(startmat)] <- 0
		endmat <- mat2
		endmat[is.na(endmat)] <- startmat[is.na(endmat)]
		diffmat <- startmat != endmat
		if (is.null(x$multfactor))
			f$niter <- multfactor * sum(diffmat)
		else
			f$niter <- x$multfactor * sum(diffmat)
		### create initial chain
		chain <- matrix(0, nrow=sum(diffmat), ncol=4)
		chain[,1] <- row(diffmat)[diffmat]
		chain[,2] <- col(diffmat)[diffmat]
		chain <- chain[sample(1:nrow(chain)),]
		chain[, 4] <- 1:nrow(chain)
		##chain<- chain ##(here you can put a known chain in (eg from
		##delphi!)
		cat(nrow(chain), '\n')
		### initialise
		pinsdel <- pinsdel/(1 - pperm)
		pdiaginsdel <- pdiaginsdel/(1 - pperm)
		iter <- 0
		##burnin
		###construct a max like object to be passed to FRAN
		f$startmat <- startmat
		f$endmat <- endmat
		f$chain <- chain
		f$accepts <-  rep(0,7)
		f$rejects <- rep(0,7)
		f$probs <- c(pinsdel, 0, pdiaginsdel)#
		f$madechain <- FALSE
		f$numm <- 20
		for (i in 1:nstart)
		{
			iter <- iter+1
			##   cat(iter,'\n')
			f <- mhstep(z$theta, f, promul, prelins)
		}
		f$madechain <- TRUE
		pinsdel <- pinsdel * (1-pperm)
		pdiaginsdel <- pdiaginsdel * ( 1-pperm)
		f$probs <- c(pinsdel, pperm, pdiaginsdel)
		f$mats <- f$mystructs <- f$mynets <- NULL
		FRANstore(f)
		z
	}

	##start of function
	scores <- NULL
	dff <- NULL
	theta <- z$theta
	##   f<- z$f
	if (INIT)
	{
		z <- mlInit(z, x, data, effects)
		## f <<-f
		return(z)
	}
	else if (TERM)
	{
		f <-  FRANstore()
		z$f <- f
		return(z)
	}
	else
	{
		f <- FRANstore()
		niter <- f$niter
		##  nactors <- nrow(f$startmat)
		promul <- promul0
		#  int <- x$int
		if (z$Phase==2)
		{
			f$accepts <-  rep(0, 6)
			f$rejects <- rep(0, 6)
			varmat <- FALSE
			## browser()
			if (z$nit == 1)## beginning of a subphase
			{
				i <- 1
				while (i <= 10 * niter)
				{
					f <- mhIntStep(theta, f, promul, prelins, int)
					i <- i + f$n
				}
			}
			Z <- vector("list", noSamples)
			i <- 1
			while (i <= niter)
			{
				f <- mhIntStep(theta, f, promul, prelins, int)
				i <- i + f$n
			}
			Z[[1]] <- f$chain
			if (noSamples > 1)
			{
				for (sample in 2:noSamples)
				{
					i <- 1
					while (i < sampInterval)
					{
						f <- mhIntStep(theta, f, promul, prelins, int)
						i <- i + f$n
					}
					Z[[sample]] <- f$chain
				}
			}
			ans <- calcgrad(theta, Z, f$startmat, varmat)
			# browser()
			f$Z <-  Z
			f$chain <- f$Z[[noSamples]]
		}
		else
		{
			varmat <- TRUE
			i <- 1
			while (i <= niter)
			{
				f <- mhIntStep(theta, f, promul, prelins, int)

				i <- i + f$n
			}
			ans <- calcgrad(theta, f$chain, f$startmat, varmat)
		}

		## browser()
		scores <- ans$sc
		dff <- ans$dff
		##cat(object.size(f),  object.size(f$chain), object.size(f$startmat), '\n')
		FRANstore(f)
		#  cat(scores,'\n')
		##browser()
		list(fra=matrix(scores, nrow=1), ntim0 = NULL, feasible = TRUE, OK = TRUE, dff=dff, accepts=f$accepts, rejects= f$rejects)

	}
}
mhIntStep <- function(theta, f, promul, prelins, int)
{
    ff <- lapply(1:int, function(x, theta, f, promul, prelins)
                 mhstep(theta, f, promul, prelins),
                 theta=theta, f=f, promul=promul,
                 prelins=prelins)
    acts <- sapply(ff, function(x)x$act)
    steptypes <- sapply(ff, function(x)x$steptype)
    if (any(acts))
    {
        firstAct <- min(which(acts))
        n <- firstAct
        f$chain <- ff[[n]]$chain
        f$act <- c(rep(FALSE, n-1), TRUE)
        f$steptype <- steptypes[1:n]
        f$accepts[steptypes[n]] <- f$accepts[steptypes[n]]+1
        f$numm <- ff[[n]]$numm
    }
    else
    {
        n <- int
        f$chain <- ff[[n]]$chain
        f$act <-  rep(FALSE, n)
        f$steptype <- steptypes
        f$numm <- ff[[n]]$numm
    }
    rejects <- f$steptype[!f$act]
    rejNos <- table(rejects)
    rejs <- as.numeric(names(rejNos))
    f$rejects[rejs] <- f$rejects[rejs] + as.vector(rejNos)
    f$n <- n
    f
}
mhstep <- function(theta, f, promul, prelins)
{
    mhtryinsertdiag<- function(i,kpos)
    {
        tmpnet<- startmat
        mychain<- chain[1:(kpos-1),,drop=FALSE]
        ## tmpnet[mychain[,1:2]]<-
        ##     1-tmpnet[mychain[,1:2]]
        for (mysub in 1:nrow(mychain))
        {
            ii <- mychain[mysub,1]
            jj <- mychain[mysub,2]
            tmpnet[ii,jj] <- 1- tmpnet[ii,jj]
        }
        diag(tmpnet)<- 0
        ps <- calcprobs(i,tmpnet,betapar,nactors)
        pr <- 1/sum(ps)
        ndiag <- nrow(chain[chain[, 1] == chain[, 2], , drop=FALSE])
        ##   ans <- kappasigmamu(nactors,nrow(chain),lambda,add1=TRUE)
        ##  mu<- ans$mu + 1/lambda/nactors
        ##  sigma2 <- ans$sigma2 + 1/lambda/lambda/nactors/nactors
        prr<- pr*lambda*nactors/(ndiag+1 )
        ##    prdelphi <- log(pr) - log(nactors) - ans$kappa - 0.5*log(sigma2) -
        ##        (1-mu)^2/2/sigma2
        ##    prr <- exp(prdelphi) * (nrow(chain)+1)/(ndiag+1) * nactors
##cat(prr,'\n')
        if (runif(1)<prr)
            accept<-TRUE
        else
            accept<-FALSE
        accept
    }
    mhtrycanceldiag<- function(i,dpos)
    {
        tmpnet<- startmat
        mychain<- chain[1:(dpos-1),,drop=FALSE]
        ## tmpnet[mychain[,1:2]]<-
        ##    1-tmpnet[mychain[,1:2]]
        for (mysub in 1:nrow(mychain))
        {
            ii <- mychain[mysub,1]
            jj <- mychain[mysub,2]
            tmpnet[ii,jj] <- 1- tmpnet[ii,jj]
        }
        diag(tmpnet)<-0
        ps<-calcprobs(i,tmpnet,betapar,nactors)
        pr<- 1/sum(ps)
        ## browser()
        ndiag<-nrow(chain[chain[,1]==chain[,2],,drop=FALSE])
        prr<- ndiag/pr/lambda/nactors
        ## ans <- kappasigmamu(nactors,nrow(chain+1),lambda,add1=TRUE)
        ## mu<- ans$mu - 1/lambda/nactors
        ## sigma2 <- ans$sigma2 - 1/lambda/lambda/nactors/nactors
        ## prdelphi <- -log(pr) + log(nactors) - ans$kappa - 0.5*log(sigma2) -
        ##     (1-mu)^2/2/sigma2
        ## prr <- exp(prdelphi) / (nrow(chain))*(ndiag) / nactors
##cat(prr,'\n')
        if (runif(1)<prr)
            accept<-TRUE
        else
            accept<-FALSE
        accept
    }
    mhinsdelperm<- function(insdel)
    {
        findCCPs <- function(mat)
        {
            if (nrow(mat) > 0)
            {
                dontuse <- c(diff(mat[, 4]), 0) == 1
                ccps <- mat[!dontuse, , drop=FALSE]
                nc <- nrow(ccps) - 1
            }
            else
            {
                ccps <- NULL
                nc <- 0
            }
            list(ccps=ccps, nc=nc)
        }
        ##fix up numm
        numm <- f$numm
        if (numm > nrow(chain))
        {
            numm <- nrow(chain)
        }
        if (numm > 40)
        {
            numm <- 40
        }
        if (numm < 2)
        {
            numm <- 2
        }
        num <- trunc(numm)
        if (!f$madechain)
        {
            num <- 0
        }
        ##  else
        ##      cat('num=', num, 'numm=', numm, '\n')
        if (insdel)
        {
            mults <- as.matrix(unique(chain[duplicated(chain[, 1:2]) &
                                          chain[, 1]!=chain[, 2], , drop=FALSE]))
            nmul <- nrow(mults)
            if (is.null(nmul))
                nmul <- 0
            if (nmul > 0 && runif(1) < promul)
            {
                ##choose one of unique duplicates
                mypair <- mults[sample(1:nmul, 1), ]
            }
            else
            {
                mypair <- sample(1:nactors, 2)
            }
            from <- mypair[1]
            to <- mypair[2]
            similar <- chain[chain[, 1] == from & chain[, 2]==to, , drop=FALSE]
            nk <- nrow(similar)
            ##nk is number of this connection
            tmp <- findCCPs(similar)
            nc <- tmp$nc
            ccps <- tmp$ccps
            ## nc<- nrow(similar[similar[,3]>0,,drop=FALSE])/2
            if (nc < 1)
            {
                ins <- TRUE
            }
            else
            {
                if (runif(1) < prelins)
                {
                    ins <- TRUE
                }
                else
                {
                    ins <- FALSE
                }
            }
            del <- !ins
            if (ins)
            {
                if (nk==0)
                {
                    kmypair<- sample(1:(nrow(chain)+1),2)
                    numav1<- nrow(chain)+1
                    numav2<- nrow(chain)
                    k1<- min(kmypair)
                    k2<- max(kmypair)
                    ##  newsub<- 1
                }
                else
                {
                    ##   if (nc>0)
                    ##   {
                    ##       ccps <- chain[chain[, 1]==from & chain[, 2]==to &
                    ##                    chain[, 3] > 0, ]
                    ##       ccpsubs <- unique(ccps[, 3])
                    ##       subslist <- 1:(max(ccpsubs) + 1)
                    ##       newsub <- min(subslist[!subslist %in% ccpsubs])
                    ##   }
                    ##   else
                    ##   {
                    ##       newsub <- 1
                    ##   }
                    samplist <- 1:nrow(chain)
                    samplist <- samplist[!(from==chain[, 1] & to==chain[, 2])]
                    samplist <- c(samplist, nrow(chain) + 1)
                    numav1 <- length(samplist)

                    k1 <- sample(samplist, 1)
                    ## find last previous and next of this kind
                    thiskind <- (1:nrow(chain))[from==chain[, 1] & to==chain[, 2]]
                    k <- max(c(0, thiskind[thiskind<k1])) + 1
                    kk <- min(c(thiskind[thiskind>k1], nrow(chain) + 1))
                    if (k>nrow(chain))   ##not possible to proceed
                    {
                        return(list(chain=chain, accept=FALSE))
                    }
                    numav2 <- kk - k

                    k2 <- sample(1:numav2, 1)
                    if (k2==k1)
                    {
                        k2 <- kk
                    }
                    if (k2 < k1)
                    {
                        kk <- k2
                        k2 <- k1
                        k1 <- kk
                    }
                }
            }
            else ##del
            {
                ##  possible<- chain[chain[,1]==from&chain[,2]==to&chain[,3]>0,]
                ##  subs<- unique(possible[,3])
                subs <- 1:nc
                if (length(subs)>1)
                    ##       kmypair<- sample(subs,1)
                    kccp <- sample(subs, 1)
                else
                    kccp <- subs
                ##   kmypair<- subs
                pair <- ccps[kccp : (kccp + 1), 4]
                ##  pair<- chain[chain[,1]==from&chain[,2]==to&
                ##              chain[,3]==kmypair,4]
                k1<- min(pair)
                k2<- max(pair)
                numav1<-nrow(chain)+1-nrow(similar)
                tempchain <- chain[-pair, ]
                thiskind <- (1:nrow(tempchain))[from == tempchain[, 1] &
                                                to == tempchain[, 2]]
                ##  thiskind<- (1:nrow(chain))[from==chain[,1]&to==chain[,2]]
                k<- max(c(0,thiskind[thiskind<k1]))
                kk<- min(c(thiskind[thiskind>k1],nrow(chain)+1))
                if (k==0 &kk>nrow(chain))
                    numav2<- numav1-1
                else
                    numav2 <- kk-k-3
                ##     numav2<-kk-k-2
            }

        } #if insdel
        else
        {
            ins<- FALSE
            del<- FALSE
            k2<- nrow(chain)+1
            k1<- sample(1:nrow(chain),1)
        }
        ##set up network to just before k1
        tmpnet<- startmat
        if (k1>2)
        {
            mychain<- chain[1:(k1-1),,drop=FALSE]
            ##  tmpnet[mychain[,1:2]]<-
            ##      1-tmpnet[mychain[,1:2]]
            for (mysub in 1:nrow(mychain))
            {
                ii <- mychain[mysub,1]
                jj <- mychain[mysub,2]
                tmpnet[ii,jj] <- 1- tmpnet[ii,jj]
            }
            diag(tmpnet)<- 0
        }
        k1mat<- tmpnet
        ##decide on permutation
        ## permute num steps but not including k2, stop before end
        ##if deleting k1, remove that one first
        if (del)
        {
            num<- min(num,k2-k1-1)
            end<- k1+num
            start<- k1+1
        }
        else
        {
            num <- min(num,k2-k1)
            end <- k1+num-1
            start <- k1
        }
        truncated <- num!=trunc(numm)
        ##calculate probs of existing chain from k1 to end or k2
        if (ins)
            thisend<- k2-1
        else if (del)
            thisend<- k2
        else
            thisend<- end
        prbefore<- rep(0,thisend-k1+1)
        for (link in k1:thisend)
        {
            i<- chain[link,1]
            j<- chain[link,2]
            ps<- calcprobs(i,tmpnet,betapar,nactors)
            prbefore[link-k1+1]<- ps[j]/sum(ps)
            if (i!=j)
                tmpnet[i,j]<- 1-tmpnet[i,j]
        }
        ##reset up network to just before k1
        tmpnet<- k1mat
        if (num > 0)
        {
            permchain<- chain[start:end,,drop=FALSE]

            myorder<- sample(1:nrow(permchain))
            permchain<- permchain[myorder,]
        }
        else
            permchain <- NULL
        ##reconstruct new chain
        if (k1>1)
            tempchain<- chain[1:(k1-1),]
        else
            tempchain<- NULL
        if (ins)
            tempchain<- rbind(tempchain,c(from,to,0,0))
        ## tempchain<- rbind(tempchain,c(from,to,newsub,0))
        tempchain<- rbind(tempchain,permchain)
        if (end<(k2-1))
            tempchain<- rbind(tempchain,chain[(end+1):(k2-1),]) ##check this
        if (ins)
            tempchain<- rbind(tempchain,c(from,to,0,0))
        ## tempchain<- rbind(tempchain,c(from,to,newsub,0))
        if (!del & k2<=nrow(chain))
            tempchain<- rbind(tempchain, chain[k2:nrow(chain),])
        if (del &k2 < nrow(chain))
            tempchain<- rbind(tempchain, chain[(k2+1):nrow(chain),])
        ##calc all probs between k1 and k2 for new chain
        if (ins)
        {
            start<- k1
            end<- k2+1
        }
        else
            if (del)
            {
                start<- k1
                end<- k2-2
            }
            else
            {
                start<- k1
                end<- k1+num-1
            }
        if (end-start+1>0)
        {
            prafter<- rep(0,end-start+1)
            for (link in start:end)
            {
                i<- tempchain[link,1]
                j<- tempchain[link,2]
                ps<- calcprobs(i,tmpnet,betapar,nactors)
                prafter[link-start+1]<- ps[j]/sum(ps)
                if (i!=j)
                    tmpnet[i,j]<- 1-tmpnet[i,j]
            }
        }
        else
            prafter<- 1
        ##now try to get the ratios
        ##  if(del) browser()
        ##  ans<- kappasigmamu(nactors,nrow(chain),lambda,add1=TRUE)
        ##  if (ins)
        ##  {
        ##      mu<- ans$mu+2/nactors/lambda
        ##      sigma2 <- ans$sigma2+2/nactors/nactors/lambda/lambda
        ##      prdelphi <- sum(log(prafter)) - sum(log(prbefore)) -
        ##          2 * log(nactors) -ans$kappa - 0.5 * log(sigma2) -
        ##              (1-mu)^2/2/sigma2
        ##      prdelphi <- exp(prdelphi)
        ##  }
        ##  if (del)
        ##  {
        ##      mu<- ans$mu - 2/nactors/lambda
        ##      sigma2 <- ans$sigma2 - 2/nactors/nactors/lambda/lambda
        ##      prdelphi <- sum(log(prafter))-sum(log(prbefore)) +
        ##          2 * log(nactors) - ans$kappa - 0.5 * log(sigma2) -
        ##              (1-mu)^2/2/sigma2
        ##      prdelphi <- exp(prdelphi)
        ##  }
        ##  if (!insdel)
        ##  {
        ##      prdelphi<- sum(log(prafter))-sum(log(prbefore))
        ##  }
        logprobrat<- 0
        if (ins)
            logprobrat<- 2*log(lambda)-log(nrow(chain)+2)-
                log(nrow(chain)+1)
        if (del)
            logprobrat<- -2*log(lambda)+log(nrow(chain))+
                log(nrow(chain)-1)
        logprobrat<- logprobrat - sum(log(prbefore))
        if (length(prafter)>0)
            logprobrat<- logprobrat +
                sum(log(prafter))
        probrat<- exp(logprobrat)
        pa<- (1-promul)/nactors/(nactors-1)
        ## tempchain[,4] <- 1:nrow(tempchain)
        if (insdel)
        {
            newsimilar <- tempchain[tempchain[, 1] == from &
                                    tempchain[, 2] == to, ,
                                    drop=FALSE]
            newnc <- findCCPs(newsimilar)$nc
        }
        if (ins)
        {
            if (nk<2)
                pp<- (promul/(nmul+1) +pa)/pa*(1-prelins)*
                    ##      numav1*numav2/2/(nc+1)
                    numav1*numav2/2/(newnc)
            else
                pp<- (1-prelins)/numav1*numav2/2/(newnc)
            ##     pp<- (1-prelins)/numav1*numav2/2/(nc+1)
            if (nc>=1) pp <- pp/prelins
        }
        else if (del)
        {
            if (nk<4)
                pp<- pa*2*nc/(promul/nmul +pa)/
                    (1-prelins)/numav1/numav2
            else
                pp<- 1/(1-prelins)/numav1/numav2*2*nc
            ## if (nc>=2)
            if (newnc >= 1)
                pp <- pp*prelins
        }
        else
        {
            pp<- 1
        }
        probrat <- probrat*pp
        ##cat('probrat ',probrat1,' ')
        ##   probrat <- prdelphi * pp
        ## cat(probrat,' ',probrat-probrat1,'\n')
        ## if (insdel & ins)
        ## {
        ##    nn <- 3
        ##  }
        ##else if (insdel & del)
        ##{
        ##    nn <- 4
        ##}
        ##else
        ##{
        ##    nn <- 5
        ##}
        ##  cat(probrat, '\n')
        accept <- FALSE
        if (probrat > 1)
        {
            accept <- TRUE
        }
        else
        {
            if (runif(1) < probrat)
            {
                accept <- TRUE
            }
        }
        if (sum(tempchain[, 3]==1) %% 2!=0)
        {
            browser()
        }
        if (accept)
        {
            chain <- tempchain
        }
        ## cat(num,f$numm,accept , insdel,truncated,accept,'\n')
        ##   if (!insdel)
        ##      browser()
        if (accept && !insdel && !truncated)
        {
            numm <- numm + 0.5
        }
        if (!accept && !insdel && !truncated)
        {
            numm <- numm - 0.5
        }
        #browser()
        list(chain=chain, accept=accept, numm=numm)
    }##end of procedure
    #########################################################################
    ##start of mhstep
    #########################################################################
    #cat('start', f$numm,'\n')
   ## print(table(f$chain))
    startmat<- f$startmat
    nactors<- nrow(startmat)
    chain<- f$chain
    lambda<- theta[1]
    betapar<- theta[-1]
    steptype<- sample(1:3,1,prob=f$probs)
   # cat(steptype,'step type\n')
    if (steptype==1)
    {
        ans<- mhinsdelperm(TRUE)
        act<- ans$accept
        chain<- ans$chain
    }
    else  if (steptype==2)
    {
        if (f$madechain)
        {
            ans<-mhinsdelperm(FALSE)
            act<- ans$accept
            chain<- ans$chain
            f$numm <- ans$numm
        }
        else
        {
            act<- NA
        }
    }
    else ##if (steptype==3)
    {
        actor <- sample(1:nactors,1)
        kpos<- sample(1:(nrow(chain)+1),1)
        act<- mhtryinsertdiag(actor,kpos)

        if (act)
        {
            if (kpos>1)
                tmp<- chain[1:(kpos-1),]
            else
                tmp<- NULL
            tmp<- rbind(tmp,c(actor,actor,0,0))
            if (kpos<=nrow(chain))
                tmp<- rbind(tmp, chain[kpos:nrow(chain),])
            ##if (sum(chain[,3]==1)%%2!=0)
            ##      browser()
            chain <- tmp
            chain[, 4] <- 1:nrow(chain)
        }
        ## }
        ## else
        ## {
        diags<-chain[chain[,1]==chain[,2],,drop=FALSE]
        ndiag<- nrow(diags)
        act<- FALSE
        if (ndiag>0)
        {
            ## browser()
            kk<- sample(1:ndiag,1)
            actor<- diags[kk,1]
            dpos<- diags[kk,4]
            ## if(chain[dpos,1]!=actor) browser()
            act <- mhtrycanceldiag(chain[dpos,1],dpos)
            if (act)
            {
                if (dpos>1)
                    tmp<- chain[1:(dpos-1),]
                else
                    tmp<- NULL
                if (dpos<nrow(chain))
                    tmp<- rbind(tmp, chain[(dpos+1):nrow(chain),])
                if (sum(chain[,3]==1)%%2!=0)
                    browser()
                chain <- tmp
            }
        }
    }
    ##    cat('act',act,'steptype',steptype,' ')
    chain[, 4] <- 1:nrow(chain)
    f$chain <- chain
    f$act <- act
    f$steptype <- steptype
    #cat(f$numm,'\n')
    f
}

calcprobs <- function(i, tmpnet, betapar, nactors)
{
    ss <- matrix(0, nrow=nactors, ncol=2)
    myout <- tmpnet[i, ]
    myin <- tmpnet[, i]
    ss[, 1] <- 1 - 2 * myout
    ss[, 2] <- (1 - 2 * myout) * myin
    ss[i, ] <- 0
    fs <- drop(ss %*% betapar)
    exp(fs)
}

doMoreUpdates <- function(z, x, n)
{
    f <- FRANstore()
    for (i in 1:n)
    {
        ans <- calcgrad(z$theta, f$Z, f$startmat, varmat=FALSE)

        fchange <- as.vector(z$gain * ans$sc %*% z$dinv)
        fchange <- ifelse(z$posj & (fchange > z$theta), z$theta * 0.5, fchange)

        z$theta <- z$theta - fchange
        z$thav <- z$thav + z$theta
        z$thavn <- z$thavn + 1
    }
    z
}
