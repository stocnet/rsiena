calclk<- function(pars,chain,mat1,mat2,startmat)
{
    lambda<- exp(pars[1])
    beta <- pars[2:3]
    nc<- nrow(chain)
    ##   cat(lambda,beta,nc,'\n')
    tmpnet<- startmat
    nactors<- nrow(startmat)
    pvec<- rep(0,nc)
    for (step in 1:nc)
    {
        i<- chain[step,1]
        j<- chain[step,2]
        ss<- matrix(0,nrow=nactors,ncol=2)
        if (!is.na(mat1[i,j])&!is.na(mat2[i,j]))
        {
            for (ii in 1:nactors)
            {
                tmp<- tmpnet[i,]
                if (ii!=i)
                    tmp[ii]<- 1-tmp[ii]
                ss[ii,1]<- sum(tmp)
                ss[ii,2]<- sum(tmp&tmp==tmpnet[,i])
            }
            fs<-drop(ss%*%beta)
            fs<- fs-fs[i]
            ps<- exp(fs)
            ps[is.na(mat1[i,])|is.na(mat2[i,])]<- 0
            ps<- ps/sum(ps)
            ##  ps<- ps*nactors/sum(ps==0) ##
            ##  cat(sum(ps),'\n')
            ##  ps[is.na(mat1[i,])|is.na(mat2[i,])]<- 1/nactors##
            pvec[step]<- ps[j]
        }
        else
            pvec[step]<-0# 1/nactors
        if (i!=j)
            tmpnet[i,j]<- 1-tmpnet[i,j]
    }
    ##   llik<-nactors*lambda - nc*log(lambda) +log(factorial(nc))-
    ##       sum(log(pvec))
    ##  cat(llik,kappa(nactors,nc,lambda)-sum(log(pvec)),'\n')
    llik <- kappa1(nactors,nc,lambda) - sum(log(pvec))
    llik
}
calcfn<- function(pars,Z,startmat,mat1,mat2)
{
    llik<- sapply(Z,function(x,pars,mat1,mat2,startmat)
                  calclk(pars,x$chain,mat1,mat2,startmat),pars=pars,
                  startmat=startmat,mat1=mat1,mat2=mat2)
    ## cat (pars,llik,'\n')
    ##    browser()
    mean(llik)
}
calcgrad<- function(pars,Z,startmat,varmat=FALSE)
{
    derivs<- function(pars,chain, varmat)
    {
        lambda<- pars[1]
        beta <- pars[2:3]
        nc<- nrow(chain)
        nactors<- nrow(startmat)
      ##  cat(lambda,beta,nc,'\n')
        tmpnet<- startmat
        sc<- rep(0,2)
        st<- rep(0,2)
        ss<- matrix(0,nrow=nactors,ncol=2)
        dff<- matrix(0,nrow=3,ncol=3)
                                        ##  browser()
        for (k in 1:nc)
        {
            i<- chain[k,1]
            j<- chain[k,2]
         ##   if (!is.na(mat1[i,j])&!is.na(mat2[i,j]))
            {
                myout<- tmpnet[i,]
                myin<- tmpnet[,i]
                ss[,1]<- 1-2*myout
                ss[,2]<- (1-2*myout)*myin
                ss[i,]<- 0
               fs<-drop(ss%*%beta)
                ps<- exp(fs)
            ##    ps[is.na(mat1[i,])|is.na(mat2[i,])]<- 0
                ps<- ps/sum(ps)
               ## ss<- sweep(ss,2,ss[i,])
               ## browser()
                wmean<- apply(ss*ps,2,sum)
                if (varmat)
                {
                    ssp<- ss*ps
                    w2mean<-sapply(1:nactors,function(x)outer(ssp[x,],ss[x,]))
                    dim(w2mean)<- c(2,2,nactors)
                    w2mean<- apply(w2mean,c(1,2),sum)
                    dff[-1,-1] <- dff[-1,-1] + outer(wmean,wmean) - w2mean
                }
                sc<- sc + ss[j,]-wmean
                st<- st + ss[j,]
            }
            if (i!=j)
                tmpnet[i,j]<- 1-tmpnet[i,j]
            ##  browser()
        }
        ## cdscore<- c(-lambda,sc)
        ## dsum[i,,]<- outer(obs[i,],cdscore)
        ## browser()
          dlambda<- nactors-nc/lambda
      ##  dlambda<- nactors*lambda - nc
        dff[1,1] <- -nc/lambda/lambda
     #   cat(dlambda,-sc,'\n')
        list(sc= c(dlambda,-sc), dff= -dff)
    }
    if (is.list(Z))
    {
       # browser()
        ##   pars[1]<- log(pars[1])
        grad <- lapply(Z, function(x, pars)
                       derivs(pars, x, varmat=FALSE), pars=pars)
        gradsc <- rowMeans(sapply(grad, function(x)x$sc))
        grad <- list(sc=gradsc, dff=NULL)
       # gradsc <- sapply(grad, function(x)x$sc)
       # gradmn <- rowMeans(gradsc)
       # ##   ##    browser()
       # ##   ##  cat('gradient',gradmn,'\n')
       # gradmn # - gradmn/length(Z)
    }
    else
        grad <- derivs(pars, Z, varmat)
    grad
}
kappa1 <- function(n, nc, lambda)
{
    nc <- nc + 1
    mu <- nc/ n / lambda
    sigma2 <- nc/ n / n / lambda / lambda
    logp <- -(1 - mu)^2 / 2 / sigma2 - log(sigma2) / 2 - log(2 * pi) / 2 -
        log(n * lambda)
    ##   pp <- dpois(nc, n * lambda, log=TRUE)
    ##   cat(logp, pp, '\n')
    -logp + nc * log(n)
}
kappasigmamu<- function(n, nc, lambda, add1=FALSE)
{
    if (add1)
    {
        nc <- nc + 1
    }
    mu <- nc / n / lambda
    sigma2 <- nc / n / n / lambda / lambda
    kappa1 <- -(1 - mu)^2 / 2 / sigma2 - log(sigma2) / 2
    list(kappa=kappa1, mu=mu, sigma2=sigma2)
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

