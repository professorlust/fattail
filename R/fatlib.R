## utilities for inference about the GPD posterior and the semiparametric joint model

require(evd)
require(KernSmooth)

## functions for GPD inference

# neg log posterior
nljoint <- function(v, xx, p=NULL){
    if(is.null(p)) p <- c(10,10,0,0)
    if(any(v <= 0)) return(1e20)
    if(v[1]>=1) return(1e20)
    ll <- -(length(xx)+1-p[3])*log(v[2]) - (1/v[1]+1)*sum(log(1+v[1]*xx/v[2])) +
        (p[1]-1)*log(v[1]) + (p[2]-1)*log(1-v[1]) - p[4]*v[2] 
    return(-ll)
}

# profile alpha map and resulting neg log post
alpha <- function(tg, xx) sapply(tg, 
	function(t) (length(xx)+1)/sum(log(1+t*xx)))
nlprof <- function(t, xx){
    a <- alpha(t, xx)
    -log(t*a) + 1/a
}

# gradient of the profile neg log post
gp <- function(t, xx) 
	(length(xx)+1)/t - (alpha(t,xx)+1)*sum(xx/(1+t*xx))

# wrap it all together
gpdMAP <- function(xx, p=NULL){	
    if(is.null(p)){
        tauhat <- try(uniroot(
                function(t) gp(t, xx=xx),
                 c(1e-10,.1), tol=1e-8)$root,
                silent=TRUE)
        if(class(tauhat)!="try-error"){
            ahat <- alpha(tauhat, xx=xx)
            if(ahat>=1) return(list(
                tau=tauhat, alpha=ahat, 
                xi=1/ahat, sigma=1/(tauhat*ahat),
                lambda=1/(tauhat*(ahat-1))))
        }
    }
    if(is.null(p)) p <- c(10,10,0,0)
    theta <- optim(c(.5,1000), 
            function(v) nljoint(v, xx=xx, p=p), control=list(reltol=1e-12))$par
    xihat <- theta[1]
    sighat <- theta[2]
    ahat <- 1/xihat
    lamhat <- sighat/(1-xihat)
	return( 
		list(tau=theta[1]/theta[2],alpha=1/theta[1],
            xi=theta[1],sigma=theta[2], 
			lambda=theta[2]/(1-theta[1])) )
}

# laplace variance approximation
gpdLV <- function(xx, mapfit, p=NULL){
    if(is.null(p)) p <- c(10,10,0,0)
    q = xx*mapfit$xi/(mapfit$lam*(1-mapfit$xi) + xx*mapfit$xi)
    v = -mapfit$lam^2/(length(xx) - p[3] + 1 + (1/mapfit$xi + 1)*sum(q^2-2*q))
    return(v)
}

## Bootstrap (for the MAP)
gpdBoot <- function(xx, B=1000, p=c(10,10,0,0)){
	lamboot <- vector(length=B)
	tauboot <- vector(length=B)
	xiboot <- vector(length=B)
	sigboot <- vector(length=B)
	predboot <- vector(length=B)
	mapfit <- gpdMAP(xx=xx, p=p)
	for(b in 1:B){  
        # cat(b, "")  	
		# generate the data
    	xb <- rgpd(length(xx), scale=mapfit$sigma, shape=mapfit$xi) 
    	bootmap <- gpdMAP(xb, p=p)
        tauboot[b] <- bootmap$tau
    	xiboot[b] <- bootmap$xi
    	lamboot[b] <- bootmap$lambda
    	sigboot[b] <- bootmap$sigma
    	predboot[b] <- rgpd(1, scale=sigboot[b], shape=xiboot[b])
 	}
 	return(data.frame(
    	lambda=lamboot,
    	xi=xiboot,
    	sigma=sigboot,
    	tau=tauboot,
    	pred=predboot
    	))
}

## Independence MH MCMC (essentially, re-weight the bootstrap)
gpdIMH <- function(xx, boots=NULL, B=1000, p=c(10,10,0,0)){
	if(is.null(boots)) 
        boots <- gpdBoot(xx, B=B, p=p)
    # the edge solutions have zero prior weight, 
    # but cause issues numerically
    edges <- which(boots$xi>0.999)
    ne <- length(edges)
    if(ne>0)
        boots[edges,] <- boots[-edges,][sample.int(B-ne,ne,replace=TRUE),]

	# fit KDE to the bootstrap density (it sometimes warns...)
	g <- suppressWarnings(
            bkde2D(cbind(boots$xi,boots$sigma),
				c(dpik(boots$xi),dpik(boots$sigma))))

	# wrapper to obtain bootstrap density value
	# v is c(xi, sigma)
	lgprop <- function(v){
    	i <- which.min(abs(v[1]-g$x1))
    	j <- which.min(abs(v[2]-g$x2))
    	return(log(g$fhat[i,j]))
	}
	# posterior evaluation
	lpost <- function(v) -nljoint(v,xx=xx, p=p)
	# run I-MH through the bootstrap sample
	pmove <- rep(1,B)

	xisamp <- boots$xi
	sigsamp <- boots$sigma
	predsamp <- boots$pred
	for(b in 2:B){
	    lgnew <- lgprop(c(xisamp[b],sigsamp[b]))
    	lgold <- lgprop(c(xisamp[b-1],sigsamp[b-1]))
    	lpnew <- lpost(c(xisamp[b],sigsamp[b]))
    	lpold <- lpost(c(xisamp[b-1],sigsamp[b-1]))
    	pmove[b] <- exp(lpnew + lgold - lpold - lgnew)
    	if(runif(1) > pmove[b]){
        	xisamp[b] <- xisamp[b-1]
        	sigsamp[b] <- sigsamp[b-1]
        	predsamp[b] <- predsamp[b-1]
    	} 
    }
    pmove[pmove>1] <- 1
    return(data.frame(
    	lam=sigsamp/(1-xisamp),
    	xi=xisamp, sigma=sigsamp,
    	pred=predsamp, pmove=pmove))
}

## simple mean and var
muAndVar <- function(v){
    n = v[1]
    y = v[2]
    y2 = v[3]
    list(n=n, mu=y/n, var=y2/n - (y/n)^2)
}

## put things together
tailInference <-  function(v, B=1000, p=c(10,10,0,0), plotit=FALSE){

    if(length(v)<2)
        return(list( 
            lamMAP=0, 
            lamMean=0, 
            lamSD=0,
            lamLaplaceSD=0,
            pmove=0,
            xi=0.5, sigma=NA))
 

    map <- gpdMAP(v, p=p)
   
    lamMAP <- map$lam
    lamLaplaceSD <- sqrt(gpdLV(v,map, p=p))

    if(B==0)
        return(list( 
            lamMAP=lamMAP, 
            lamMean=lamMAP, 
            lamSD=lamLaplaceSD,
            lamLaplaceSD=lamLaplaceSD,
            pmove=0,
            xi=map$xi, sigma=map$sigma))

    boot <- gpdBoot(v, B=B, p=p)
    post <- try( gpdIMH(v, boot, B=B, p=p) )
    if(class(post)=="try-error") post <- boot
    xi <- map$xi

    lamMean <- mean(post$lam)
    lamSD <- sqrt(var(post$lam))

    if(plotit)try(
    {
        lg <- seq(-4,4,length=100)*lamSD + lamMean
        par(mfrow=c(1,2),mai=c(.6,.6,.3,.2))
        # approx and full posterior
        plot(lg, dnorm(lg, lamMAP, lamLaplaceSD), 
            bty="n", type="l", col=2, xlab="lambda",ylab="density")
        lines(bkde(post$lam),col=4)
        lines(bkde(boot$lam),col=4,lty=2)
        legend("topleft", col=c(2,4,4),lty=c(1,1,2), lwd=2, 
            legend=c("Laplace","indep-MH","bootstrap"), bty="n")
        # posterior predictive
        qqplot(log(v), log(post$pred),ylab="log(posterior predicted v)",bty="n",
            pch=21,bg=8,cex=.75)
        abline(a=0,b=1)
    })

    return(list( 
        lamMAP=lamMAP, 
        lamMean=lamMean, 
        lamSD=lamSD,
        lamLaplaceSD=lamLaplaceSD,
        pmove=mean(post$pmove),
        xi=map$xi, sigma=map$sigma))
}


meanInference <- function(z, u, 
        zstats=NULL, tail=NULL, 
        B=1000, p=c(10,10,0,0), plotit=FALSE)
{
    v <- z[z>=u]-u
    n <- length(v)

    if(is.infinite(u)){
        return(list(
            n=0, m = length(z), xi=NA, sigma=NA,
            Emu=mean(z), SDmu=sd(z)/sqrt(length(z)), 
            EmuLaplace=NA, SDmuLaplace=NA))
    }

    if(is.null(zstats)){
        m <- length(z)-n
        zsum <- sum(z[z<u])
        z2sum <- sum(z[z<u]^2)
    }
    else{
        m <- zstats["nguid"] - n
        zsum <- zstats["gmb"] - sum(v+u)
        z2sum <- zstats["gmb2"] - sum( (v+u)^2 )
    }
    N <- m+n

    if(is.null(tail)) 
        tail <- tailInference(v, B, p, plotit)

    Emu <- zsum/(m+n) + n*(u + tail$lamMean)/(m+n)
    Vmu <- (z2sum - 2*Emu*zsum + Emu^2 
        + n*(u + tail$lamMean - Emu)^2
        + 2*n^2*(m+n-0.5)*tail$lamSD^2/(m+n))/
            ((m+n)*(m+n+1))  
    SDmu <- sqrt(Vmu)

    EmuLaplace <- zsum/(m+n) + n*(u + tail$lamMAP)/(m+n)
    VmuLaplace <- (z2sum - 2*EmuLaplace*zsum + EmuLaplace^2 
        + n*(u + tail$lamMAP - EmuLaplace)^2
        + 2*n^2*(m+n-0.5)*tail$lamLaplaceSD^2/(m+n))/
            ((m+n)*(m+n+1))  
    SDmuLaplace <- sqrt(VmuLaplace)

    return(list(
        n=n, m=m, 
        xi=tail$xi, sigma=tail$sigma,
        Emu=Emu, SDmu=SDmu, 
        EmuLaplace=EmuLaplace, SDmuLaplace=SDmuLaplace))
}


capmean <- function(z, u, zstats){
    vi <- which(z>u)
    n <- length(vi)
    zstats <- unlist(zstats)
    zsum <- zstats["gmb"] - sum(z[vi]) + n*u
    z2sum <- zstats["gmb2"] - sum(z[vi]^2) + n*u^2
    return( muAndVar( c(zstats["nguid"], zsum, z2sum) ) )
}