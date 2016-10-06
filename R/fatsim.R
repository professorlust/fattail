source("R/fatlib.R")

n <- 100000
xi <- 0.5
z <- c(rexp(n, 1/rexp(n)),rgpd(n,scale=1*rexp(n),shape=xi))
print(mu <- 0.5*(1+1/(1-xi)))
print(mean(z))
hist(z)

ugrid <- c(1,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,5000,10000,50000,Inf)
xigrid <- c(.25,.5,.75,.95)
library(parallel)
cl <- makeCluster(4)
clusterExport(cl, "ugrid")

getR <- function(xi, n=100000){
	source("R/fatlib.R")
	z <- c(rexp(n, 1/rexp(n)),rgpd(n,scale=1,shape=xi))
	R <- c()
	for(u in ugrid){
		R <- rbind(R, meanInference(z,u, p=c(5,5,0,0)))
	}
	return(as.data.frame(R))
}

RR <- parLapply(cl, xigrid, getR)

perf <- function(j){
	xi <- xigrid[j]
	mu <- 0.5*(1+1/(1-xi))
	err <- unlist(RR[[j]][,"Emu"])-mu
	ratio <- unlist(RR[[j]][,"sigma"])/unlist(RR[[j]][,"xi"])/ugrid
	return(data.frame(err=err, ratio=ratio))
}

PP <- lapply(1:length(xigrid), perf)

save(RR,PP,file="fatsim.rda")
