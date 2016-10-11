source("R/fatlib.R")
library(parallel)
cl <- makeCluster(4,outfile="")
ugrid <- c(1,10,20,30,40,50,60,70,80,90,
	100,200,300,400,500,600,700,800,900,
	1000,2000,5000,10000,50000,100000,Inf)
xigrid <- c(.2,.3,.4,.5,.6,.7,.8,.9)
clusterExport(cl, "ugrid")

getR <- function(xi, n=100000){
	source("R/fatlib.R")
	print(xi)
	z <- c(rexp(n, 0.1),rexp(n, 0.1)+rgpd(n,scale=10,shape=xi))
	R <- c()
	for(u in ugrid){
		print(u)
		R <- rbind(R, meanInference(z,u, p=c(4,4,0,0)))
	}
	return(as.data.frame(R))
}

perf <- function(j){
	xi <- xigrid[j]
	mu <- 10 + 5/(1-xi)
	err <- unlist(RR[[j]][,"Emu"])-mu
	ratio <- unlist(RR[[j]][,"sigma"])/unlist(RR[[j]][,"xi"])/ugrid
	return(data.frame(err=err, ratio=ratio))
}

RR <- parLapply(cl, xigrid, getR)

PP <- lapply(1:length(xigrid), perf)

save(RR,PP,file="fatsim.rda")

stopCluster(cl)