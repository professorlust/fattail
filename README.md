# fattail

R code implementation of [Scalable semiparametric inference for the means of heavy-tailed distributions](https://arxiv.org/abs/1602.08066)

All the code you need to apply the bootstrap or posterior inference algorithms is in `R/fatlib.R`.  The other files implement our simulation study.

Below is a demonstration of posterior inference.  See comments for detail.

	 source("R/fatlib.R")
	 
	 ## simulate some heavy tailed data
	 N <- 1000
	 z <- c(rexp(N, 0.1), rexp(N, 0.1)+rgpd(N,scale=10,shape=1/2))
	 
	 ## run the inference algorithm
	 u <- 50 # threshold, you need to choose
	 fit <- meanInference(z, u=u)
	 # (if you get warnings from gpdLV, its because you don't have 
	 # enough data above u for a good Laplace approx; Emu and SDmu still OK)

	 names(fit) # uses variable names from the paper

	 # this should be in the neighborhood of u; if not, choose a higher threshold
	 fit$sigma/fit$xi

	 # posterior mean and SD on the full distribution mean
	 fit$Emu # mean
	 fit$SDmu # SD

	 # Full posterior inference:
	 # * lamSamp is the posterior sample of means for threshold exceedances
	 # * The mean below threshold has approx normal posterior
	 # * probability of surpassing threshold has posterior Beta(n, m)
	 B <- length(fit$lamSamp)
	 below_u_mean <- rnorm(B, mean(z[z<u]), sd(z[z<u])/sqrt(fit$m) )
	 prob_exceed <- rbeta(B, fit$n, fit$m)

	 # combine to get the sample for distribution mean
	 muSamp <- below_u_mean*(1-prob_exceed) + (u+fit$lamSamp)*prob_exceed

	 # plot full sample against a Gaussian with posterior mean and var
	 mm <- seq(min(muSamp),max(muSamp),length=100)
	 hist(muSamp, col=8, freq=FALSE, main="posterior for mu", xlab="mu")
	 lines(mm, dnorm(mm, fit$Emu, fit$SDmu), lwd=2, col="dodgerblue")

This should produce a nice picture of the posterior _for_ the distribution mean.  

![image](https://cloud.githubusercontent.com/assets/4673515/25775366/71c49a48-3271-11e7-90c5-804830f1bbf9.png)

The histogram is a full posterior sample, and the blue line shows a Gaussian parametrized by the posterior moments.
