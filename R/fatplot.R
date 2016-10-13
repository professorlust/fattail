ugrid <- c(0,1,10,20,30,40,50,60,70,80,90,
	100,200,300,400,500,600,700,800,900,1000)
xigrid <- c(.25,.5,.75)

load("fatsim.rda")

err <- sapply(PP, function(m) m[,"err"]^2)
ratio <- sapply(PP, function(m)  m[,"ratio"])
ratio[1,] <- NA

uf <- ugrid
uf[1] <- 1/10 # replace zero for log plotting
cols <- c("blue", "mediumvioletred", "orange")

pdf("fatsim.pdf",width=6,height=6)
par(mai=c(.9,.9,.2,.9))
plot(uf, rep(1,length(uf)), type="n", bty="n", log="xy",
	ylim=range(err), xlab="cutoff", ylab="squared error (solid)")
for(j in 1:3){
	lines(uf, err[,j], col=cols[j], lwd=2)
}
par(new=T)
plot(uf, rep(1,length(uf)), type="n", bty="n", log="xy", axes=FALSE,
	ylim=range(na.omit(ratio)), xlab=NA, ylab=NA)
axis(side=4, at=c(1,10,100,1000))
abline(h=1, col=8)
mtext(side = 4, line = 3, 'ratio (dashed)')
for(j in 1:3){
	jna <- is.na(ratio[,j])
	lines(uf[!jna], ratio[!jna,j], col=cols[j], lwd=2, lty=2)
}
legend("topleft", h=T, title="tail index",  legend=xigrid, col=cols, bty="n", lwd=4)
dev.off()