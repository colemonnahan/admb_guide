
source("C:/Users/Cole/Dropbox/Research/R Functions/ADMB_Covariance_Functions.R")
source("../functions.R")
library(ellipse)
library(coda)

## age structured exmaple
setwd("blue_age")
load(".RData")
make.png=TRUE

## First run one and make sure it converges, for plotting
age.full.list <- run.mcmc("age", se.scale=1, thin.every=1, mcsave=1000,
                          burn.in=1000, mcscale=TRUE, Nmcmc=5000000, cor.max=.5,
                          verbose=T)
age.full <- age.full.list$df
age.full.mcmc <- age.full.list$mcmc[seq(1, nrow(age.full.list$mcmc),
                                              by=1),]
age.mle <- read.admbFit("age")
names(age.full.mcmc) <- age.mle$names[-5]
age.names <- age.mle$names[-5]
acf(age.full.mcmc)
plot(age.full.mcmc[,4], type="l")

if(make.png) png("../Plots/age_pairs.png", 9,6, "in", res=300)
pairs(age.full.mcmc, lower.panel=my.panel.age, diag.panel=panel.hist,
      upper.panel=panel.cor)
if(make.png) dev.off()

## Now run them across a sequence
se.scale.seq <- seq(.005,4, len=30)
Nmcmc <- 200000
thin.every <- 100
burn.in <- 1000
age.list <- list()
for(ss in 1:length(se.scale.seq)){
    age.list[[ss]] <- run.mcmc("age", se.scale=se.scale.seq[ss], thin.every=thin.every,
                               burn.in=burn.in, mcscale=FALSE, Nmcmc=Nmcmc,
                               verbose=F)
    print(ss)
}
age <- do.call(rbind, age.list)
se.best <- 16 #which.max(age$efs4)
age.default.list <-  run.mcmc("age", se.scale=1, thin.every=thin.every,
                            burn.in=burn.in, mcscale=TRUE, Nmcmc=Nmcmc,
                            verbose=T)
age.mcrb.list <- run.mcmc("age", se.scale=1, thin.every=thin.every,
                            burn.in=burn.in, mcscale=TRUE, Nmcmc=Nmcmc,
                            verbose=T, mcrb=5)
age.mcprobe.list <- run.mcmc("age", se.scale=1, thin.every=thin.every,
                            burn.in=burn.in, mcscale=TRUE, Nmcmc=Nmcmc,
                            verbose=T, mcgrope=.4)
age.best.list <- run.mcmc("age", se.scale=se.scale.seq[se.best], thin.every=thin.every,
                              burn.in=burn.in, mcscale=FALSE, Nmcmc=Nmcmc,
                              verbose=T)
age.relaxed.list <- run.mcmc("age", se.scale=se.scale.seq[se.best], thin.every=thin.every,
                              burn.in=burn.in, mcscale=FALSE, Nmcmc=Nmcmc, cor.max=.95,
                              verbose=T)
age.default <- age.default.list$df
age.mcrb <- age.mcrb.list$df
age.mcprobe <- age.mcprobe.list$df
age.best <- age.best.list$df
age.relaxed <- age.relaxed.list$df
age.best/age.default

for(param in 1:4){
    if(make.png)
        png(paste0("../Plots/age_performances_",age.names[param],".png"), 5, 6,
            units="in", res=300)
    par(mfrow=c(3,1), mar=c(3,2,3,.1))
    acf(age.default.list$mcmc[,param], lwd=3, main=NA)
    mtext("Autocorrelation", line=.5)
    xx <- acf(age.best.list$mcmc[,param], plot=F)
    points(x=xx$lag-.5, y=xx$acf, col="red", type="h", lwd=3)
    legend("topright", legend=c("Default", "Ratio-Optimized"), lty=1, lwd=3,
           col=c("black","red"))
    ylim=c(.9, 1.2)*range(age.best.list$mcmc[,param])
    plot(age.default.list$mcmc[,param], type="l", lwd=.5,ylim=ylim, xlab=NULL, ylab=NULL)
    mtext("Trace of Default MCMC", line=.5)
    plot(age.best.list$mcmc[,param], col=rgb(1,0,0,1), lwd=.5, type="l",ylim=ylim,
         xlab="Index", ylab=NULL)
    mtext("Trace of Ratio-Optimized MCMC", line=.5)
    if(make.png) dev.off()
}

if(make.png) png("../Plots/age_comparisons.png", 5,6, "in", res=300)
par(mfrow=c(3,1), mar=c(3,4,1.5,.1), oma=c(2,0,0,0))
ylim=c(0,.2)
with(age,plot(se.scale, ratio, ylim=c(0,1), pch=16, type="b", ann=F))
mtext("Standard Error Multiplier", 1, 2.5)
mtext("Acceptance Ratio", 2,2.5)
abline(h=age.default$ratio)
text(x=2, y=age.default$ratio, labels="Default Ratio", pos=3, cex=1.25)
with(age,plot(se.scale, efs4, type="b", ylim=ylim, ann=F, pch=16))
mtext("Standard Error Multiplier", 1, 2.5)
mtext("% Effective Sample Size", 2,2.5)
abline(h=age.default$efs4)
abline(h=age.mcrb$efs4, lty=2, col=col.mcrb)
abline(h=age.mcprobe$efs4, lty=2, col=col.mcprobe)
legend("topright", legend=c("Default", "mcrb", "mcprobe"), bty="n",
       col=c("black", col.mcrb, col.mcprobe), lty=c(1,2,2))
with(age,plot(ratio, efs4, type="p", xlim=c(0,1), ylim=ylim, ann=F, pch=16))
mtext("Acceptance Ratio", 1,2.5)
mtext("% Effective Sample Size", 2,2.5)
abline(h=age.mcrb$efs4, lty=2, col=col.mcrb)
abline(h=age.mcprobe$efs4, lty=2, col=col.mcprobe)
abline(h=age.default$efs4)
abline(v=age.default$ratio)
if(make.png) dev.off()




## Repeat but changing the correlation
cor.scale.seq <- seq(.1,.99, len=30)
age2.list <- list()
for(ss in 1:length(cor.scale.seq)){
    age2.list[[ss]] <- run.mcmc("age", se.scale=se.scale.seq[se.best], thin.every=thin.every,
                               burn.in=burn.in, mcscale=FALSE, Nmcmc=Nmcmc,
                                cor.max=cor.scale.seq[ss], verbose=F)
    print(ss)
}

age2 <- do.call(rbind, age2.list)
cor.best <- which.max(age2$efs4)
age2.default.list <-  age.default.list
age2.best.list <- run.mcmc("age", se.scale=se.scale.seq[se.best], thin.every=thin.every,
                              burn.in=burn.in, mcscale=FALSE, Nmcmc=Nmcmc,
                              verbose=T,  cor.max=cor.scale.seq[cor.best])
age2.default <- age2.default.list$df
age2.best <- age2.best.list$df
age2.best/age.default


for(param in 1:4){
    if(make.png)
        png(paste0("../Plots/age2_performances_",age.names[param],".png"), 5, 6,
            units="in", res=300)
    par(mfrow=c(3,1), mar=c(2,2,2,.1))
    acf(age2.default.list$mcmc[,param], lwd=3, main=NA)
    mtext("Autocorrelation", line=.5)
    xx <- acf(age2.best.list$mcmc[,param], plot=F)
    points(x=xx$lag-.5, y=xx$acf, col="red", type="h", lwd=3)
    legend("topright", legend=c("Default", "Correlation-Optimized"), lty=1, lwd=3,
           col=c("black","red"))
    ylim=c(.9, 1.2)*range(age2.best.list$mcmc[,param])
    plot(age2.default.list$mcmc[,param], type="l", lwd=.5,ylim=ylim, xlab=NULL, ylab=NULL)
    mtext("Trace of Default MCMC", line=.5)
    plot(age2.best.list$mcmc[,param], col=rgb(1,0,0,1), lwd=.5, type="l",ylim=ylim,
         xlab="Index", ylab=NULL)
    mtext("Trace of Correlation-Optimized MCMC", line=.5)
    if(make.png) dev.off()
}


if(make.png) png("../Plots/age2_comparisons.png", 5,6, "in", res=300)
par(mfrow=c(3,1), mar=c(3,4,1.5,.1), oma=c(2,0,0,0))
ylim=c(0,1)
with(age2,plot(-cor.scale.seq, ratio, ylim=c(0,.2), pch=16, type="b", ann=F))
mtext("Correlation", 1, 2.5)
mtext("Acceptance Ratio", 2,2.5)
abline(h=age2.default$ratio)
abline(v=age.mle$cor[1,4])
text(x=age.mle$cor[1,4], y=.15, labels="Original Correlation", pos=4)
text(x=2, y=age2.default$ratio, labels="Default Ratio", pos=3, cex=1.25)
with(age2,plot(-cor.scale.seq, efs4, type="b", ylim=ylim, ann=F, pch=16))
mtext("Correlation", 1, 2.5)
mtext("% Effective Sample Size", 2,2.5)
abline(v=age.mle$cor[1,4])
text(x=age.mle$cor[1,4], y=.9, labels="Original Correlation", pos=4)
abline(h=age2.default$efs4)
abline(h=age.mcrb$efs4, lty=2, col=col.mcrb)
abline(h=age.mcprobe$efs4, lty=2, col=col.mcprobe)
legend("topright", legend=c("Default", "mcrb", "mcprobe"), bty="n",
       col=c("black", col.mcrb, col.mcprobe), lty=c(1,2,2))
with(age2,plot(ratio, efs4, type="p", xlim=c(0,.1), ylim=ylim, ann=F, pch=16))
mtext("Acceptance Ratio", 1,2.5)
mtext("% Effective Sample Size", 2,2.5)
abline(h=age.mcrb$efs4, lty=2, col=col.mcrb)
abline(h=age.mcprobe$efs4, lty=2, col=col.mcprobe)
abline(h=age.default$efs4)
abline(v=age.default$ratio)
if(make.png) dev.off()

## Run special cases, why are these so bad??
test1 <- run.mcmc("age", se.scale=se.scale.seq[se.best], thin.every=thin.every,
                               burn.in=burn.in, mcscale=FALSE, Nmcmc=Nmcmc,
                                cor.max=cor.scale.seq[15], verbose=T)
plot(test1$mcmc[,4], type="l")
tail(test1$mcmc[,4])
test2 <- run.mcmc("age", se.scale=1, thin.every=thin.every,
                               burn.in=burn.in, mcscale=TRUE, Nmcmc=Nmcmc,
                                cor.max=cor.scale.seq[15], verbose=T)


## test <- run.mcmc("age", se.scale=2, thin.every=1, mcsave=2,
##                           burn.in=1, mcscale=FALSE, Nmcmc=1000,
##                           verbose=T)
## points(test$mcmc[,1])

