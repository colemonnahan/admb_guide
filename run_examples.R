## Make sure these are updated
library(devtools)
library(roxygen2)
library(coda)
library(matrixcalc)
library(R2admb)

## document("admbtools")
load_all("admbtools")
## dev_help("pairs_admb")                  # how to see doc
## devtools::create("admbtools")

## Demonstrate run_admb_mcmc and pairs_admb. The former runs chains, and
## the later is similar to pairs() but works specifically for ADMB model
## fits.
simple1 <- run_admb_mcmc("simple", "simple", Nout=1000, mcsave=1,
                         burn.in=1, verbose=TRUE)
pairs_admb(admb_mcmc=simple1)
pairs_admb(admb_mcmc=simple1,  diag="trace")
simple2 <- run_admb_mcmc("simple", "simple", Nout=1000, mcsave=100,
                         burn.in=1, verbose=TRUE)
pairs_admb(admb_mcmc=simple2)
pairs_admb(admb_mcmc=simple2,  diag="trace")
## Run one with user supplied covariance
cov(simple2$mcmc)
cov.user <- matrix(c(.05, .2, .2, .9), nrow=2)
simple3 <- run_admb_mcmc("simple", "simple", Nout=1000, mcsave=100,
                         burn.in=1, cov.user=cov.user, verbose=TRUE)
pairs_admb(admb_mcmc=simple3)
## Explore hybrid option
simple.hy1 <- run_admb_mcmc("simple", "simple", Nout=100, mcsave=1,
                         burn.in=1, verbose=TRUE, hybrid=TRUE, hynstep=50,
                            hyeps=.1)
pairs_admb(admb_mcmc=simple.hy1,  diag="trace")


## Run more of the examples. This finance one seems to have covariance
## estimation issues
setwd('examples')
write.table(x=c(1,1,1,1), file='finance/phases.dat', row.names=FALSE,
            col.names=FALSE)
finance1 <- run_admb_mcmc('finance', 'finance', Nout=1000, mcsave=10,
                          burn.in=5)
pairs_admb(finance1)
write.table(x=c(1,1,1,-1), file='finance/phases.dat', row.names=FALSE,
            col.names=FALSE)
finance2 <- run_admb_mcmc('finance', 'finance', Nout=1000, mcsave=10,
                          burn.in=5)
pairs_admb(finance2)
cov.user <- cov(finance2$mcmc)
finance3 <- run_admb_mcmc('finance', 'finance', Nout=1000, mcsave=10,
                          burn.in=5, cov.user=cov.user)
pairs_admb(finance3)

finance4 <- run_admb_mcmc('finance', 'finance', Nout=1000, mcsave=1,
                          burn.in=5, cov.user=cov.user, hybrid=TRUE,
                          hynstep=20, hyeps=.1)
pairs_admb(finance4)


chem-eng <- run_admb_mcmc('chem-eng', 'chem-eng', Nout=1000, mcsave=1000, burn.in=5)
pairs_admb(chem-eng)
setwd('..')


## ------------------------------------------------------------
## Explore what happens when the bound approaches the MLE
setwd("examples")
Nout <- 1000
mcsave <-  100
posterior.list <- fit.list <- list()
bhat <- 4.0782
bound.seq <- seq(bhat*.99, bhat*1.01, len=50)
for(i in 1:length(bound.seq)){
    write.table(x=bound.seq[i], file="simple/bounds.txt", row.names=FALSE,
                col.names=FALSE)
    temp <- run_admb_mcmc(model.path="simple", model.name="simple", Nout=Nout,
                     mcsave=mcsave, burn.in=1, verbose=F,
                     init.pin=c(0,0), mcseed=i)
    fit.list[[i]] <- temp$mle
    posterior.list[[i]] <- temp
}
cors <- unlist(lapply(fit.list, function(x) x$cor[1,2]))
std <-  do.call(rbind, lapply(fit.list, function(x) x$std[1:2]))
est <-  do.call(rbind, lapply(fit.list, function(x) x$est[1:2]))
## Explore what happens when the bound approaches the MLE but use posfun to
## put a soft limit on the upper bound -- does it improve it??
posterior.posfun.list <- fit.posfun.list <- list()
for(i in 1:length(bound.seq)){
    write.table(x=bound.seq[i], file="simple_posfun/bounds.txt", row.names=FALSE,
                col.names=FALSE)
    temp <- run_admb_mcmc(model.path="simple_posfun",
                     model.name="simple_posfun",  Nout=Nout,
                     mcsave=mcsave, burn.in=1, verbose=F,
                     init.pin=c(0,0), mcseed=i)
    fit.posfun.list[[i]] <- temp$mle
    posterior.posfun.list[[i]] <- temp$mcmc
}
cors.posfun <- unlist(lapply(fit.posfun.list, function(x) x$cor[1,2]))
std.posfun <-  do.call(rbind, lapply(fit.posfun.list, function(x) x$std[1:2]))
est.posfun <-  do.call(rbind, lapply(fit.posfun.list, function(x) x$est[1:2]))


par(mfrow=c(2,3))
plot(bound.seq, cors, type='b', ylim=range(c(cors, cors.posfun)))
lines(bound.seq, cors.posfun, type='b', col=2, pch=16)
plot(bound.seq, est[,1], type='b', ylim=range(c(est[,1], est.posfun[,1])))
lines(bound.seq, est.posfun[,1], type='b', col=2, pch=16)
abline(v=bhat)
plot(bound.seq, est[,2], type='b', ylim=range(c(est[,2], est.posfun[,2])))
lines(bound.seq, est.posfun[,2], type='b', col=2, pch=16)
abline(v=bhat, h=bhat)
plot(bound.seq, std[,1], type='b')
lines(bound.seq, std.posfun[,1], type='b', col=2, pch=16)
plot(bound.seq, std[,2], type='b')
lines(bound.seq, std.posfun[,2], type='b', col=2, pch=16)

for(i in 1:length(bound.seq))
    pairs_admb(posterior.list[[i]], diag="acf",
               limits=list(c(1,3), c(-1,4.1)))
for(i in 1:length(bound.seq))
    admb.pairs(posterior.posfun.list[[i]], diag="acf", fits=fit.posfun.list[[i]],
               limits=list(c(1,3), c(-1,4.1)))


## ------------------------------------------------------------
## OLD CODE -- from earlier developmental versions
## ## age.ctl contains inputted values for the phases of the 4 parameters, in
## ## order of: K, r, S0, Splus. Show progression of adding estimation for
## ## these
## write.table(x=c(1,-1,-1,-1), file="age/age.ctl", row.names=F, col.names=F)
## mcmc1 <- run_admb_mcmc("age", "age", Nout=500, mcsave=500, burn.in=10)
## write.table(x=c(1,1,-1,-1), file="age/age.ctl", row.names=F, col.names=F)
## mcmc2 <- run_admb_mcmc("age", "age", Nout=500, mcsave=500, burn.in=10)
## write.table(x=c(1,1,1,-1), file="age/age.ctl", row.names=F, col.names=F)
## mcmc3 <- run_admb_mcmc("age", "age", Nout=500, mcsave=500, burn.in=10)
## write.table(x=c(1,1,1,1), file="age/age.ctl", row.names=F, col.names=F)
## mcmc4 <- run_admb_mcmc("age", "age", Nout=500, mcsave=500, burn.in=10)
## pairs_admb(admb_mcmc=mcmc1$mcmc, fits=mcmc1$fit)
## pairs_admb(admb_mcmc=mcmc2$mcmc, fits=mcmc2$fit)
## pairs_admb(admb_mcmc=mcmc3$mcmc, fits=mcmc3$fit)
## pairs_admb(admb_mcmc=mcmc4$mcmc, fits=mcmc4$fit)


## ## ------------------------------------------------------------
## ## explore SIR with simple model
## system("simple/simple.exe -mcmc 10 -mcsave 1")
## system("simple/simple.exe -mceval")
## sir <- read.csv("simple/MCMCreport.csv", header=T)
## if(file.exists(xx <- "simple/MCMCreport.csv"))
##     file.remove(xx)
## psv <- file("simple/simple.psv", "rb")
## nparams <- readBin(psv, "integer", n=1)
## mcmc <- matrix(readBin(psv, "numeric", n=nparams*12),
##                ncol=nparams, byrow = TRUE)
## close(psv)
## psv <- file("age/age.psv", "rb")
## nparams <- readBin(psv, "integer", n=1)
## mcmc <- matrix(readBin(psv, "numeric", n=nparams*1000),
##                ncol=nparams, byrow = TRUE)
## close(psv)
## pairs_admb(mcmc)
## draw.priors <- function(n){cbind(runif(n, -2, max=5), runif(n, 2, max=6))}
## set.seed(1)
## priors <- t(matrix(t(draw.priors(5)), nrow=1, byrow=F))
## psv <- file("simple/simple.psv", "wb")
## writeBin(nparams, psv)
## writeBin(as.vector(priors), psv)
## close(psv)
## xx <- c(3.7881 , 0 ,4.05181 , 0.00379198 ,4.31553 , 0 ,4.57924 , 0.00758395 ,4.84296 ,
## 0.0113759 ,5.10667 , 0.0227519 ,5.37039 , 0.0568796 ,5.6341 , 0.0492957 ,5.89782
## , 0.140303 ,6.16153 , 0.0910074 ,6.42525 , 0.0910074 ,6.68896 , 0.216143
## ,6.95267 , 0.170639 ,7.21639 , 0.348862 ,7.4801 , 0.546045 ,7.74382 , 0.329902
## ,8.00753 , 0.356446 ,8.27125 , 0.238894 ,8.53496 , 0.318526 ,8.79868 , 0.227519
## ,9.06239 , 0.257854 ,9.32611 , 0.200975 ,9.58982 , 0.0455037 ,9.85354 ,
## 0.0417117 ,10.1173 , 0.0151679 ,10.381 , 0.00379198 ,10.6447 , 0)
## mat <- matrix(xx, ncol = 2, byrow = T)
## plot(mat)
## sum(mat[,2])
## ## ------------------------------------------------------------
