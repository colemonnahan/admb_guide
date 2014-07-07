
library(devtools)
library(roxygen2)
library(coda)
library(R2admb)
simple.fit <- read_admb("simple/simple")
document("admbtools")
load_all("admbtools")
dev_help("pairs_admb")                  # how to see doc
## devtools::create("admbtools")

## Demonstrate run.mcmc and pairs_admb. The former runs chains, and the
## later is similar to pairs() but works specifically for ADMB model fits.
simple1 <- run.mcmc("simple", "simple", Nout=1000, mcsave=1, burn.in=1,
                    verbose=TRUE)
pairs_admb(posterior=simple1$mcmc, fits=simple1$fit)
pairs_admb(posterior=simple1$mcmc, fits=simple1$fit, diag="trace")
simple2 <- run.mcmc("simple", "simple", Nout=1000, mcsave=100, burn.in=1, verbose=TRUE)
pairs_admb(posterior=simple2$mcmc, fits=simple2$fit)
pairs_admb(posterior=simple2$mcmc, fits=simple2$fit, diag="trace")


## age.ctl contains inputted values for the phases of the 4 parameters, in
## order of: K, r, S0, Splus. Show progression of adding estimation for
## these
write.table(x=c(1,-1,-1,-1), file="age/age.ctl", row.names=F, col.names=F)
mcmc1 <- run.mcmc("age", "age", Nout=500, mcsave=500, burn.in=10)
write.table(x=c(1,1,-1,-1), file="age/age.ctl", row.names=F, col.names=F)
mcmc2 <- run.mcmc("age", "age", Nout=500, mcsave=500, burn.in=10)
write.table(x=c(1,1,1,-1), file="age/age.ctl", row.names=F, col.names=F)
mcmc3 <- run.mcmc("age", "age", Nout=500, mcsave=500, burn.in=10)
write.table(x=c(1,1,1,1), file="age/age.ctl", row.names=F, col.names=F)
mcmc4 <- run.mcmc("age", "age", Nout=500, mcsave=500, burn.in=10)
pairs_admb(posterior=mcmc1$mcmc, fits=mcmc1$fit)
pairs_admb(posterior=mcmc2$mcmc, fits=mcmc2$fit)
pairs_admb(posterior=mcmc3$mcmc, fits=mcmc3$fit)
pairs_admb(posterior=mcmc4$mcmc, fits=mcmc4$fit)

## Use simple model to explore what happens as the MLE approaches the bound

setwd("simple")
write.table(x=10, file="bounds.txt", row.names=FALSE, col.names=FALSE)



## explore SIR with simple model
system("simple/simple.exe -mcmc 10 -mcsave 1")
system("simple/simple.exe -mceval")
sir <- read.csv("simple/MCMCreport.csv", header=T)
if(file.exists(xx <- "simple/MCMCreport.csv"))
    file.remove(xx)
psv <- file("simple/simple.psv", "rb")
nparams <- readBin(psv, "integer", n=1)
mcmc <- matrix(readBin(psv, "numeric", n=nparams*12),
               ncol=nparams, byrow = TRUE)
close(psv)


psv <- file("age/age.psv", "rb")
nparams <- readBin(psv, "integer", n=1)
mcmc <- matrix(readBin(psv, "numeric", n=nparams*1000),
               ncol=nparams, byrow = TRUE)
close(psv)

pairs_admb(mcmc)


draw.priors <- function(n){cbind(runif(n, -2, max=5), runif(n, 2, max=6))}

set.seed(1)
priors <- t(matrix(t(draw.priors(5)), nrow=1, byrow=F))

psv <- file("simple/simple.psv", "wb")
writeBin(nparams, psv)
writeBin(as.vector(priors), psv)
close(psv)

xx <- c(3.7881 , 0 ,4.05181 , 0.00379198 ,4.31553 , 0 ,4.57924 , 0.00758395 ,4.84296 ,
0.0113759 ,5.10667 , 0.0227519 ,5.37039 , 0.0568796 ,5.6341 , 0.0492957 ,5.89782
, 0.140303 ,6.16153 , 0.0910074 ,6.42525 , 0.0910074 ,6.68896 , 0.216143
,6.95267 , 0.170639 ,7.21639 , 0.348862 ,7.4801 , 0.546045 ,7.74382 , 0.329902
,8.00753 , 0.356446 ,8.27125 , 0.238894 ,8.53496 , 0.318526 ,8.79868 , 0.227519
,9.06239 , 0.257854 ,9.32611 , 0.200975 ,9.58982 , 0.0455037 ,9.85354 ,
0.0417117 ,10.1173 , 0.0151679 ,10.381 , 0.00379198 ,10.6447 , 0)

mat <- matrix(xx, ncol = 2, byrow = T)
plot(mat)
sum(mat[,2])
