## Make sure these are updated
library(devtools)
library(roxygen2)
library(coda)
library(matrixcalc)
library(R2admb)
## document("admbtools")
load_all("admbtools")
dev_help("pairs_admb")                  # how to see doc
## global plotting settings
width <- 7; height <- 5

## Demonstrate run_admb_mcmc and pairs_admb. The former runs chains, and
## the later is similar to pairs() but works specifically for ADMB model
## fits.
simple1 <- run_admb_mcmc("simple", "simple", Nout=1000, mcsave=1,
                         burn.in=1, verbose=TRUE)
pairs_admb(admb_mcmc=simple1)
dev.copy2pdf(width=width, height=height,file="Plots/simple1.pdf")
pairs_admb(admb_mcmc=simple1,  diag="trace")
dev.copy2pdf(width=width, height=height,file="Plots/simple1_trace.pdf")

simple2 <- run_admb_mcmc("simple", "simple", Nout=1000, mcsave=100,
                         burn.in=1, verbose=TRUE)
pairs_admb(admb_mcmc=simple2)
dev.copy2pdf(width=width, height=height,file="Plots/simple2.pdf")
pairs_admb(admb_mcmc=simple2,  diag="trace")
dev.copy2pdf(width=width, height=height,file="Plots/simple2_trace.pdf")

### ------------------------------------------------------------
### The hybrid (Hamiltonian) algorithm

## Demonstrate the hybrid option with simple chain
simple.hy1 <-
    run_admb_mcmc("simple", "simple", Nout=1000, mcsave=1,
                  burn.in=1, hybrid=TRUE, verb=FALSE,
                  hynstep=20, hyeps=.1)
pairs_admb(simple.hy1)
dev.copy2pdf(width=width, height=height,file="Plots/simple_hy1.pdf")

## Make a grid of different parameters of the hybrid to show the impact on
## performance
hynstep.seq <- c(10, 100)
hyeps.seq <- c(.05, .5)
hy.grid <- expand.grid(hynstep.seq, hyeps.seq)
labs <- paste0("hyeps=", hy.grid[,2], "; hynstep=", hy.grid[,1])
simple.hy <- simple.hy2 <- list()
a.seq <- seq(1.5,2.5, len=50)
b.seq <- seq(2,6, len=50)
param.grid <- expand.grid(a=a.seq, b=b.seq)
a.mle <- 1.90909098475
b.mle <- 4.07817738582
yobs <- c(1.4, 4.7, 5.1, 8.3, 9.0, 14.5, 14.0, 13.4, 19.2, 18)
xobs <- c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8)
U <- function(q) log(sum((q[1]*xobs+q[2]-yobs)^2))*length(yobs)/2
NLL <- sapply(1:nrow(param.grid), function(x)
                         U(c(param.grid$a[x], param.grid$b[x])))
NLL <- matrix(NLL, nrow=50)


## Make leapfrog algorithm examples
setwd("examples/simple")
leapfrog.list <- list()
write.table(x=c(2,3), file="init.pin", row.names=F, col.names=F)
for(k in 1:nrow(hy.grid)){
    system(paste("simple -mcmc 1 -hybrid -noest -nohess -hyeps",
                 hy.grid[k,2], "-hynstep", hy.grid[k,1], "-ainp init.pin"))
    leapfrog.list[[k]] <- read.table("leapfrog.csv", sep=",")
}
setwd("../..")
par(mfrow=c(2,2), mar= .5*c(1,1,4,1), oma=c(2.5,2.5,0,0), cex.axis=.8)
for(k in 1:4){
    contour(a.seq, b.seq, NLL, axes=F, col=gray(.5))
    points(a.mle, b.mle, col='red', pch=16)
    leapfrog <- leapfrog.list[[k]][-(1:2),1:2]
    n <- nrow(leapfrog)
    arrows(leapfrog[-n,1], leapfrog[-n,2], leapfrog[-1,1],
           leapfrog[-1,2], length=.05)
    points(leapfrog[c(1,n),], pch=c(16, 1), cex=1.5)
    mtext(labs[k], line=.5)
    if(k %in% c(1,3)) {axis(2, mgp=c(1,.5,0), tck=-.02)
                       mtext("a", side=2, line=1.5, cex=1.1)}
    if(k %in% c(3,4)) {axis(1, mgp=c(1,.5,0), tck=-.02)
                       mtext("b", side=1, line=1.5, cex=1.1)}
    box(col=gray(.5))
}
dev.copy2pdf(width=width, height=height,file="Plots/hybrid_grid_trace.pdf")

## Run some cases with same arguments, but different seeds to show the
## impact of momentum variables.
setwd("examples/simple")
leapfrog.list <- list()
write.table(x=c(2,3), file="init.pin", row.names=F, col.names=F)
set.seed(3)
seeds <- sample(1:1000, 4)
for(k in 1:4){
    system(paste("simple -mcmc 1 -hybrid -noest -nohess -hyeps",
                 .25, "-hynstep", 10, "-ainp init.pin -mcseed", seeds[k]))
    leapfrog.list[[k]] <- read.table("leapfrog.csv", sep=",")
}
setwd("../..")
dev.off()
par(mfrow=c(1,1), mar=c(3,3, .5, .5), mgp=c(1.5, .25, 0), cex.axis=.8 ,tck=-.01)
contour(a.seq, b.seq, NLL, axes=TRUE, col=gray(.5), xlab="a", ylab="b")
points(a.mle, b.mle, col='red', pch=16)
for(k in 1:4){
    leapfrog <- leapfrog.list[[k]][-(1:2),1:2]
    n <- nrow(leapfrog)
    arrows(leapfrog[-n,1], leapfrog[-n,2], leapfrog[-1,1],
           leapfrog[-1,2], length=.05)
    points(leapfrog[c(1,n),], pch=c(16, 1), cex=1.5)
}
box(col=gray(.5))
dev.copy2pdf(width=width, height=height,file="Plots/hybrid_seeds.pdf")

## Run one longer to get good ACF values
for(k in 1:nrow(hy.grid)){
    simple.hy[[k]] <-
        run_admb_mcmc("examples/simple", "simple", Nout=1000, mcsave=1,
                      burn.in=1, hybrid=TRUE, verb=FALSE,
                      hynstep=hy.grid[k,1], hyeps=hy.grid[k,2])
}
## Run one shorter to show the step sizes
for(k in 1:nrow(hy.grid)){
    simple.hy2[[k]] <-
        run_admb_mcmc("examples/simple", "simple", Nout=10, mcsave=1,
                      burn.in=1, hybrid=TRUE, verb=FALSE,
                      hynstep=hy.grid[k,1], hyeps=hy.grid[k,2])
}
par(mfrow=c(2,2), mar= .5*c(1,1,1,1), oma=c(2.5,2.5,0,0), cex.axis=.8)
for(k in 1:4){
    with(simple.hy[[k]]$mcmc, acf(a, ylim=c(-1,1), main=NA,
                                  axes=FALSE, ann=FALSE))
    mtext(labs[k], line=-1.5)
    if(k %in% c(1,3)) {axis(2, mgp=c(1,.5,0), tck=-.02)
                       mtext("acf", side=2, line=1.5, cex=1.1)}
    if(k %in% c(3,4)) {axis(1, mgp=c(1,.5,0), tck=-.02)
                       mtext("lag", side=1, line=1.5, cex=1.1)}
    box(col=gray(.5))
}
dev.copy2pdf(width=width, height=height,file="Plots/hybrid_grid_acf.pdf")
par(mfrow=c(2,2), mar= .5*c(1,1,1,1), oma=c(2.5,2.5,0,0), cex.axis=.8)
for(k in 1:4){
    n <- nrow(simple.hy2[[1]]$mcmc)
    plot(0, type='n', xlim=c(1.5,2.5), ylim=c(2.5, 5.5),
         main=NA, axes=FALSE, ann=FALSE)
    with(simple.hy2[[k]]$mcmc, arrows(a[-n], b[-n], a[-1], b[-1], len=.05))
    mtext(labs[k], line=-1.5)
    if(k %in% c(1,3)) {axis(2, mgp=c(1,.5,0), tck=-.02)
                       mtext("b", side=2, line=1.5, cex=1.1)}
    if(k %in% c(3,4)) {axis(1, mgp=c(1,.5,0), tck=-.02)
                       mtext("a", side=1, line=1.5, cex=1.1)}
    box(col=gray(.5))
}
dev.copy2pdf(width=width, height=height,file="Plots/hybrid_grid_trace2.pdf")
## End of hybrid
### ------------------------------------------------------------


## A tougher posterior is the logistic, since it is non-multivariate normal
## as expected by the MH algorithm.
logistic.mh <- run_admb_mcmc('examples/logistic', 'logistic', Nout=1000, mcsave=100,
                          burn.in=5, mcscale=TRUE)
pairs_admb(admb_mcmc=logistic.mh,  diag="acf")
dev.copy2pdf(width=width, height=height,file="Plots/logistic_mh.pdf")
logistic.mh2 <- run_admb_mcmc('examples/logistic', 'logistic', Nout=1000, mcsave=100,
                          burn.in=5, mcscale=TRUE, cov.user=cov(logistic.mh$mcmc))
pairs_admb(admb_mcmc=logistic.mh2,  diag="acf")
dev.copy2pdf(width=width, height=height,file="Plots/logistic_mh2.pdf")
## Try the hybrid algorithm
logistic.hy <- run_admb_mcmc('examples/logistic', 'logistic', Nout=1000, mcsave=1,
                          burn.in=1, hybrid=TRUE, hyeps=.05, hynstep=100)
pairs_admb(admb_mcmc=logistic.hy,  diag="acf")
dev.copy2pdf(width=width, height=height,file="Plots/logistic_hy.pdf")
logistic.hy2 <- run_admb_mcmc('examples/logistic', 'logistic', Nout=1000, mcsave=1,
                          burn.in=1, hybrid=TRUE, hyeps=.05, hynstep=100,
                              cov.user=cov(logistic.mh$mcmc))
pairs_admb(admb_mcmc=logistic.hy2,  diag="acf")
dev.copy2pdf(width=width, height=height,file="Plots/logistic_hy2.pdf")

## Run more of the examples as tests
write.table(x=c(1,1,1,1), file='examples/finance/phases.dat', row.names=FALSE,
            col.names=FALSE)
finance1 <- run_admb_mcmc('examples/finance', 'finance', Nout=1000, mcsave=10,
                          burn.in=5)
pairs_admb(finance1)
dev.copy2pdf(width=width, height=height,file="Plots/finance.pdf")

chem_eng <- run_admb_mcmc('examples/chem-eng', 'chem-eng', Nout=500,
                          mcsave=1000, burn.in=1, mcscale=TRUE)
cov.user <- cov(chem_eng$mcmc)
pairs_admb(chem_eng)
dev.copy2pdf(width=width, height=height,file="Plots/chem_eng.pdf")
pairs_admb(chem_eng, "trace")
cov.user <- cov(chem_eng$mcmc)
chem_eng2 <- run_admb_mcmc('examples/chem-eng', 'chem-eng', Nout=500,
                          mcsave=1000, burn.in=1, mcscale=TRUE,
                           cov.user=cov.user)
pairs_admb(chem_eng2)
chem_eng3 <- run_admb_mcmc('examples/chem-eng', 'chem-eng', Nout=500,
                          mcsave=1000, burn.in=1, mcscale=TRUE, mcrb=2)
pairs_admb(chem_eng3)
pairs_admb(chem_eng3, 'trace')

chem_eng_hy <- run_admb_mcmc('examples/chem-eng', 'chem-eng', Nout=500,
                          mcsave=1, burn.in=1, hybrid=TRUE, hyeps=50,
                           hynstep=100, mcscale=TRUE, cov.user=cov.user)
pairs_admb(chem_eng_hy)
pairs_admb(chem_eng_hy, 'trace')

## End of file
### ------------------------------------------------------------
