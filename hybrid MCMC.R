
## Plot a simple example of Hamiltonian dynamics
qt <- function(t, r=1, a=0) r*cos(a+t)
pt <- function(t, r=1, a=0) -r*sin(a+t)
H <- function(t, r, a)
    qt(t=t, r=r, a=a)^2/2 + pt(t=t, r=r, a=a)^2/2
t.seq <- seq(pi/2,10, len=1000)
plot(t.seq, q.seq <- qt(t.seq, r=-1), type="l")
lines(t.seq, p.seq <- pt(t.seq, r=-1), col="red")
plot(0,0, xlim=c(-1,1), ylim=c(-1,1), type='n')
## transformation is a rotation
arrows(0,0, qt(0), pt(0))
arrows(0,0, qt(.5), pt(.5))
## H is conserved
H(0,1,0)
H(1,1,0)

p <- c(2,5)
M <- matrix(c(2,0,0,2.4), nrow=2)
## Looks like a typo in 2.23, 2 should be a 1
t(p) %*% solve(M) %*% p
p^2 %*% (1/diag(M))


## Reproduce figure 1; numerical solutions
## Use Euler's method
par(mfrow=c(2,2))
Nsim <- 21
eps <- .3
psim <- qsim <- rep(NA, Nsim)
qsim[1] <- 0
psim[1] <- 1
for(n in 2:Nsim){
    psim[n] <- psim[n-1] - eps*qsim[n-1]
    qsim[n] <- qsim[n-1] + eps*psim[n-1]
}
plot(q.seq, p.seq, type='l', lwd=2, ylim=c(-2,2), xlim=c(-2,2), col=gray(.5))
lines(qsim, psim, type='b', pch=16)
## Euler's modified method
psim <- qsim <- rep(NA, Nsim)
qsim[1] <- 0
psim[1] <- 1
for(n in 2:Nsim){
    psim[n] <- psim[n-1] - eps*qsim[n-1]
    qsim[n] <- qsim[n-1] + eps*psim[n]
}
plot(q.seq, p.seq, type='l', lwd=2, ylim=c(-2,2), xlim=c(-2,2), col=gray(.5))
lines(qsim, psim, type='b', pch=16)
## The "leapfrog" method
psim <- qsim <- rep(NA, Nsim)
qsim[1] <- 0
psim[1] <- 1
for(n in 2:Nsim){
    ptemp <- psim[n-1] - (eps/2)*qsim[n-1]
    qsim[n] <- qsim[n-1] + eps*ptemp
    psim[n] <- ptemp - (eps/2)*qsim[n]
}
plot(q.seq, p.seq, type='l', lwd=2, ylim=c(-2,2), xlim=c(-2,2), col=gray(.5))
lines(qsim, psim, type='b', pch=16)
## The "leapfrog" method with bigger eps
eps <- 1.2
psim <- qsim <- rep(NA, Nsim)
qsim[1] <- 0
psim[1] <- 1
for(n in 2:Nsim){
    ptemp <- psim[n-1] - (eps/2)*qsim[n-1]
    qsim[n] <- qsim[n-1] + eps*ptemp
    psim[n] <- ptemp - (eps/2)*qsim[n]
}
plot(q.seq, p.seq, type='l', lwd=2, ylim=c(-2,2), xlim=c(-2,2), col=gray(.5))
lines(qsim, psim, type='b', pch=16)


## The HMC function from the paper
HMC <- function(U, grad.U, eps, L, current.q){
    q <- current.q
    p <- rnorm(length(q), 0,1)
    current.p <- p
    ## Make a half step
    p <- p-eps*grad.U(q)/2
    ## alternate full steps for position and momentum
    for (i in 1:L){
        q <- q+ eps*p
        if(i!=L) p <- p-eps*grad.U(q)
    }
    ## half step for momentum at the end
    p <- p-eps*grad.U(q)/2
    ## negate p to make proposal symmetric
    p <- -p
    current.U <- U(current.q)
    current.K <- sum(current.p^2)/2
    proposed.U <- U(q)
    proposed.K <- sum(p^2)/2
    ## Return a list of the current state, proposed, and acceptance
    NLL <- current.U
    ## accept or reject
    if(runif(1) < exp(current.U-proposed.U+current.K-proposed.K))
        final.q <- q
    else final.q <- current.q
    results <- list(current.p=current.p, current.q=current.q, proposed.q=q,
                    q=final.q, NLL=NLL)
    return(results)
}

## Try to get a standard normal working
U <- function(q) q^2/2
grad.U <- function(q) q
Nsim <- 1000
qsim <- rep(NA, len=Nsim)
qsim[1] <- 0                            # initialize at 0
for(n in 2:Nsim){
    qsim[n] <-
        HMC(U, grad.U, eps=.1, L=5, current.q=qsim[n-1])
}
plot(qsim, type='l', ylim=c(-3,3))

## Try a bivariate normal
cor <- matrix(c(1,.5,.5,1), nrow=2)
cor.inv <- solve(cor)
U <- function(q)  t(q) %*% cor.inv %*% q
## Do I suck at math? This probably has a clean matrix formula.
grad.U <- function(q) c(2*q[1]*cor.inv[1,1]+2*q[2]*cor.inv[1,2],
                        2*q[2]*cor.inv[2,2]+2*q[1]*cor.inv[2,1])
Nsim <- 100
qsim <- matrix(NA, nrow=Nsim, ncol=nrow(cor) )
qsim[1,] <- c(-15, -15)
temp <- list()
for(n in 2:Nsim){
    temp[[n]] <- HMC(U, grad.U, eps=.25, L=5, current.q=qsim[n-1,])
    qsim[n,] <- temp[[n]]$q
}
## clean up and organize the results
q <- do.call(rbind, lapply(temp, function(x) x$q))
proposed.q <- do.call(rbind, lapply(temp, function(x) x$proposed.q))
current.q <- do.call(rbind, lapply(temp, function(x) x$current.q))
current.p <- do.call(rbind, lapply(temp, function(x) x$current.p))
NLL <- do.call(rbind, lapply(temp, function(x) x$NLL))
accepted <- proposed.q[,1]==q[,1]

plot(q, type='p', ylim=c(-3,3), xlim=c(-3,3))
plot(proposed.q, type='p', ylim=c(-3,3), xlim=c(-3,3), col=ifelse(accepted,
                                                       "black", "red"))
plot(current.q, type='p', ylim=c(-3,3), xlim=c(-3,3))
plot(qsim)

## Recreate the simple example in R for producing example plots
                                        # observed Y values
a.mle <- 1.90909098475
b.mle <- 4.07817738582

U <- function(q) log(sum((q[1]*xobs+q[2]-yobs)^2))*length(yobs)/2
grad.U <- function(q){
    a <- q[1]; b <- q[2]
    ypred <- a*xobs+b
    dadf <- length(yobs)/2/sum((ypred-yobs)^2)*sum(2*a*(ypred-yobs))
    dbdf <- length(yobs)/2/sum((ypred-yobs)^2)*sum(2*(ypred-yobs))
    return(c(dadf, dbdf))
}
HMC.simple <- function(hyeps, hynstep, current.q){
    ## Define data and the energy functions
    yobs <- c(1.4, 4.7, 5.1, 8.3, 9.0, 14.5, 14.0, 13.4, 19.2, 18)
    xobs <- c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8)
    q <- current.q
    p <- rnorm(length(q), 0,1)
    current.p <- p
    ## Make a half step
    p <- p-hyeps*grad.U(q)/2
    ## alternate full steps for position and momentum
    leapfrog <- matrix(NA, nrow=hynstep+1, ncol=3)
    leapfrog[1,] <- c(q, U(q))
    for (i in 1:hynstep){
        q <- q+ hyeps*p
        leapfrog[i+1,] <- c(q, U(q))
        if(i!=hynstep) p <- p-hyeps*grad.U(q)
    }
    ## half step for momentum at the end
    p <- p-hyeps*grad.U(q)/2
    ## negate p to make proposal symmetric
    p <- -p
    current.U <- U(current.q)
    current.K <- sum(current.p^2)/2
    proposed.U <- U(q)
    proposed.K <- sum(p^2)/2
    ## Return a list of the current state, proposed, and acceptance
    NLL <- current.U
    ## accept or reject
    if(runif(1) < exp(current.U-proposed.U+current.K-proposed.K))
        final.q <- q
    else final.q <- current.q
    results <- list(leapfrog=leapfrog, proposed.q=q, q=final.q, NLL=NLL)
    return(results)
}

## Make a grid of the likelihood surface
a.seq <- seq(1,5, len=50)
b.seq <- seq(2,6, len=50)
param.grid <- expand.grid(a=a.seq, b=b.seq)
NLL <- sapply(1:nrow(param.grid), function(x)
                         U(c(param.grid$a[x], param.grid$b[x])))
NLL <- matrix(NLL, nrow=50)
contour(a.seq, b.seq, NLL)
points(a.mle, b.mle, col='red', pch=16)

simple.hy1 <- HMC.simple(.1, 10, c(2,4))
plot(simple.hy1$leapfrog)
plot(simple.hy1$NLL)

Nsim <- 10
qsim <- matrix(NA, nrow=Nsim, ncol=nrow(cor) )
qsim[1,] <- c(2, 4)
temp <- list()
for(n in 2:Nsim){
    temp[[n]] <- HMC(U, grad.U, eps=.25, L=5, current.q=qsim[n-1,])
    qsim[n,] <- temp[[n]]$q
}
## clean up and organize the results
q <- do.call(rbind, lapply(temp, function(x) x$q))
proposed.q <- do.call(rbind, lapply(temp, function(x) x$proposed.q))
current.q <- do.call(rbind, lapply(temp, function(x) x$current.q))
current.p <- do.call(rbind, lapply(temp, function(x) x$current.p))
NLL <- do.call(rbind, lapply(temp, function(x) x$NLL))
accepted <- proposed.q[,1]==q[,1]

plot(q, type='p', ylim=c(-3,3), xlim=c(-3,3))
plot(proposed.q, type='p', ylim=c(-3,3), xlim=c(-3,3), col=ifelse(accepted,
                                                       "black", "red"))
plot(current.q, type='p', ylim=c(-3,3), xlim=c(-3,3))
plot(qsim)
