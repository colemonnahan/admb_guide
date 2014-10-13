## NOTE: I copied this from another folder and didn't change any paths so
## it won't run without modification. Used to create the image of
## mcprobe. CCM 10/13/2014.
stop("do not source this file")

## Run ADMB mcmcs with the mcgrope and plot the results
setwd("simple_mcgrope")
name <- "MCMCreport.csv"

system("simple")

## Run it once without mcgrope
if(file.exists(name)) file.remove(name)
system("simple -noest -nohess -mcmc 1000000 -mcsave 1 -mcnoscale -mcseed 100")
system("simple.exe -mceval")
mcmc <- read.csv(name, header=F)[-c(1,2),2]
N <- length(mcmc)
## split chain by what was propsoed and what was accepted
results.nogrope <- cbind(mcmc[1:(N/2)], mcmc[(N/2+1):N])
## throw away the burnin
results.nogrope <- results.nogrope[-(1:500),]
results.nogrope <- cbind(results.nogrope, results.nogrope[,1]==results.nogrope[,2])
ratio.nogrope <- mean(results.nogrope[,3])
## Calculate what is proposed
proposed.nogrope <- results.nogrope[-1,1]-results.nogrope[-nrow(results.nogrope),2]
## Run it once with mcgrope
if(file.exists(name)) file.remove(name)
system("simple -noest -nohess -mcmc 1000000 -mcsave 1 -mcnoscale -mcgrope  -mcseed 100")
system("simple.exe -mceval")
mcmc <- read.csv(name, header=F)[-c(1,2),2]
N <- length(mcmc)
## split chain by what was propsoed and what was accepted
results.grope <- cbind(mcmc[1:(N/2)], mcmc[(N/2+1):N])
## throw away the burnin
results.grope <- results.grope[-(1:500),]
results.grope <- cbind(results.grope, results.grope[,1]==results.grope[,2])
ratio.grope <- mean(results.grope[,3])
## Calculate what is proposed
proposed.grope <- results.grope[-1,1]-results.grope[-nrow(results.grope),2]
## Run it again with mcgrope
if(file.exists(name)) file.remove(name)
system("simple -noest -nohess -mcmc 1000000 -mcsave 1 -mcnoscale -mcgrope .49 -mcseed 100")
system("simple.exe -mceval")
mcmc <- read.csv(name, header=F)[-c(1,2),2]
N <- length(mcmc)
## split chain by what was propsoed and what was accepted
results.grope2 <- cbind(mcmc[1:(N/2)], mcmc[(N/2+1):N])
## throw away the burnin
results.grope2 <- results.grope2[-(1:500),]
results.grope2 <- cbind(results.grope2, results.grope2[,1]==results.grope2[,2])
ratio.grope2 <- mean(results.grope2[,3])
## Calculate what is proposed
proposed.grope2 <- results.grope2[-1,1]-results.grope2[-nrow(results.grope2),2]

pow <- 1
temp1 <- (cbind(proposed.nogrope, proposed.grope, proposed.grope2)<0)*-2+1
temp2 <- cbind(abs(proposed.nogrope)^pow, abs(proposed.grope)^pow, abs(proposed.grope2)^pow)
temp <- temp1*temp2


nogrope.col <- gray(.3)
grope.col <- "red"
grope2.col <- "blue"
png("../mcprobe_example.png", width=8, height=4, units="in", res=300)
par(mfrow=c(1,2), cex.axis=.8, mar=c(2,2,2,1), oma=c(0,2,2,0),
    mgp=c(1.5,.2,0), tck=-.02)
boxplot(temp,
        col=c(nogrope.col, grope.col, grope2.col),
        pch=16, cex=.35, names=c("Standard", "-mcprobe", "-mcprobe .49"),
        main="Boxplots")
mtext("Proposed Parameter Values", side=2, line=1.25)
qqnorm(proposed.nogrope, col=nogrope.col, cex=.35, pch=16, ylim=c(-5,5),
       ylab=NA, xlab=NA)
qqline(proposed.nogrope)
points(qqnorm(proposed.grope, plot.it=F), pch=16, col=grope.col, cex=.35)
points(qqnorm(proposed.grope2, plot.it=F), pch=16, col=grope2.col, cex=.35)
legend("topleft", legend=c("Standard", "-mcprobe", "-mcprobe .49"), col=c(nogrope.col,
                                                    grope.col, grope2.col), pch=16, bty="n")
mtext(paste0(prettyNum(N/2, scientific=F,big.mark=","), "Proposed Values"), outer=T, line=.5,
      cex=1.25)
dev.off()


