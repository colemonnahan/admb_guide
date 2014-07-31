#' Examine MCMC output from an ADMB model.
#'
#' This function is useful for checking covergence and posterior
#' properties.
#'
#' @param posterior Dataframe containing the MCMC output, as read in using
#' function \link{\code{run.mcmc}}
#' @param diag What type of plot to include on the diagonal, options are
#' 'acf' which plots the autocorrelation function \code{acf}, 'hist' shows
#' marginal posterior histograms, and 'trace' the trace plot.
#' @param acf.ylim If using the acf function on the diagonal, specify the y
#' limit. The default is c(-1,1).
#' @param ymult A vector of length ncol(admb_mcmc) specifying how much room to
#' give when using the hist option for the diagonal. For use if the label
#' is blocking part of the plot.
#' @param limits A list containing the ranges for each parameter to use in
#' plotting.
#' @return Produces a plot, and returns nothing.
#' @author Cole Monnahan
#' @export
pairs_admb <- function(admb_mcmc, diag=c("acf","hist", "trace"),
                       acf.ylim=c(-1,1), ymult=NULL, axis.col=gray(.5),
                       limits=NULL, ...){
    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    diag <- match.arg(diag)
    posterior.names <- names(admb_mcmc$mcmc)
    posterior <- admb_mcmc$mcmc
    mle <- admb_mcmc$mle
    n <- NCOL(posterior)
    if(n==1) warning("This function is only meaningful for >1 parameter")
    if(is.null(ymult)) ymult <- rep(1.3, n)
    ## If no limits given, calculate the max range of the posterior samples and
    ## parameter confidence interval
    if(is.null(limits)){
        limits <- list()
        for(i in 1:n){
            limit.temp <- with(mle, coefficients[i]+c(-1,1)*1.96*se[i])
            ## multiplier for the ranges, adjusts the whitespace around the
            ## plots
            min.temp <- min(posterior[,i], limit.temp[1])
            max.temp <- max(posterior[,i], limit.temp[2])
            margin <- .15*(max.temp-min.temp)
            limits[[i]] <- c(min.temp-margin, max.temp+margin)
        }
    }
    par(mfrow=c(n,n), mar=0*c(.1,.1,.1,.1), yaxs="i", xaxs="i", mgp=c(.25, .25,0),
        tck=-.02, cex.axis=.65, col.axis=axis.col, oma=c(2, 2, 2,.5))
    temp.box <- function() box(col=axis.col, lwd=.5)
    ## Row and col here are not the posterior, but the matrix of pairwise
    ## combinations
    for(row in 1:n){
        for(col in 1:n){
            ## Diagonal, so add user choice
            if(row==col){
                if(diag=="hist"){
                    h <- hist(posterior[,row], plot=F)
                    ## Annoyingling you can't pass NULL to xlim in hist. So
                    ## have to split up for two cases depending on limits.
                    if(is.null(limits)){
                        hist(posterior[,row], axes=F, freq=FALSE, ann=F,
                             ylim=c(0, ymult[row]*max(h$density)),
                             col=gray(.8), border=gray(.5))
                    } else {
                        ## Else use the user provided limits
                        hist(posterior[,row], axes=F, freq=FALSE, ann=F,
                             ylim=c(0, ymult[row]*max(h$density)),
                             col=gray(.8), border=gray(.5), xlim=limits[[row]])
                    }
                    temp.box()
                } else if(diag=="acf") {
                    acf(posterior[,row], axes=F, ann=F, ylim=acf.ylim)
                    ## legend("topright", bty='n', legend=NA,
                    ##        title=sprintf("EFS=%.3f", 100*admb_mcmc$diag$efsize[row],2))
                    temp.box()
                } else if(diag=="trace") {
                    plot(x=posterior[,row], lwd=.5, col=gray(.5), type="l", axes=F,
                         ann=F, ylim=limits[[row]])
                    temp.box()
                }
                mtext(posterior.names[row], line=-2)
            }
            ## If lower triangle add scatterplot
            if(row>col){
                par(xaxs="r", yaxs="r")
                plot(x=posterior[,col], y=posterior[,row], axes=FALSE, ann=FALSE,
                     pch=ifelse(NROW(posterior)>=5000,".", 1), col=1,
                     xlim=limits[[col]], ylim=limits[[row]])
                ## Add bivariate 95% normal levels for both the MLE
                ## estimated covariance, but also the user supplied cov.user
                points(x=mle$coefficients[col], y=mle$coefficients[row],
                       pch=16, cex=1, col=2)
                ## Get points of a bivariate normal 95% confidence contour
                ellipse.temp <- ellipse::ellipse(x=mle$cor[col, row],
                                       scale=mle$se[1:mle$npar][c(col, row)],
                                       centre= mle$coefficients[c(col, row)], npoints=1000,
                                       level=.95)
                lines(ellipse.temp , lwd=1.5, lty=1, col="red")
                if(!is.null(mle$cov.user)){
                    se.user <- sqrt(diag(mle$cov.user))
                    cor.user <- mle$cov.user/(se.user %o% se.user)
                    lines(ellipse::ellipse(
                        x=cor.user[col, row],
                        scale=se.user[c(col, row)],
                        centre= mle$coefficient[c(col, row)], npoints=1000,
                        level=.95) , lwd=1.5, lty=1, col="blue")
                }
                par(xaxs="i", yaxs="i")
                temp.box()
            }
            if(row<col){
                ## If upper triangle add text showing the empirical correlation
                plot(0,0,type="n", xlim=c(0,1), ylim=c(0,1), axes=F,ann=F)
                temp.cor <- round(cor(posterior[,c(row,col)])[1,2],2)
                ## Set a minimum limit for this, so they're still
                ## visible, but still a function of correlation. This
                ## might need to be dynamic with n.
                legend("center", legend=NA, title=temp.cor,
                       cex=(3*abs(temp.cor)+.25)*.5, bty='n')
                ## text(.5,.5, labels=temp.cor, cex=(3*abs(temp.cor)+.5)*.9,
                ##      col=1)
                temp.box()
            }
            ## Add special cases of axes on the ends
            if(row==n) {
                par( mgp=c(.05, ifelse(col %% 2 ==0, 0, .5),0) )
                axis(1, col=axis.col, lwd=.5)
            }
            if(col==1 & row >1) {
                par( mgp=c(.05, ifelse(row %% 2 ==1, .15, .65),0) )
                axis(2, col=axis.col, lwd=.5)
            }
            if(col==1 & row ==1){
                par( mgp=c(.05, ifelse(row %% 2 ==1, .15, .65),0) )
                axis(2, col=axis.col, lwd=.5)
            }

        }
    }
}

