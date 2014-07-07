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
#' @param fits A list containing the fitted information, read in using xxx
#' @param acf.ylim If using the acf function on the diagonal, specify the y
#' limit. The default is c(-1,1).
#' @param ymult A vector of length ncol(posterior) specifying how much room to
#' give when using the hist option for the diagonal. For use if the label
#' is blocking part of the plot.
#' @param limits A list containing the ranges for each parameter to use in
#' plotting.
#' @export
admb.pairs <- function(posterior, diag=c("acf","hist", "trace"), fits=NULL,
                       acf.ylim=c(-1,1), ymult=NULL, axis.col=gray(.5),
                       limits=NULL, ...){
    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    diag <- match.arg(diag)
    posterior.names <- names(posterior)
    posterior <- as.matrix(posterior)
    n <- NCOL(posterior)
    if(n==1) warning("This function is only meaningful for >1 parameter")
    if(is.null(ymult)) ymult <- rep(1.3, n)
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
                    hist(posterior[,row], axes=F, freq=FALSE, ann=F,
                         ylim=c(0, ymult[row]*max(h$density)),
                         col=gray(.8), border=gray(.5), xlim=limits[[row]])
                }
                    temp.box()
                } else if(diag=="acf") {
                    acf(posterior[,row], axes=F, ann=F, ylim=acf.ylim)
                    ## If the fit is  given, add effective sample sizes
                    if(!is.null(fits)){
                        legend("topright", bty='n', legend=NA,
                               title=sprintf("EFS=%.3f", 100*fits$efsize[row],2))                    }
                    temp.box()
                } else if(diag=="trace") {
                    plot(x=posterior[,row], lwd=.5, col=gray(.5), type="l", axes=F,
                         ann=F, ylim=limits[[row]])
                    temp.box()
                }
                legend("top", bty='n', legend=NA, title=posterior.names[row], ...)
            }
            if(row>col){
                ## If lower triangle add scatterplot
                par(xaxs="r", yaxs="r")
                plot(x=posterior[,col], y=posterior[,row], axes=FALSE, ann=FALSE,
                     pch=ifelse(NROW(posterior)>=5000,".", 1), col=1,
                     xlim=limits[[col]], ylim=limits[[row]])
                ## Add bivariate 95% normal levels for both the MLE
                ## estimated covariance, but also the user supplied cor.mat
                if(!is.null(fits)){
                    points(x=fits$est[col], y=fits$est[row],
                           pch=16, cex=1, col=2)
                    lines(ellipse::ellipse(x=fits$cor[col, row],
                                  scale=fits$std[1:fits$nopar],
                                  centre= fits$est[c(col, row)], npoints=1000,
                                  level=.95) , lwd=1.5, lty=1, col="red")
                    if(!is.null(fits$cor.mat))
                        lines(ellipse::ellipse(x=fits$cor.mat[col, row],
                                               scale=fits$std[1:fits$nopar],
                                               centre= fits$est[c(col, row)], npoints=1000,
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
                       cex=(3*abs(temp.cor)+.5)*.9, bty='n')
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
        }
    }
}

