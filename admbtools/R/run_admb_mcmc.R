#' Run an MCMC using an ADMB model, return (1) the posterior draws, MLE
#' fits and covariance/correlation matrices, and some MCMC convergence
#' diagnostics using CODA.
#'
#' @param model.path (Character) A path to the folder containing the model. NULL
#' indicates the current folder.
#' @param mode.name (Character) The name of the model executable. A character string,
#' without '.exe'.
#' @param Nout (Integer) The number of draws after thinning and burn in.
#' @param mcsave (Integer) Controls thinning of samples. Save every mcsave
#' value, such that 1 corresponds to keeping all draws, and 100 saving
#' every 100th draw.
#' @param burn.in (Integer) How many samples to discard from the beginning
#' of the chain, *after* thining. The burn in period (i.e., the first
#' burn.in*mcsave draws) should be at least large enough to cover dynamic
#' scaling.
#' @param cor.user (Numeric matrix) A manually defined correlation matrix (in bounded space)
#' to use in the Metropolis-Hastings algorithm.
#' @param init.pin (Numeric vector) A vector of initial values, which are written to file
#' and used in the model via the -mcpin option.
#' @param se.scale (Numeric) A value which scales all of the variances from
#' the MLE fit. A value of 1 indicates to use the estimated variances.
#' @param mcscale (Logical) Whether to use the mcscale option, which
#' dynamically scales the covariance matrix for efficient acceptance
#' ratios.
#' @param mcseed (Integer) Which seed (integer value) to pass ADMB. Used
#' for reproducibility.
#' @param mcrb (Integer) Which value to use in the rescale bounded
#' algorithm. Must be an integer from 1-9. The default NULL value disables
#' this feature. See the vignette for more information on this algorithm
#' and how to best use it.
#' @param mcdiag (Logical) Whether to use the \code{mcdiag} feature. This
#' uses an identity matrix for the covariance matrix.  #' @param mcprobe
#' Which value to use in the probing algorithm. The default NULL value
#' disables this feature. See the vignette for more information on this
#' algorithm and how to best use it.
#' @param verbose (Logical) Whether to print ADMB warnings and other
#' information. Useful for testing and troubleshooting.
#' @param extra.args (Character) A string which is passed to ADMB at
#' runtime. Useful for passing additional arguments to the model
#' executable.
#' @return Returns a list containing (1) the posterior draws, (2) and
#' object of class 'admb', read in using the results read in using
#' \code{read_admb}, and (3) some MCMC convergence diagnostics using CODA.
run_admb_mcmc <- function(model.path, model.name, Nout, mcsave, burn.in,
                     cor.user=NULL, init.pin=NULL, se.scale=NULL,
                     mcscale=FALSE,  mcseed=NULL, mcrb=NULL, mcdiag=FALSE,
                     mcprobe=NULL, verbose=TRUE, extra.args=NULL){
    ## This function runs an ADMB model MCMC, burns and thins, calculates
    ## effective sizes, and returns stuff depending on verbose.
    ## browser()
    iterations <- (Nout+burn.in)*mcsave
    if(verbose) print(paste("Run started at", round(Sys.time())))
    if(iterations <1) stop(paste0("Iterations too low: ", iterations))
    if(!is.null(model.path)) {
        wd.old <- getwd(); on.exit(setwd(wd.old))
        setwd(model.path)
    }
    ## Run to get MLE and covariance matrix
    system(model.name, ignore.stdout=T)
    ## Grab original admb fit and metrics
    mle <- read_admb(model.name)
    ## If user provided covar matrix, write it to file and save to results
    if(!is.null(cor.user)){
        if(!is.positive.definite(cor.user))
            stop("Invalid cor.user matrix, not positive definite")
        write.admb.cov(cor.user)
        mle$cor.user <- cor.user
        cor <- cor.user
    }
    ## otherwise use the estimated one
    else {
        cor <- with(mle, cor[1:npar, 1:npar])
        mle$cor.user <-  NULL
    }
    ## If init values not provided, start from the MLEs (see below for why
    ## necessary)
    if(is.null(init.pin)) init.pin <- mle$coefficients[1:mle$npar]
    ## Write the starting values to file
    write.table(file="init.pin", x=init.pin, row.names=F, col.names=F)
    ## make the argument and run the MCMC
    cmd <- paste(model.name,"-mcmc",iterations, "-nohess -noest")
    if(!is.null(mcseed)) cmd <- paste(cmd, "-mcseed", mcseed)
    cmd <- paste(cmd, "-mcsave",mcsave)
    if(mcscale==FALSE) cmd <- paste(cmd, "-mcnoscale")
    if(!is.null(mcrb)) cmd <- paste(cmd, "-mcrb",mcrb)
    if(!is.null(mcprobe)) cmd <- paste(cmd, "-mcprobe",mcprobe)
    if(mcdiag==TRUE) cmd <- paste(cmd, "-mcdiag")
    if(!is.null(extra.args)) cmd <- paste(cmd, extra.args)
    ## Always using a init.pin file b/c need to use -nohess -noest so that
    ## the cor.user can be specified and not overwritten. HOwever, this
    ## feature then starts the mcmc chain from the initial values instead
    ## of the MLEs. So let the user specify the init values, or specify the
    ## MLEs manually
    cmd <- paste(cmd, "-mcpin init.pin")
    if(verbose) print(cmd)
    ## Scale the covariance matrix
    if(!is.null(se.scale)) SetScale(se.scale) # change the covariance matrix
    system(cmd, ignore.stdout=T)
    system(paste(model.name, "-mceval -noest -nohess"), ignore.stdout=T)
    ## Read in the MCMC results
    psv <- file(paste0(model.name, ".psv"), "rb")
    nparams <- readBin(psv, "integer", n=1)
    mcmc <- matrix(readBin(psv, "numeric", n=nparams*(Nout+burn.in)), ncol=nparams,
                   byrow=TRUE)
    close(psv)
    mcmc <- as.data.frame(mcmc)
    if(!is.null(mle)){
        names(mcmc) <- names(with(mle, coefficients[1:npar]))
    }
    ## Remove the 'burn in' specified by user
    if(burn.in>0) mcmc <- mcmc[-(1:burn.in),]
    ## Run effective sample size calcs from CODA, metric of convergence
    efsize <- data.frame(t(effectiveSize(mcmc)/NROW(mcmc)))
    names(efsize) <- paste0(names(mcmc), "_efs")
    results <- list(mcmc=mcmc, mle=mle)
    results$diag <- list(efsize=efsize)
    class(results) <- 'admb_mcmc'
    return(results)
}
