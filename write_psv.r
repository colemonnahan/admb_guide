


model.path <- 'mcr_test/simple_long'
model.name <- 'simple'
Nout <- 100
mcsave <- 100
burn.in <- 1

mcr

 cov.user=NULL; init.pin=NULL; se.scale=NULL;
                          mcscale=FALSE;  mcseed=NULL; mcrb=NULL; mcdiag=FALSE;
                          mcprobe=NULL; verbose=TRUE; extra.args=NULL;
                          hybrid=FALSE; hyeps=NULL; hynstep=NULL
########################################################
#Function
## This function runs an ADMB model MCMC, burns and thins, calculates
  ## effective sizes, and returns stuff depending on verbose.
  ## browser()
  iterations <- (Nout+burn.in)*mcsave

  wd.old <- getwd(); on.exit(setwd(wd.old))
  setwd(model.path)
  ## Run to get MLE and covariance matrix
  system(paste0('./', model.name), ignore.stdout=T)
  ## Grab original admb fit and metrics
  mle <- read_admb(model.name)
  ## If user provided covar matrix, write it to file and save to results
  if(!is.null(cov.user)){
    cor.user <- cov.user/ sqrt(diag(cov.user) %o% diag(cov.user))
    if(!is.positive.definite(x=cor.user))
      stop("Invalid cov.user matrix, not positive definite")
    write.admb.cov(cov.user)
    mle$cov.user <- cov.user
  } else {
    ## otherwise use the estimated one
    mle$cov.user <-  NULL
  }
  ## Write the starting values to file. Always using a init.pin file b/c
  ## need to use -nohess -noest so that the cov.user can be specified and
  ## not overwritten. HOwever, this feature then starts the mcmc chain
  ## from the initial values instead of the MLEs. So let the user specify
  ## the init values, or specify the MLEs manually
  if(is.null(init.pin)) init.pin <- mle$coefficients[1:mle$npar]
  write.table(file="init.pin", x=init.pin, row.names=F, col.names=F)
  ## Separate the options by algorithm, first doing the shared arguments
  cmd <- paste(paste0('./', model.name), "-mcmc",iterations)
  ## If user written one, make sure not to overwrite it
  if(!is.null(cov.user)) cmd <- paste(cmd, "-nohess")
  cmd <- paste(cmd, "-mcpin init.pin")
  
  if(!is.null(extra.args)) cmd <- paste(cmd, extra.args)
  if(!is.null(mcseed)) cmd <- paste(cmd, "-mcseed", mcseed)
  if(mcdiag==TRUE) cmd <- paste(cmd, "-mcdiag")
  if(!is.null(mcrb)) cmd <- paste(cmd, "-mcrb",mcrb)
  ## Those options for the standard MH algorithm
  if(!hybrid){
    cmd <- paste(cmd, "-mcsave",mcsave)
    if(mcscale==FALSE) cmd <- paste(cmd, "-mcnoscale")
    if(!is.null(mcprobe)) cmd <- paste(cmd, "-mcprobe",mcprobe)
  } else {
    ## The hybrid options
    if(mcsave!=1)
      stop("mcsave option is incompatible with the hybrid algorithm (fixed at 1 internally)")
    if(hyeps<=0) stop("hyeps must be positive number")
    if(hynstep<=0) stop("hynstep must be positive integer")
    cmd <- paste(cmd, "-hybrid -hyeps",hyeps,"-hynstep",hynstep)
  }
  ## The command is constructed.
  if(verbose) print(cmd)
  ## Scale the covariance matrix
  if(!is.null(se.scale)) SetScale(se.scale) # change the covariance matrix
  ## Run it and get results
  system(cmd, ignore.stdout=T)
  system(paste(model.name, "-mceval -noest -nohess"), ignore.stdout=T)
  psv <- file(paste0(model.name, ".psv"), "rb")
  nparams <- readBin(psv, "integer", n=1)
  mcmc <- matrix(readBin(psv, "numeric", n=nparams*(Nout+burn.in)), ncol=nparams,
                 byrow=TRUE)
  close(psv)
  

  mcmc <- as.data.frame(mcmc)
  

  if(!is.null(mle)){
    names(mcmc) <- names(with(mle, coefficients[1:npar]))
  }
  ## If mcrb is used, read  that in for plotting. It's in the corrtest file.
  if(!is.null(mcrb)){
    L <- readLines('corrtest')
    st <- grep("modified S", L)[1]+1
    en <- grep("S* modified S", L)-1
    L <- gsub("^\\s+|\\s+$", "", L[st:en])
    cov.mcrb <- do.call(rbind, lapply(strsplit(L, " "), as.numeric))
    cor.mcrb <- cov.mcrb/(sqrt(diag(cov.mcrb) %o% diag(cov.mcrb)))
    if(!is.positive.definite(cor.mcrb))
      warning("the modified mcrb covariance matrix read in was not positive definite")
    else mle$cov.user <- cov.mcrb
  }
  ## Remove the 'burn in' specified by user
  if(burn.in>0) mcmc <- mcmc[-(1:burn.in),]
  ## Run effective sample size calcs from CODA, metric of convergence
  efsize <- data.frame(t(effectiveSize(mcmc)/NROW(mcmc)))
  names(efsize) <- paste0(names(mcmc), "_efs")
  results <- list(mcmc=mcmc, mle=mle)
  results$diag <- list(efsize=efsize)
  class(results) <- 'admb_mcmc'
  



