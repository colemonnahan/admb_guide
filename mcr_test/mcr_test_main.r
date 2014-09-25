# setwd('/Users/peterkuriyama/School/Research/admb_guide')
setwd('/Users/peterkuriyama/School/Research/admb_guide/mcr_test/mcr_check')

## Make sure these are updated
library(devtools)
library(knitr)
library(roxygen2)
library(coda)
library(matrixcalc)
library(ellipse)
library(R2admb)
library(coda)

check_mcr <- function(model.name = 'simple', 
                      Nout,
                      mcsave,
                      burn.in = 0,
                      mcseed, 
                      Nout.mcr,
                      mcsave.mcr){

  ########################################
  #Run 1 - 1/2 Run
  iterations <- (Nout + burn.in) * mcsave

  #Creat Cmd for MAC
  cmd <- paste0('./', model.name)
  cmd <- paste(cmd, '-mcmc', iterations, '-mcsave', mcsave, 
               '-mcseed', 30)
  print(cmd)
  #Run admb command
  system(cmd, ignore.stdout = TRUE)

  #Copy .hst File
  # system('cp -i simple.hst simple1.hst')

  #Extract MCMC draws and save locally
  psv <- file("simple.psv", "rb")
  nparams <- readBin(psv, "integer", n=1)
  mcmc1 <- matrix(readBin(psv, "numeric", n= 2 * (Nout + burn.in)), ncol=2,
                   byrow=TRUE)
  close(psv)

  ########################################
  #Run 2 - 1/2 Run

  # Nout.mcr <- 50
  # mcsave.mcr <- 50

  #MCR Runs
  iterations.mcr <- (Nout.mcr + burn.in) * mcsave.mcr
  print(iterations.mcr)

  cmd <- paste0('./', model.name)
  cmd <- paste(cmd, '-mcmc', iterations.mcr, '-mcsave', mcsave.mcr,
               '-mcr', '-nosdmcmc')

  print(cmd)
  system(cmd, ignore.stdout = TRUE)

  # system('cp -i simple.hst simple2.hst')

  psv <- file("simple.psv", "rb")
  nparams <- readBin(psv, "integer", n=1)

  mcmc2 <- matrix(readBin(psv, "numeric", n = 2*(2*(Nout.mcr + burn.in))), ncol=2,
                   byrow=TRUE)
  close(psv)

  nrow(mcmc2) - nrow(mcmc1)


  mcmc1 <- rbind(mcmc1, matrix(0, nrow = (nrow(mcmc2) - nrow(mcmc1)), 
    ncol = 2))
  colnames(mcmc1) <- c('a-first', 'b-first')
  colnames(mcmc2) <- c('a-resumed', 'b-resumed')
  mcmc.check <- cbind(mcmc1, mcmc2)
  mcmc.check <- round(mcmc.check, digits = 3)
  #Print things
  # print(mcmc.check)
  write.csv(mcmc.check, file = 'mcmc_check.csv')
  # par(mfcol = c(1, 2))
  # acf(mcmc.check[, 1])
  # acf(mcmc.check[, 3])
  return(mcmc.check)
}


out <- check_mcr(Nout = 10, mcsave = 20, mcseed = 30,
          Nout.mcr = 10, mcsave.mcr = 50)


# kable(out)






# model.name <- 'simple'
# Nout <- 50
# mcsave <- 50
# burn.in <- 0
# mcseed <- 30
# Nout.mcr <- 50
# mcsave.mcr <- 50










# ########################################
# #Run 3 - Full Length

# iterations.long <- iterations + iterations.mcr

# cmd <- paste0('./', model.name)
# cmd <- paste(cmd, '-mcmc', iterations.long, '-mcsave', mcsave,
#              '-mcseed', mcseed)

# system(cmd)

# system('cp -i simple.hst simple3.hst')

# psv <- file("simple.psv", "rb")
# nparams <- readBin(psv, "integer", n=1)

# mcmc3 <- matrix(readBin(psv, "numeric", n = 500), ncol=2,
#                  byrow=TRUE)

# #Expand MCMC1
# mcmc1 <- rbind(mcmc1, matrix(0, nrow = (nrow(mcmc2) - nrow(mcmc1)), 
#   ncol = 2))
# colnames(mcmc1) <- c('first half a', 'first half b')

# colnames(mcmc2) <- c('second half a', 'second half b')

# colnames(mcmc3) <- c('full a', 'full b') 

# #Write Final
# check <- cbind(mcmc1, mcmc2, mcmc3)
# write.csv(check, file = 'mcmc_check.csv')



# # # mcmc3 <- matrix(readBin(psv, "numeric", n = (Nout + burn.in + 30) * 2), ncol=2,
# # #                  byrow=TRUE)

# # mcmc3[1:131, ]

# # (mcmc2)


# # close(psv)



# # short <- test_mcr(model.path = 'mcr_test/simple_long',
# #                   model.name = 'simple',
# #                   Nout = 100,
# #                   mcsave = 100,
# #                   burn.in = 1,
# #                   mcseed = 30,
# #                   mcr = NULL)

# # long <- test_mcr(model.path = 'mcr_test/simple_long',
# #                   model.name = 'simple',
# #                   Nout = 100,
# #                   mcsave = 100,
# #                   burn.in = 1,
# #                   mcseed = 30,
# #                   mcr = 1e5)
