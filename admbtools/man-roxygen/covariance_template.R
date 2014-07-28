#' @details Read and write ADMB \code{admodel.hes} and \code{admodel.cov}
#' files. This function reads in all of the information contained in the
#' admodel.hes file. It returns a list w/ the number of parameters, the
#' Hessian matrix (in unbounded space), the hybrid bounded flag (determines
#' the bounding function), and the scales, which are the second derivatives
#' of the bounding functions evaluated at the MLE (used to convert from
#' unbounded to bounded space).
#'
#' @param model.path The directory path to the folder containing the model.
#' @param cov.user (Numeric matrix) A covariance matrix (in bounded space)
#' to write to file. This matrix will be read in by ADMB and can be used in
#' subsequent calculations.
#' @author Cole Monnahan
#' @details The admodel.hes file is written by ADMB during the optimization
#' phase. The Hessian matrix contains the second derivatives evaluated at
#' the maximum likelihood estimates. The inverse of the Hessian is the
#' covariance matrix used to determine standard errors and in the MCMC
#' algorithm. This file is written even if the matrix is not invertible
#' (not positive definite).
#' @return The \code{write} functions return nothing, but have the side
#' effect of changing ADMB files (including creating a backup of the
#' original). The \code{get} functions return a list of the number of
#' parameters (\code{num.pars}), the covariance or Hessian matrix in
#' unbounded space (\code{cov} or \code{hes}), the hybrid bounded flag
#' (\code{hybrid_bounded_flag}) which controls which bounding functions are
#' used, and the \code{scales}, which are the second derivatives of the
#' bounding functions evaluated at the MLE.

