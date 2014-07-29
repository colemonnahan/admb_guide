#' The bounding functions, with their inverses and derivatives.
#'
#' @details If \code{hbf=0} then the function is minb+(maxb-minb)/(1+exp(-x)),
#' otherwise it is minb+(maxb-minb)*(.5*sin(x*pi/2)+.5).
#'
#' @param x The input value.
#' @param minb The minimum bound.
#' @param maxb The maximum bound.
#' @param hbf The hybrid bounded flag. Determines which function is used.
#' @return Numeric function value.
#' @author Cole C. Monnahan
#' @section Note: These functions are named as in the ADMB source code.
