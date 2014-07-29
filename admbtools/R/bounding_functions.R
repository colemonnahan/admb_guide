#' The bounding function.
#'
#' @template bounding_template
boundp <- function(x, minb, maxb, hbf=0){
    ## The internal transformations used in ADMB depending on the value of the
    ## Hybrid_bounded_flag (hbf) value.
    if(hbf==1)
        result <- minb+(maxb-minb)/(1+exp(-x))
    else if(hbf==0)
        result <- minb+(maxb-minb)*(.5*sin(x*pi/2)+.5)
    else stop("Invalid hbf value, should be 0 or 1")
    return(result)
}

#' Inverse bounding transformation function used by ADMB.
#'
#' @template bounding_template
boundpin <- function(x, minb, maxb, hbf) {
    ## The inverse of the transformation
    if(hbf==1)
        result <- -log( (maxb-x)/(x-minb) )
    else if(hbf==0)
        result <- asin(2*(x-minb)/(maxb-minb)-1)/(pi/2)
    else stop("Invalid hbf value, should be 0 or 1")
    return(result)
}

#' Derivative of the bounding transformation function used by ADMB.
#'
#' @template bounding_template
ndfboundp <- function(x, minb, maxb, hbf) {
    ## The derivative used to find the "scales"
    if(hbf==1)
        result <- (maxb-minb)*exp(-x)/(1+exp(-x))^2
    else if(hbf==0)
        result <- (maxb-minb)*.5*pi/2*cos(x*pi/2)
    else stop("Invalid hbf value, should be 0 or 1")
    return(result)
}
