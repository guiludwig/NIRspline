#' nirSpline.default
#' 
#' Default action for the nirSpline function
#' 
#' @param wavelengths 
#' @param absorbance
#' @param concentration
#' @param norder
#'
#' @export
#' @return List of two items 
#'   \item{betaMatrix}{TODO}
#'   \item{fittedMatrix}{TODO}
#'
#' @examples
#' require(fda)
#' data(tecator)
#' with(tecator, nirSpline(seq(850, 1050, length = 100), training[,1:100], 
#'                         training[,123:125]/100))
#' 
#' @references
#'  Ramsay, J. O. (2006) Functional Data Analysis. New York: Springer.
#'  
#'  Dias, R., Garcia, N., Ludwig, G. and Saraiva, M. (2014) Aggregated functional 
#'  data model for Near-Infrared Spectroscopy calibration and prediction, 
#'  Journal of Applied Statistics
#'
#' @seealso \code{\link{tecator}}
#' @keywords Aggregated Functional Data Analysis
nirSpline.default <- function(wavelengths, absorbance, concentration, norder = 4){
  I <- dim(absorbance)[1] # Number of components
  n <- dim(absorbance)[2] # Number of wavelength samples
  C <- dim(concentration)[2] # Number of analytes
  if(dim(concentration)[1] != I) stop("Provide complete concentration values 
                                      for all absorbance curves.")
  K <- n-norder
  absbasis <- fda::create.bspline.basis(rangeval = range(wavelengths), 
                                        norder = norder, nbasis = K)
  B <- fda::getbasismatrix(wavelengths, basisobj = absbasis)
  X <- c(as.numeric(t(absorbance)), rep(0,n))    
  Y <- cbind(rep(1,I), concentration)
  Y <- rbind(Y, c(0, rep(1,C)))
  R <- fda::getbasispenalty(absbasis)
  
  D.s <- kronecker(Y,B)
  DD.s <- crossprod(D.s,D.s) # efficient computation
  crossy <- crossprod(D.s,X)
  p <- dim(Y)[2]
  N <- length(X)
  
  GCV.2 <- function(lambda, y = X, ycross = crossy, DD = DD.s, 
                    D = D.s, P = p, RR = R){
    H.t <- solve(DD + kronecker(diag(P), lambda*RR), DD) # efficient trace
    trace <- sum(diag(H.t))
    adjust <- D%*%solve(DD + kronecker(diag(P),lambda*RR), ycross)
    gcv <- N*sum((y-adjust)^2)/(length(y)-trace)
    return(gcv)
  }
  lambda <- optimize(GCV.2, c(10^(-5),10^5))$minimum
  
  coefficients <- solve(DD.s + kronecker(diag(p),lambda*R), crossprod(D.s, X))
  betaMatrix <- matrix(coefficients, nrow = K)
    for(c in 1:C){
    betaMatrix[, c+1] <- betaMatrix[, c+1] + betaMatrix[, 1]
  }
  
  fitted <- D.s%*%coefficients
  fitted <- as.numeric(fitted[c(1:n-1-length(X))]) # remove zeros
  fittedMatrix <- matrix(fitted, ncol = n, byrow = TRUE)
  
  ret <- list(betaMatrix = betaMatrix, fittedMatrix = fittedMatrix)
  class(ret) <- "nirSpline"
  
  return(ret)
}