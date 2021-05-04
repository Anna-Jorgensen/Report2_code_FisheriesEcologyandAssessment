#' Virtual Population Analysis
#'
#' Run simple VPA model, with no tuning index.
#'
#' @param C catch-at-age matrix.
#' @param M natural mortality rate, either a scalar or a vector of same length as the
#'        number of ages.
#' @param Fterm terminal F, either a scalar or a vector of same length as the
#'        number of ages.
#' @param Fages number of ages to calculate F for second-oldest age. For example, if
#'        the catch-at-age matrix contains ages up to 15 and \code{Fages=5},
#'        then F for the ages 14 and 15 will be set as the average of ages 9-13.
#'
#' @return
#' List containing the model results (\code{N}, \code{F}, \code{Z}), as well as
#' the input objects passed by the user (\code{C}, \code{M}, \code{Fterm},
#' \code{Fages}).
#'
#' @export

vpa = function(C, M, Fterm, Fages)
{
  ## Prepare matrices
  C = as.matrix(C)
  T = nrow(C)
  A = ncol(C)
  N = F = Z = matrix(NA_real_, nrow=T, ncol=A, dimnames=dimnames(C))

  ## Check size of M
  if(length(M)==1) M=rep(M, A)
  
  ## Set F in terminal year
  F[T,] = Fterm
  Z[T,] = F[T,] + M

  ## Calculate N in terminal year
  N = C*Z / (F*(1-exp(-Z))) #Eq 3

  ## Calculate N and F up to terminal year,
  ## assuming F[oldest] = avg(preceding ages)
  for(t in (T-1):1)
  {
    for(a in 1:(A-1))
    {
  	  N[t,a] = N[t+1,a+1] * exp(M[a]) + C[t,a] * exp(M[a]/2) #Eq 2
      F[t,a] = log(N[t,a] / N[t+1,a+1]) - M[a] #Eq 1
    }
    F[t,A] = mean(F[t,A-1-(1:Fages)])
    Z[t,] = F[t,] + M
    N[t,A] = C[t,A]*Z[t,A] / (F[t,A]*(1-exp(-Z[t,A]))) #Eq 3
  }

  list(N=N, F=F, Z=Z, C=C, M=M, Fterm=Fterm, Fages=Fages)
}
