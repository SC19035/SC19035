#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C functions.
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats runif rnorm
#' @useDynLib SC19035
#' @examples
#' \dontrun{
#' ts <- microbenchmark(
#' rwC = rwMetropolisC(0.5, x0, N=1000),
#' rwR = rw.Metropolis(0.5, x0, N=1000))
#' print(summary(ts)[, c(1,3,5,6)])
#' 
#' }
NULL


#' @title A rw.Metropolis sampler using R
#' @description A rw.Metropolis sampler using R
#' @param N the number of samples
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @param x0 the random numbers
#' @param sigma the sigma
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' sigma<-c(0.1,0.2)
#' x0<-200
#' rw.Metropolis(sigma[1], x0, N=1000)
#' }
#' @export
rw.Metropolis <- function(sigma, x0, N=1000) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= exp(-((abs(y)) - (abs(x[i-1])))))
      x[i] <- y
    else {
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(list(x = x, k = k))
}
