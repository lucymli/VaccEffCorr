frailty_function <- function (p0, delta, N) {
  num.fully.protected <- floor(p0*N)
  susceptibility <- c(rep(0, num.fully.protected), 
                      rgamma(N-num.fully.protected, 1/delta, 1/delta))
  return(susceptibility)
}