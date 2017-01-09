frailty_function <- function (p0, delta, N) {
  num.fully.protected <- floor(p0*N)
  susceptibility <- c(rep(0, num.fully.protected), 
                      rgamma(N-num.fully.protected, 1/delta, 1/delta))
  return(susceptibility)
}

simulate.data <- function (
  N, # even number of participants. Split between the two arms 50:50
  ntypes, # number of serotypes in the vaccine
  ntot, # total number of serotypes including those not in vaccines
  nswabs, # number of times nasopharyngeal swabs were taken
  times, # time intervals between swabs (length=nswabs)
  lambda, # rate of acquisition (length=ntot)
  mu, # rate of clearance (length=ntot)
  p0, # probability of full protection (length=ntypes)
  thetaSI, # relative rate of acquisition (length=ntypes)
  thetaIS, # relative rate of clearance (length=ntypes)
  frailtySI, # frailty parameter that determines individual-level heterogeneity in preventing acquisition
  frailtyIS, # frailty parameter that determines individual-level heterogeneity in clearance rates
  interaction # rate at which the current  
) {
  ### For testing within the function
  # load("simulate.data.test.RData")
  ###
  require(parallel)
  require(Matrix)
  patient.data <- data.frame(ID=1:N, vacc=rep(c(TRUE, FALSE), each=N/2))
  # The rates of clearance or infection by another serotype
  rates.matrix <- replicate(ntot, lambda*interaction)*(1-diag(ntot))
  rates.matrix <- rbind(c(-sum(lambda), lambda),
                        cbind(mu, rates.matrix))
  # Individual heterogeneities in rates of acquisition(frailtySI)/clearance(frailtyIS) after being vaccinated
  frailtySI.mat <- mapply(frailty_function, p0, frailtySI, N/2)
  frailtyIS.mat <- mapply(frailty_function, 0, frailtyIS, N/2)
  ###
  ### Check that the frailty matrices have been simulated correctly
  # require(ggplot2); ggplot(reshape2::melt(frailtyIS.mat)) + geom_tile(aes(x=Var2, y=max(Var1)-Var1+1, colour=value))
  ###
  # Simulate data for each patient
  patient.data <- cbind(patient.data, do.call(rbind, lapply(1:N, function (patient_i) {
    rates.matrix_i <- rates.matrix
    # Modify rates for vaccinated individuals
    if (patient.data$vacc[patient_i]) {
      rates.matrix_i[2:(ntot+1), 2:(ntypes+1)] <- 
        sweep(rates.matrix[2:(ntot+1), 2:(ntypes+1)], MARGIN=2, thetaSI * frailtySI.mat[patient_i, ], `*`)
      rates.matrix_i[2:(ntypes+1), 1] <- rates.matrix_i[2:(ntypes+1), 1] * frailtyIS.mat[patient_i, ] * thetaIS
    }
    diag(rates.matrix_i)[-1] <- -rowSums(rates.matrix_i[-1, -1])
    # Construct vector for observations
    observations <- rep(0, nswabs)
    multinom.p <- expm(rates.matrix_i*1e0)[1, ]
    # Simulate first observation
    observations[1] <- which(rmultinom(N, 1, multinom.p)==1, arr.ind = TRUE)[patient_i, 1]
    P <- expm(rates.matrix_i*times[1])
    for (time.step in 2:nswabs) {
      observations[time.step] <- sample(ntot+1, 1, prob=P[observations[time.step-1], ])
      P <- expm(rates.matrix_i*(times[time.step]-times[time.step-1]))
    }
    return (observations-1)
  })))
  return(patient.data)
}