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
  # load("simulate.data.test.RData")
  # source("functions.R")
  require(parallel)
  require(Matrix)
  patient.data <- data.frame(ID=1:N, vacc=rep(c(TRUE, FALSE), each=N/2))
  rates.matrix <- replicate(ntot, lambda*interaction)*(1-diag(ntot))
  rates.matrix <- rbind(c(-sum(lambda), lambda),
                        cbind(mu, rates.matrix))
  frailtySI.mat <- mapply(frailty_function, p0, frailtySI, N/2)
  frailtyIS.mat <- mapply(frailty_function, 0, frailtyIS, N/2)
  # require(ggplot2); ggplot(reshape2::melt(frailtyIS.mat)) + geom_tile(aes(x=Var2, y=max(Var1)-Var1+1, colour=value))
  patient.data <- cbind(patient.data, do.call(rbind, lapply(1:N, function (patient_i) {
    rates.matrix_i <- rates.matrix
    if (patient.data$vacc[patient_i]) {
      rates.matrix_i[2:(ntot+1), 2:(ntypes+1)] <- 
        sweep(rates.matrix[2:(ntot+1), 2:(ntypes+1)], MARGIN=2, thetaSI * frailtySI.mat[patient_i, ], `*`)
      rates.matrix_i[2:(ntypes+1), 1] <- rates.matrix_i[2:(ntypes+1), 1] * frailtyIS.mat[patient_i, ] * thetaIS
    }
    diag(rates.matrix_i)[-1] <- -rowSums(rates.matrix_i[-1, -1])
    observations <- rep(0, nswabs)
    multinom.p <- expm(rates.matrix_i*1e0)[1, ]
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

sim.data <- simulate.data(
  N=1000, ntypes=13, ntot=20, nswabs=5, times=c(365/12, 6*365/12, 7*365/12, 365, 18*365/12),
  lambda=c(0.005,0.0037,0.002,0.0025,0.003,0.0019,0.0015,0.002,0.002,0.0027,0.0018,0.0019,0.002, # vaccine types
           0.003689, 0.002932, 0.003226, 0.002592, 0.003269, 0.002542, 0.003142)*10, # nonvaccine types
  mu=c(0.0333,0.0333,0.0286,0.0286,0.0286,0.0286,0.0222,0.0333,0.0333,0.0333,0.025,0.0143,0.0286, #vaccine types
       0.02914, 0.02736, 0.0272, 0.03303, 0.03161, 0.03369, 0.03009), #nonvaccine types
  p0=c(0.46, 0.3, 0.41, 0.16, 0.175, 0.287, 0.435, 0.178, 0.368, 0.217, 0.389, 0.372, 0.219),
  thetaSI=c(0.632, 0.476, 0.518, 0.617, 0.496, 0.53, 0.54, 0.682, 0.602, 0.639, 0.555, 0.618, 0.507),
  thetaIS=1/c(0.854, 1, 0.723, 0.919, 0.859, 1, 1, 0.891, 1, 0.911, 1, 0.918, 0.84),
  frailtySI=0.2, frailtyIS=0.8, interaction=0.4)

# Compare the number of samples testing positive for vaccine types
Reduce('+', lapply(sim.data[, -1:-2], function (x) c(sum(x[1:(N/2)]<=(ntypes+1)), sum(x[(N/2+1):N]<=(ntypes+1)))))



