frailty_function <- function (p0, delta, N) {
  if ((p0 < 0) | (p0 > 1)) stop("p0 must be between 0 and 1", call. = TRUE)
  num.fully.protected <- floor(p0*N)
  susceptibility <- rep(0, num.fully.protected)
  if (num.fully.protected<N) {
    susceptibility <- c(susceptibility, rgamma(N-num.fully.protected, 1/delta, 1/delta))
  }
  return(susceptibility)
}

is_fully_protected <- function (dataset, type.ids) {
  apply(dataset[, -1:-2], 1, function (x) !any(x %in% (type.ids)))
}

reorder_by_protection <- function (dataset, ntypes) {
  protected <- is_fully_protected(dataset, 1:ntypes)
  dataset[protected, ] <- dataset[protected, ][order(apply(dataset[protected, -1:-2], 1, function (x) sum(x > 0))), ]
  dataset[!protected, ] <- dataset[!protected, ][order(apply(dataset[!protected, -1:-2], 1, function (x) sum(x > 0))), ]
  return (dataset)
}
simulate_data <- function (
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
  interaction, # rate at which the current
  ab.noise # noisiness in antibody levels 
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
  patient.data <- cbind(patient.data, do.call(rbind, mclapply(1:N, function (patient_i) {
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
  }, mc.cores=detectCores())))
  random.frailties <- rgamma(N/2, 1/frailtySI, 1/frailtySI)
  antibody.measurements <- t(sapply(1:length(random.frailties), function (i) {
    out <- random.frailties[i]+rnorm(ntypes, -1.5, ab.noise)
    protected.serotypes <- which(i<=floor(p0*N/2))
    if (length(protected.serotypes) > 0) {
      out[protected.serotypes] <- out[protected.serotypes] + log(abs(rnorm(length(protected.serotypes), 1.5)))
    }
    return (out)
  }))
  # antibody.measurements <- rbind(antibody.measurements, matrix(NA, nrow=N/2, ncol=ntypes))
  colnames(antibody.measurements) <- paste0("IgG.", 1:ntypes)
  # patient.data <- data.frame(patient.data, antibody.measurements)
  return(list(vdata=patient.data[patient.data$vacc, ], 
              nvdata=patient.data[!patient.data$vacc, ],
              abdata=antibody.measurements))
}

get_mean_positive <- function (dataset, ntypes) {
  out.positive <- sapply(by(dataset[, -1:-2]>0, dataset$vacc, function (x) rowSums(x)/ncol(x)), mean)
  out.positive.nv <- sapply(by(dataset[, -1:-2]>ntypes, dataset$vacc, function (x) rowSums(x)/ncol(x)), mean)
  out.positive.v <- out.positive-out.positive.nv
  out <- rbind(vtypes=out.positive.v, nvtypes=out.positive.nv)
  colnames(out) <- sapply(colnames(out), ifelse, "vacc", "not.vacc")
  return (out)
}

get_ts_positive <- function (dataset, ntypes) {
  # mat should be TRUE/FALSE
  ts.positive <- do.call(cbind, by(dataset[, -1:-2]>0, dataset$vacc, function (x) colSums(x)/nrow(x)))
  ts.positive.nv <- do.call(cbind, by(dataset[, -1:-2]>ntypes, dataset$vacc, function (x) colSums(x)/nrow(x)))
  ts.positive.v <- ts.positive-ts.positive.nv
  out <- rbind(data.frame(time=1:nrow(ts.positive.v), serotype="v", ts.positive.v),
        data.frame(time=1:nrow(ts.positive.nv), serotype="nv", ts.positive.nv))
  out <- reshape2::melt(out, id.vars=c("time", "serotype"), variable.name="vacc", value.name="mean.positive")
  out$vacc <- as.logical(gsub(".", "", out$vacc, fixed=TRUE))
  return (out)
}

summarize_data <- function (dataset, ntypes, include.graphics=TRUE) {
  mean.positive <- get_mean_positive(dataset, ntypes)
  ts.positive <- get_ts_positive(dataset, ntypes)
  if(include.graphics) {
    require(ggplot2)
    P <- ggplot(ts.positive) + theme_bw() +
      geom_line(aes(x=time, y=mean.positive, colour=serotype)) + 
      ylab(expression(paste("Proportion carrying ", italic("S. pneumoniae")))) +
      xlab("Swab #") +
      scale_colour_manual("Serotype", values=c("#d4006f", "#4a3d6b"), labels=c("Vaccine", "Non-vaccine")) +
      facet_grid(.~vacc, labeller=as_labeller(c(`TRUE`="Vaccinated", `FALSE`="Not Vaccinated")))
  }
  return (list(mean.num.samples=mean.positive, timeseries=ts.positive, plot=P))
}

get_agr_SI <- function (lambda, thetaSI) {
  theta.SI.agr <- sum(lambda*thetaSI)/sum(lambda)
  return (theta.SI.agr)
}

get_agr_IS <- function (lambda, mu, thetaSI, thetaIS, ntypes) {
  ntot <- length(lambda)
  B <- thetaSI*lambda[1:ntypes]/(thetaSI*lambda[1:ntypes]+mu[1:ntypes])
  A <- thetaIS*mu[1:ntypes]*B
  D <- lambda[(ntypes+1):ntot]/(lambda[(ntypes+1):ntot]+mu[(ntypes+1):ntot])
  E <- mu[(ntypes+1):ntot]*D
  theta.IS.agr <- sum(A)/sum(B) * sum(D)/sum(E)
  return (theta.IS.agr)
}


test_ablevels_protection <- function (abdata, p0, N, types.ids, return.p=FALSE) {
  p.values <- sapply(types.ids, function (i) {
    protected.selection <- 1:floor(p0[i]*N)
    out <- c()
    out["t.test"] <- t.test(abdata[protected.selection, i], abdata[-protected.selection, i])$p.value
    out["ks.test"] <- ks.test(abdata[protected.selection, i], abdata[-protected.selection, i])$p.value
    return (out)
  })
  if (return.p) return(p.values)
  else return (p.values < 0.05)
}

plot_abdata <- function (abdata, p0, N) {
  abdata.long <- reshape2::melt(abdata)
  abdata.long$protected <- 
    unlist(lapply(p0, function (x) c(rep("full", floor(x*N)), rep("partial", N-floor(x*N)))))
  require(ggplot2)
  ggplot(abdata.long) + theme_bw() +
    geom_density(aes(x=value, fill=protected), alpha=.4) +
    facet_grid(Var2~.) +
   xlab("IgG Concentration")
}



