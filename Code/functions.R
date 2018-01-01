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

simulate_data <- function (
  Nv, Nnv, # numbers of participants in the vaccine/nonvaccine arms
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
  ab.mean,  # mean antibody level (log10)
  ab.var # variance in antibody level
) {
  ### For testing within the function
  # load("simulate.data.test.RData")
  ###
  require(parallel)
  require(Matrix)
  N <- Nv + Nnv
  patient.data <- data.frame(ID=1:(N), vacc=c(rep(TRUE, Nv), rep(FALSE, Nnv)))
  # The rates of clearance or infection by another serotype
  interactions <- c(rep(interaction[1], ntypes), rep(interaction[2], ntot-ntypes))
  rates.matrix <- t(replicate(ntot, lambda*interactions))
  rates.matrix <- rbind(c(-sum(lambda), lambda),
                        cbind(mu, rates.matrix))
  for (i in 1:nrow(rates.matrix)) rates.matrix[i,i] <- -sum(rates.matrix[i,-i])
  # # Individual heterogeneities in rates of acquisition(frailtySI)/clearance(frailtyIS) after being vaccinated
  # frailtySI.mat <- mapply(frailty_function, p0, frailtySI, Nv)
  # frailtyIS.mat <- mapply(frailty_function, 0, frailtyIS, Nv)
  # Simulate antibody data:
  antibody.measurements <- data.frame(t(sapply(1:Nv, function (i) rnorm(ntypes, rnorm(1, ab.mean, ab.mean), sqrt(ab.var)))))
  colnames(antibody.measurements) <- paste0("IgG.", 1:ntypes)
  ###
  ### Check that the frailty matrices have been simulated correctly
  # require(ggplot2); ggplot(reshape2::melt(frailtyIS.mat)) + geom_tile(aes(x=Var2, y=max(Var1)-Var1+1, colour=value))
  ###
  # Simulate data for each patient
  patient.data <- cbind(patient.data, do.call(rbind, lapply(1:N, function (patient_i) {
    # Construct vector for observations
    observations <- rep(0, nswabs)
    multinom.p <- expm(rates.matrix*60)[1, ]
    # Simulate first observation
    observations[1] <- which(rmultinom(N, 1, multinom.p)==1, arr.ind = TRUE)[patient_i, 1]
    rates.matrix_i <- rates.matrix
    # Modify rates for vaccinated individuals
    if (patient.data$vacc[patient_i]) {
      rates.matrix_i[2:(ntot+1), 2:(ntypes+1)] <- 
        sweep(rates.matrix[2:(ntot+1), 2:(ntypes+1)], MARGIN=2, unlist(p0*exp(thetaSI*antibody.measurements[patient_i, ])), `*`)
      # rates.matrix_i[2:(ntypes+1), 1] <- rates.matrix_i[2:(ntypes+1), 1] * thetaIS #* frailtyIS.mat[patient_i, ]
    }
    rates.matrix_i[rates.matrix_i>1] <- 1
    rates.matrix_i[rates.matrix_i<0] <- 0
    for (i in 1:nrow(rates.matrix_i)) rates.matrix[i,i] <- -sum(rates.matrix_i[i,-i])
    P <- expm(rates.matrix_i*times[1])
    for (time.step in 2:nswabs) {
      observations[time.step] <- sample(ntot+1, 1, prob=P[observations[time.step-1], ])
      P <- expm(rates.matrix_i*(times[time.step]-times[time.step-1]))
    }
    return (observations-1)
  })))
  # antibody.measurements <- t(sapply(1:(N/2), function (i) {
  #   out <- (max(frailtySI.mat)-(frailtySI.mat[i, 1:ntypes]+rnorm(ntypes, 5, ab.noise)))/2+1
  #   return (out)
  # }))
  # patient.data <- data.frame(patient.data, antibody.measurements)
  # if (ab.noise!=0) antibody.measurements <- antibody.measurements * rnorm(length(antibody.measurements), sd=ab.noise)
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

plot_timeseries <- function (dataset, ntypes) {
  V <- apply(dataset$vdata[, -1:-2], 2, function (x) c(sum(x %in% 1:ntypes), sum(x > ntypes), sum(x==0)))
  NV <- apply(dataset$nvdata[, -1:-2], 2, function (x) c(sum(x %in% 1:ntypes), sum(x > ntypes), sum(x==0)))
  V.df <- data.frame(time=rep(1:ncol(V), nrow(V)), count=c(t(V)), serotype=rep(c("Vaccine type", "Non-vaccine type", "No S.pneumo"), each=ncol(V)), set="Vaccine")
  NV.df <- data.frame(time=rep(1:ncol(NV), nrow(NV)), count=c(t(NV)), serotype=rep(c("Vaccine type", "Non-vaccine type", "No S.pneumo"), each=ncol(NV)), set="No Vaccine")
  combined.df <- rbind(V.df, NV.df)
  require(ggplot2)
  ggplot(combined.df) + theme_bw() +
    geom_bar(aes(x=time, y=count, fill=set), stat="identity", position="dodge", alpha=.85) +
    facet_grid(serotype~.)
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



mcmc_infer <- function (params, params.sd, dataset, mcmc.options, 
                        file.out=paste0("output", gsub("[^a-zA-Z0-9]", "", Sys.time()), ".txt"), 
                        mcmc.program.dir="VaccInfer/VaccInfer",
                        program.name="mcmc_infer",
                        input.file="input.txt",
                        run.program=FALSE) {
  selected.params <- params[c("lambda", "mu", "thetaSI", "p0", "interaction")]
  tot_params <- length(unlist(selected.params))#-params$ntypes
  param.vec <- paste(gsub(" ", "", format(unlist(selected.params), scientific=FALSE)), collapse="\n")
  pre.vec <- paste(c(vaccN=params$Nv, unvaccN=params$Nnv, n_vt=params$ntypes, n_nvt=params$ntot-params$ntypes,
    total_params=tot_params, n_swabs=params$nswabs, unlist(mcmc.options), file.out 
    #n_param_blocks=length(selected.params), sapply(selected.params, length)
    ), collapse="\n")
  param.sd.vec <- paste(unlist(params.sd), collapse="\n")
  data.vec <- paste(unlist(lapply(c(lapply(dataset[1:2], `[`, , -1:-2), dataset[3]), c)), collapse="\n")
  swab.times.vec <- paste(params$times, collapse="\n")
  input.vec <- paste(pre.vec, param.vec, param.sd.vec, data.vec, swab.times.vec, sep="\n")
  cat(input.vec, file=input.file)
  command <- paste0(mcmc.program.dir, "/./", program.name, " ", input.file)
  if (run.program) {
    output <- system(command, intern=TRUE, wait=TRUE)
  }
  return (command)
}


get.ab.lambda.reduction <- function (ab.values, mcmc.output, selected, selected.names=NULL, as.percent=FALSE) {
  p0 <- mcmc.output[, paste0("p0", selected-1)]
  thetaSI <- mcmc.output[, paste0("thetaSI", selected-1)]
  reductions <- 1-(p0-thetaSI*log(ab.values[2]))/(p0-thetaSI*log(ab.values[1]))
  reductions[reductions>1] <- 1
  if (length(selected) > 1) {
    MAP <- apply(reductions, 2, get.MAP)
    bounds <- apply(reductions, 2, EpiGenR::hpd, show.median=FALSE)
    results <- cbind(MAP, t(bounds))
  } else {
    results <- c(get.MAP(reductions), EpiGenR::hpd(reductions, show.median=FALSE))
  }
  results[results>1] <- 1
  if (as.percent) {
    if (length(selected) > 1) {
      results <- apply(results, 2, scales::percent)
      if (!is.null(selected.names)) rownames(results) <- selected.names
    } else {
      results <- scales::percent(results)
    }
  }
  return (results)
}

get_first_serotype <- function (x) {
  output <- strsplit(unlist(x), ",", fixed=TRUE) %>% sapply(`[`, 1)
  output[is.na(output)] <- ""
  return (output)
}

extract_serotype <- function (swabs, ignore=NULL) {
  require(magrittr)
  subswabs <- swabs
  if (!is.null(ignore)) subswabs <- swabs[ignore*(-1)]
  serotypes <- get_first_serotype(subswabs)
  sapply(serotypes, function (z) {
    NEW <- TRUE
    if (!is.null(ignore)) {
      if (z %in% (strsplit(swabs[ignore], ",", fixed=TRUE)%>%unlist)) NEW <- FALSE
    }
    if (NEW) return (z)
    else return ("")
  })
}


get_prevalence <- function (swabs, ignore=1:3) {
  require(magrittr)
  apply(swabs, 1, function (x) {
    old.serotypes <- extract_serotype(swabs=x[ignore])
    new.serotypes <- extract_serotype(swabs=x, ignore=ignore)
    unique(sapply(new.serotypes, function (z) ifelse(z %in% old.serotypes, "", z)))
  }) %>% 
    unlist %>%
    table() %>%
    sort(decreasing=TRUE)
}

get_serotype_id <- function (value, key.table) {
  if (length(value) > 1) return(sapply(value, get_serotype_id, key.table))
  else {
    if (nchar(value)==0) return (0)
    if (any(key.table[, 2]==value)) index <- which(key.table[, 2]==value)
    else index <- nrow(key.table)
    return (key.table[index, 1])
  }
}

get_serotype_name <- function (id, key.table) {
  if (length(id) > 1) return(sapply(id, get_serotype_name, key.table))
  else {
    if (id==0) return ("")
    index <- which(key.table$key==id)
    return(key.table[index, 2])
  }
}
