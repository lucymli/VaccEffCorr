#setwd("~/Dropbox (hsph.harvard.edu)/Research/Projects/Vaccine/VaccEffCorr/Code")
source("functions.R")
library(magrittr)
library(parallel)
ntypes = 13
ntot = 31
day.units = 1
nswabs = 3
load("parse_data_new.RData")
sim.params <- list(Nv=700, Nnv=0, ntypes=ntypes, ntot=ntot, nswabs=nswabs,
                   times=(c(6, 7, 12)[1:nswabs])*30/day.units,
     lambda=c(7.183e-05, 0.00083921, 0.00022521, 0.00042459, 0.00019507, 0.00080566, 0.00044165, 2.908e-05, 0.00012329, 1.33e-05, 0.0009141, 3.583e-05, 0.00117243, 0.00061351, 0.00054494, 0.00033651, 0.0004877, 0.00046708, 0.00046319, 0.00048696, 0.00035486, 0.00027702, 0.00018986, 0.00021988, 0.00026736, 0.00041915, 0.00025374, 0.00024932, 0.00023734, 0.00010649, 0.00167282), #c(rlnorm(ntypes, log(0.0003696), 2e-1), rlnorm(ntot-ntypes, log(0.0003696), 2e-1))*day.units,
     mu=serotypes.key$mu, #rlnorm(ntot, log(0.016696), 2e-1)*day.units,
     thetaSI=-rnorm(ntypes, runif(1, 0, 1), 0.06),
     thetaIS=rep(1, ntypes),#1/rnorm(ntypes, 1, 0.06)[1:ntypes],
     p0=sapply(rnorm(ntypes, runif(1, 0, 1), 0.3), function (x) {
       if (x < 0) return (0)
       if (x > 1) return (1)
       return (x)
     }),
     frailtySI=0.2, frailtyIS=0.2, interaction=c(1, 1),
     ab.mean=0.3, ab.var=0.45)

reps <- 24
simulated.data <- mclapply(1:reps, function (i) do.call(simulate_data, sim.params), mc.cores=4)
simulated.prev <- lapply(simulated.data, function (x) table(unlist(x$vdata[-1:-2])) %>% prop.table) %>% do.call(what=rbind) %>% apply(2, summary)
data.prev <- table(unlist(DATA$vdata)) %>% prop.table
#sim.params$lambda <- sim.params$lambda * data.prev[-1]/simulated.prev[4, 2:32]


sim.params.vec <- unlist(sim.params[-1:-5])
sim.params.sd <- list(lambda=rep(mean(sim.params$lambda)*0.01, ntot), 
                      mu=rep(mean(sim.params$mu)*0.01, ntot), 
                      #thetaSI=0.3, thetaIS=0.1, p0=0.2, frailtySI=0.2, frailtyIS=0.2, 
                      interaction=0.2)

if (FALSE) {
  lambda.mu.combos <- expand.grid(lambda=c(1e-4, 5e-4, 7.5e-4, 1e-3, 2.5e-3, 5e-3),
                                  mu=c(1e-3, 5e-3, 1e-2, 5e-2, 1e-1))
  lambda.mu.combos[apply(lambda.mu.combos, 1, function (x) x[1]/x[2])<=1, ]
  lambda.mu.sims <- mclapply(1:nrow(lambda.mu.combos), function (i) {
    x <- unlist(lambda.mu.combos[i, ])
    sim.params$lambda <- abs(c(rnorm(ntypes, x[1], x[1]/2), rnorm(ntot-ntypes, x[1], x[1]/2))*day.units)
    sim.params$mu <- abs(rnorm(ntot, x[2], x[2]/2)*day.units)
    simulated.curve <- rowSums(sapply(1:10, function (j) {
      sim.data <- list(param=sim.params, data=do.call(simulate_data, sim.params))
      simulated.curve <- colSums(sim.data$data$vdata[, -1:-2]>0)/nrow(sim.data$data$vdata)
    }))/10
    data.frame(curve=paste0("sim_lamb_", x[1], "_mu_", x[2]), x=1:6, y=simulated.curve)
  }, mc.cores=4)
  data.curve <- colSums(combined.data$PCV13$vdata[, -1:-2]>0)/nrow(combined.data$PCV13$vdata)
  library(ggplot2)
  lambda.mu.df <- rbind(do.call(rbind, lambda.mu.sims), data.frame(curve="data", x=1:6, y=data.curve))
  lambda.mu.df$type <- sapply(strsplit(as.character(lambda.mu.df$curve), "_"), `[`, 1)
  ggplot(lambda.mu.df) + theme_bw() + 
    geom_line(aes(x=x, y=y, colour=curve, size=type)) +
    scale_size_manual(values=c(3, 1))# +
  #scale_color_manual(values=c(rep("violetred3", nrow(lambda.mu.combos)), "mediumblue"))
  
  # The lambda and mu combos that give the most similar prevalence over time:
  best.fit.lambda.mu <- lambda.mu.combos[which.min(sapply(lambda.mu.sims, function (sim) sum((sim$y-data.curve)^2))), ]
}

best.fit.lambda.mu <- c(lambda=0.00075, mu=0.05)


num.datasets <- 24
ab.noise.levels <- c(0, 0.01, 0.05, 0.1, 0.5, 1)
set.seed(2100)
sim.data.all.serotypes <- lapply(1:num.datasets, function (i) {
  sim.params$lambda <- abs(c(rnorm(ntypes, best.fit.lambda.mu[1], best.fit.lambda.mu[1]/2), 
                             rnorm(ntot-ntypes, best.fit.lambda.mu[1], best.fit.lambda.mu[1]/2))*day.units)
  ratio <- best.fit.lambda.mu[1]/best.fit.lambda.mu[2]
  sim.params$mu <- sim.params$lambda*rnorm(length(sim.params$lambda), ratio, ratio/10)#abs(rnorm(ntot, 1e1, 1e1/2)*day.units)
  sim.params$thetaSI <- -rnorm(ntypes, runif(1, 0, 5), 0.06)
  sim.params$p0 <- sapply(rnorm(ntypes, runif(1, 0, 1), 0.3), function (x) {
    if (x < 0) return (0)
    if (x > 1) return (1)
    return (x)
  })
  sim.data <- list(param=sim.params, data=do.call(simulate_data, sim.params)) # takes around 10.5 seconds per simulation
  data.reps <- lapply(ab.noise.levels, function (ab.noise.level) {
    sim.data$data$abdata <- sim.data$data$abdata + rnorm(length(sim.data$data$abdata), sd=ab.noise.level)
    return (sim.data)
  })
  return (data.reps)
})
# set.seed(2000)
# sim.data.ind.var <- mclapply(ab.noise.levels, function (ab.noise.level) {
#   sim.params$lambda <- c(rlnorm(ntypes, log(0.0003696), 2e-1), rlnorm(ntot-ntypes, log(0.0003696), 2e-1))*day.units
#   sim.params$mu <- rlnorm(ntot, log(0.016696), 2e-1)*day.units
#   sim.params$thetaSI <- rep(rnorm(ntypes, 1.5, 0.3), ntypes)
#   sim.params$p0 <- rep(rnorm(1, 0.2505, 0.06), ntypes)
#   sim.data <- list(param=sim.params, data=do.call(simulate_data, sim.params))
#   data.reps <- lapply(ab.noise.levels, function (ab.noise.level) {
#     sim.data$data$abdata <- sim.data$data$abdata + rnorm(length(sim.data$data$abdata), sd=ab.noise.level)
#     return (sim.data)
#   })
#   return (data.reps)
# }, mc.cores=4)


mcmc_options <- list(niter="500000", sample_every=100, adapt_optimal=0.23, adapt_every=10,
                  adapt_until="5000", use_mean_ab=0)
if ("sim.data.RData" %in% list.files()) {
  system("stat '/Users/Lucy/Dropbox (hsph.harvard.edu)/Research/Projects/Vaccine/VaccEffCorr/Code/sim.data.RData'", intern=TRUE) -> tmp
  newname <- paste0("sim.data", strsplit(tmp, '"')[[1]][4], ".RData")
  system(paste0("mv sim.data.RData ", newname))
}
save.image(file="sim.data.RData")
# 
# lapply(seq_along(sim.data.all.serotypes), function (i) {
#   lapply(seq_along(sim.data.all.serotypes[[i]]), function (j) {
#     mcmc_infer (sim.data.all.serotypes[[i]][[j]]$param, sim.params.sd, 
#                 sim.data.all.serotypes[[i]][[j]]$data, mcmc_options, 
#                 mcmc.program.dir="VaccInfer/VaccInfer", program.name="mcmc_infer",
#                 run.program = FALSE, 
#                 input.file=paste0("sim_data/all_sero_rep_", i, "_abnoise_", ab.noise.levels[j], ".txt"),
#                 file.out=paste0("sim_data/output/all_sero_rep_", i, "_abnoise_", ab.noise.levels[j], ".txt"))
#     mcmc_options$use_mean_ab=1
#     mcmc_infer (sim.data.all.serotypes[[i]][[j]]$param, sim.params.sd, 
#                 sim.data.all.serotypes[[i]][[j]]$data, mcmc_options, 
#                 mcmc.program.dir="VaccInfer/VaccInfer", program.name="mcmc_infer",
#                 run.program = FALSE, 
#                 input.file=paste0("sim_data/meanab_rep_", i, "_abnoise_", ab.noise.levels[j], ".txt"),
#                 file.out=paste0("sim_data/output/meanab_rep_", j, "_abnoise_", ab.noise.levels[j], ".txt"))
#   })
# })
# lapply(seq_along(sim.data.ind.var), function (i) {
#   lapply(seq_along(sim.data.ind.var[[i]]), function (j) {
#     mcmc_infer (sim.data.ind.var[[i]][[j]]$param, sim.params.sd, 
#                 sim.data.ind.var[[i]][[j]]$data, mcmc_options, 
#                 mcmc.program.dir="VaccInfer/VaccInfer", program.name="mcmc_infer",
#                 run.program = FALSE, 
#                 input.file=paste0("sim_data/ind_var_rep_", j, "_abnoise_", ab.noise.levels[i], ".txt"),
#                 file.out=paste0("output/ind_var_rep_", j, "_abnoise_", ab.noise.levels[i], ".txt"))
#   })
# })

# mcmc_infer (sim.params, sim.params.sd, sim.data, mcmc_options,
#             mcmc.program.dir="VaccInfer/VaccInfer", program.name="mcmc_infer",
#             run.program = FALSE)




if (FALSE) {
  library(ggplot2)
  
  cor(sim.data$abdata)
  # plot_timeseries(sim.data, sim.params$ntypes)
  # plot(data.frame(sim.data$abdata))
  
  # t test to determine if there is a signficant difference in the mean IgG concentrations
  # between fully protected individuals and partially protected individuals
  test_ablevels_protection(sim.data$abdata, sim.params$p0, sim.params$N/2, 1:sim.params$ntypes, return.p=TRUE) < 0.05
  # test_ablevels_protection(sim.data$abdata, sim.params$p0, sim.params$N/2, 9, return.p = TRUE)
  
  # Plot the distributions of IgG concentrations
  # plot_abdata(sim.data$abdata, sim.params$p0, sim.params$N/2)
  
  ab.levels.df <- data.frame(ab=sim.data$abdata[, 1], apply(sim.data$vdata[, -1:-2], 2, function (x) c("not infected", "infected by vt", "infected by nvt")[findInterval(x, c(1, ntypes+1))+1]))
  ab.levels.df <- reshape2::melt(ab.levels.df, id.vars="ab", variable.name="time", value.name="group")
  
  # ggplot(ab.levels.df, aes(x=factor(time), y=ab, fill=group)) + theme_bw() +
  #   geom_point(pch=21, position=position_jitterdodge()) +
  #   geom_boxplot(aes(x=factor(time), y=ab, fill=group), alpha=.5) +
  #   facet_grid(type~.) +
  #   ylab("Antibody levels") + xlab("Swab #")
  
  ### Real Data
  ab.data <- readRDS("/Users/Lucy/Dropbox (hsph.harvard.edu)/Research/Projects/Vaccine/VaccEffCorr/Documents/Pfizer/final_version/data.rds")
  rbind(colSums(subset(ab.data, rtrtn=="PCV13")[, paste0("nvt.v", 1:8)])/nrow(subset(ab.data, rtrtn=="PCV13")),
        colSums(subset(ab.data, rtrtn=="PCV13")[, paste0("vt.v", 1:8)])/nrow(subset(ab.data, rtrtn=="PCV13")))
  ### Simulated Data
  do.call(cbind, lapply(lapply(split(ab.levels.df$group, ab.levels.df$time), table), prop.table))
  
  
  ggplot(data.frame(ab=ab.data$mean.igg.v4, 
                    infected=factor(nchar(as.character(ab.data$sertyp.v4))!=0)), 
         aes(x=infected, y=ab, fill=infected)) + theme_bw() +
    geom_point(pch=21, position=position_jitterdodge()) + 
    geom_boxplot(alpha=.4) + scale_y_log10()
  ggplot(ab.levels.df, aes(x=factor(time), y=ab, fill=group)) + theme_bw() +
    geom_point(pch=21, position=position_jitterdodge()) +
    geom_boxplot(aes(x=factor(time), y=ab, fill=group), alpha=.5) +
    ylab("Antibody levels") + xlab("Swab #") +
    scale_y_log10()
}



