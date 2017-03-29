#setwd("~/Dropbox (hsph.harvard.edu)/Research/Projects/Vaccine/VaccEffCorr/Code")
source("functions.R")
ntypes = 13
ntot = 31
day.units = 1
nswabs = 6
sim.params <- list(Nv=200, Nnv=200, ntypes=ntypes, ntot=ntot, nswabs=nswabs,
                   times=(c(2, 7, 12, 13, 18, 24)[1:nswabs])*30/day.units,
     lambda=c(rlnorm(ntypes, log(0.0003696), 2e-1), rlnorm(ntot-ntypes, log(0.0003696), 2e-1))*day.units,
     mu=rlnorm(ntot, log(0.016696), 2e-1)*day.units,
     thetaSI=rnorm(ntypes, 0.2505, 0.06),
     thetaIS=rep(1, ntypes),#1/rnorm(ntypes, 1, 0.06)[1:ntypes],
     p0=rep(0, ntypes),#abs(rnorm(ntypes, 0.5, 0.05)),
     frailtySI=0.2, frailtyIS=0.2, interaction=0.99)
sim.params.vec <- unlist(sim.params[-1:-5])
sim.params.sd <- list(lambda=rep(mean(sim.params$lambda)*0.01, ntot), 
                      mu=rep(mean(sim.params$mu)*0.01, ntot), 
                      #thetaSI=0.3, thetaIS=0.1, p0=0.2, frailtySI=0.2, frailtyIS=0.2, 
                      interaction=0.2)
sim.params$ab.noise <- runif(sim.params$ntypes, 0.0, 0.05)

theta.SI.agr <- with(sim.params, get_agr_SI(lambda[1:ntypes], thetaSI))
theta.IS.agr <- with(sim.params, get_agr_IS(lambda, mu, thetaSI, thetaIS, ntypes))

# variability <- sapply(20170109:20170139, function (x) {
#   set.seed(x)
#   sim.data <- do.call(simulate_data, sim.params)
#   sim.data.summ <- summarize_data(sim.data, sim.params$ntypes)
#   colSums(sim.data.summ$mean.num.samples)
# })
# apply(variability, 1, sd)
set.seed(2100)
sim.data <- do.call(simulate_data, sim.params)
mcmc_options <- list(niter="100000", sample_every=50, adapt_optimal=0.23, adapt_every=10,
                  adapt_until="5000")

mcmc_infer (sim.params, sim.params.sd, sim.data, mcmc_options, 
            mcmc.program.dir="VaccInfer/VaccInfer", program.name="mcmc_infer",
            run.program = FALSE)

save.image(file="sim.data.RData")



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



