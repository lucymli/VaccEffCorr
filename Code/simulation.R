#setwd("~/Dropbox (VERG)/Research/Projects/Vaccine/VaccEffCorr/Code")
source("functions.R")
ntypes = 3
ntot = 5
sim.params <- list(N=1000, ntypes=ntypes, ntot=ntot, nswabs=5,
                   times=c(365/12, 6*365/12, 7*365/12, 365, 18*365/12),
     lambda=rlnorm(ntot, log(0.0026696), 2e-1),#c(0.005,0.0037,0.002,0.0025,0.003,0.0019,0.0015,0.002,0.002,0.0027,0.0018,0.0019,0.002, # vaccine types
            #  0.003689, 0.002932, 0.003226, 0.002592, 0.003269, 0.002542, 0.003142), # nonvaccine types
     mu=rlnorm(ntot, log(0.026696), 2e-1),#c(0.0333,0.0333,0.0286,0.0286,0.0286,0.0286,0.0222,0.0333,0.0333,0.0333,0.025,0.0143,0.0286, #vaccine types
        #  0.02914, 0.02736, 0.0272, 0.03303, 0.03161, 0.03369, 0.03009)/10, #nonvaccine types
     thetaSI=c(0.632, 0.476, 0.518, 0.617, 0.496, 0.53, 0.54, 0.682, 0.602, 0.639, 0.555, 0.618, 0.507)[1:ntypes],
     thetaIS=1/c(0.854, 1, 0.723, 0.919, 0.859, 1, 1, 0.891, 1, 0.911, 1, 0.918, 0.84)[1:ntypes],
     p0=rnorm(ntypes, 0.4, 0.1),
     frailtySI=0.2, frailtyIS=0.2, interaction=0.4)
sim.params.vec <- unlist(sim.params[-1:-5])
sim.params.sd <- c(lambda=0.005, mu=0.005, thetaSI=0.3, thetaIS=0.1, p0=0.2,
                   frailtySI=0.2, frailtyIS=0.2, interaction=0.2)
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
save.image(file="sim.data.RData")
cor(sim.data$abdata)
plot_timeseries(sim.data, sim.params$ntypes)
plot(data.frame(sim.data$abdata))

# t test to determine if there is a signficant difference in the mean IgG concentrations
# between fully protected individuals and partially protected individuals
test_ablevels_protection(sim.data$abdata, sim.params$p0, sim.params$N/2, 1:sim.params$ntypes, return.p=TRUE) < 0.05
# test_ablevels_protection(sim.data$abdata, sim.params$p0, sim.params$N/2, 9, return.p = TRUE)

# Plot the distributions of IgG concentrations
plot_abdata(sim.data$abdata, sim.params$p0, sim.params$N/2)

ab.levels.df <- do.call(rbind, lapply(1:ntypes, function(typei) {
  df.ts <- apply(sim.data$vdata[, -1:-2], 2, function (x) {
    df1 <- data.frame(value=sim.data$abdata[x==typei, typei])
    df2 <- data.frame(value=sim.data$abdata[x!=typei, typei])
    if (nrow(df1)==0) return (cbind(df2, group="not infected"))
    else if (nrow(df2)==0) return (cbind(df1, group="infected"))
    else return (rbind(cbind(df1, group="infected"), cbind(df2, group="not infected")))
  })
  do.call(rbind, mapply(cbind, df.ts, time=1:length(df.ts), type=paste("serotype", typei), SIMPLIFY = FALSE))
}))

ggplot(ab.levels.df, aes(x=factor(time), y=value, fill=group)) + theme_bw() +
  geom_point(pch=21, position=position_jitterdodge()) +
  geom_boxplot(aes(x=factor(time), y=value, fill=group), alpha=.5) +
  facet_grid(type~.) +
  ylab("Antibody levels") + xlab("Swab #")
ggplot(ab.levels.df, aes(x=factor(time), y=value, fill=group)) + theme_bw() +
  geom_point(pch=21, position=position_jitterdodge()) +
  geom_boxplot(aes(x=factor(time), y=value, fill=group), alpha=.5) +
  ylab("Antibody levels") + xlab("Swab #")

