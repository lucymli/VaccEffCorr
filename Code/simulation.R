source("functions.R")

sim.params <- list(N=1000, ntypes=13, ntot=20, nswabs=5, times=c(365/12, 6*365/12, 7*365/12, 365, 18*365/12),
     lambda=rlnorm(20, log(0.0026696), 2e-1),#c(0.005,0.0037,0.002,0.0025,0.003,0.0019,0.0015,0.002,0.002,0.0027,0.0018,0.0019,0.002, # vaccine types
            #  0.003689, 0.002932, 0.003226, 0.002592, 0.003269, 0.002542, 0.003142), # nonvaccine types
     mu=rlnorm(20, log(0.0026696), 2e-1),#c(0.0333,0.0333,0.0286,0.0286,0.0286,0.0286,0.0222,0.0333,0.0333,0.0333,0.025,0.0143,0.0286, #vaccine types
        #  0.02914, 0.02736, 0.0272, 0.03303, 0.03161, 0.03369, 0.03009)/10, #nonvaccine types
     p0=c(0.16, 0.3, 0.141, 0.16, 0.175, 0.287, 0.135, 0.178, 0.308, 0.217, 0.319, 0.322, 0.219)*2.5,
     thetaSI=c(0.632, 0.476, 0.518, 0.617, 0.496, 0.53, 0.54, 0.682, 0.602, 0.639, 0.555, 0.618, 0.507),
     thetaIS=1/c(0.854, 1, 0.723, 0.919, 0.859, 1, 1, 0.891, 1, 0.911, 1, 0.918, 0.84),
     frailtySI=0.2, frailtyIS=0.2, interaction=0.4)
sim.params$ab.noise <- runif(sim.params$ntypes, 0.2, 0.6)

theta.SI.agr <- with(sim.params, get_agr_SI(lambda[1:ntypes], thetaSI))
theta.IS.agr <- with(sim.params, get_agr_IS(lambda, mu, thetaSI, thetaIS, ntypes))

# variability <- sapply(20170109:20170139, function (x) {
#   set.seed(x)
#   sim.data <- do.call(simulate_data, sim.params)
#   sim.data.summ <- summarize_data(sim.data, sim.params$ntypes)
#   colSums(sim.data.summ$mean.num.samples)
# })
# apply(variability, 1, sd)
set.seed(1000)
sim.data <- do.call(simulate_data, sim.params)
#save(sim.data, file="sim.data.RData")
cor(sim.data$ab.data)
