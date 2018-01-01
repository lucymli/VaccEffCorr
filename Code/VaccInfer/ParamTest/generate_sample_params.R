params <- list(output_file_name="output.txt",
     n_vtypes=13,
     n_tot=31,
     n_params=102,
     n_ind=200,
     )
param.table <- data.frame(
  name = c(unlist(lapply(c("lambda", "mu"), paste0, 1:params$n_tot)), 
                 paste0("competition", 1:(params$n_vtypes+1)),
                 unlist(lapply(c("ab1type", "ab2type"), paste0, 1:params$n_vtypes))),
  value = c(c(7.183e-05, 0.00083921, 0.00022521, 0.00042459, 0.00019507, 0.00080566, 0.00044165, 2.908e-05, 0.00012329, 1.33e-05, 0.0009141, 3.583e-05, 0.00117243, 0.00061351, 0.00054494, 0.00033651, 0.0004877, 0.00046708, 0.00046319, 0.00048696, 0.00035486, 0.00027702, 0.00018986, 0.00021988, 0.00026736, 0.00041915, 0.00025374, 0.00024932, 0.00023734, 0.00010649, 0.00167282),
            c(0.0242198721780671, 0.00864753078964656, 0.0224452970411401, 0.0143958576110952, 0.0201246742665162, 0.0113002639907398, 0.0147408927700328, 0.0322174470516195, 0.0244905119533151, 0.0242198721780671, 0.00807621461190785, 0.0242198721780671, 0.0169726776048194, 0.0193000146098081, 0.0113524994247617, 0.0184617301502869, 0.0190479978037859, 0.013832595114345, 0.0242198721780671, 0.0181975936725522, 0.0130995077129364, 0.0242198721780671, 0.0242198721780671, 0.0242198721780671, 0.0242198721780671, 0.0242198721780671, 0.0153342255462987, 0.010748053487011, 0.0165794951260709, 0.0242198721780671, 0.0242198721780671),
            rep(1, params$n_vtypes+1),
            rep(0.5, params$n_vtypes),
            rep(-1, params$n_vtypes)),
  transform = c(rep(1, params$n_tot), rep("inv", params$n_tot), rep(1, params$n_vtypes*3+1)),
  min=c(rep(0, params$n_tot*2), rep(0, params$n_vtypes+1), rep(0, params$n_vtypes), rep(-2, params$n_vtypes)),
  max=c(rep(0.1, params$n_tot), rep(1, params$n_tot), rep(1, params$n_vtypes+1), rep(1, params$n_vtypes), rep(0, params$n_vtypes)),
  prior=c(rep("unif", params$n_tot), rep("norm", params$n_tot), rep("unif", params$n_vtypes+1), rep("unif", params$n_vtypes*2))
)

param.table$prior1 <- param.table$min
param.table$prior2 <- param.table$max

param.table[(1:params$n_tot)+params$n_tot, "prior1"] <- 
  c(10.36713, 24.26204, 12.29865, 18.40787, 20.91705, 11.64595, 13.85112, 23.01528, 18.8834, 10.36713, 20.97358, 10.36713, 21.51257, 18.66987, 25.94047, 31.0417, 21.45221, 15.14618, 10.36713, 15.84848, 36.82739, 10.36713, 10.36713, 10.36713, 10.36713, 10.36713, 20.94337, 92.32918, 13.7052, 10.36713, 10.36713)
param.table[(1:params$n_tot)+params$n_tot, "prior2"] <- 
c(41.28841, 115.63995, 44.55276, 69.46443, 49.69025, 88.49351, 67.8385, 31.03908, 40.83214, 41.28841, 123.82038, 41.28841, 58.91822, 51.81343, 88.08633, 54.1661, 52.49896, 72.29301, 41.28841, 54.95232, 76.33875, 41.28841, 41.28841, 41.28841, 41.28841, 41.28841, 65.2136, 93.0401, 60.31547, 41.28841, 41.28841)

output.file.name <- "test_params.txt"
cat(paste(sapply(seq_along(params), function (i) paste(names(params)[[i]], paste(params[[i]], collapse=" "))), collapse="\n"), file=output.file.name)
cat("\nParamTable\n", file=output.file.name, append=TRUE)
invisible(sapply(apply(param.table, 1, paste, collapse=" "), cat, sep="\n", file=output.file.name, append=TRUE))



