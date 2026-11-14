# load("c:/Users/jorda/OneDrive/Documents/Research/NIH/JLC-HMM/Output/Routputs/RoutputsNHANESOakes/JMHMMRandInitReducedFitMix5Seed33len96_oakes_profiled_se.rda")
load("c:/Users/jorda/OneDrive/Documents/Research/NIH/JLC-HMM/Output/Routputs/RoutputsNHANESOakes/JMHMMLoadInFitMix5Seed159len96_oakes_profiled_se.rda")


o <- to_save$oakes

I <- 0.5 * (o$I_obs + t(o$I_obs))
ee <- eigen(I, symmetric = TRUE)

floor_rel <- 1e-6
floor_abs <- floor_rel * max(abs(ee$values))
lambda_repaired <- pmax(ee$values, floor_abs)

V_repaired <- ee$vectors %*%
  diag(1 / lambda_repaired) %*%
  t(ee$vectors)

class_idx <- which(
  o$parameter_map$block == "survival" #& 
    #grepl("^class_", o$parameter_map$param_name)
)

se_class_repaired <- sqrt(diag(V_repaired)[class_idx])
se_class_repaired
beta <- to_save$est_params$beta_vec[-1]

# results <- data.frame(
#   class = paste0("class_", 2:5),
#   beta = beta,
#   se = se_class_repaired,
#   HR = exp(beta),
#   lower = exp(beta - qt(.975,32) * se_class_repaired),
#   upper = exp(beta + qt(.975,32)* se_class_repaired)
# )

to_save$oakes$model_based_survival_se <- se_class_repaired
to_save$oakes$model_based_survival_vcov <- V_repaired[
  class_idx, class_idx, drop = FALSE
]

save(to_save,file = "JMHMMLoadInFitMix5Seed159len96_oakes_profiled_se.rda")
