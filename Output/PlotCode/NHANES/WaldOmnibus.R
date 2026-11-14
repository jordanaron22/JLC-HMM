jkn <- readRDS("Output/parse_nhanes_jkn.rds")

beta <- jkn$reference_estimates$joint
V    <- jkn$vcov$joint  


W <- drop(t(beta) %*% solve(V, beta))
Fstat <- W / 4
p_F <- pf(Fstat, df1 = 4, df2 = 32, lower.tail = FALSE)
