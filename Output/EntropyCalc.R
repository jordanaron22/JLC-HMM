re_prob <- to_save[[2]][[10]]


#for application
#just 2 percent of inviduals have a max prob < .99 
sum(apply(re_prob,1,max) < .99)/nrow(re_prob)

mean_entropy <- ifelse(re_prob > 0, re_prob * log(re_prob), 0) |>
  rowSums() |> 
  mean() * -1

mean_entropy

load("c:/Users/jorda/OneDrive/Documents/Research/NIH/JLC-HMM/Data/JMHMMNoSurvFitMix5Seed.rda")
re_prob5_ns <- to_save[[2]][[10]]
assignment5_ns <- apply(re_prob5_ns,1,which.max)

load("c:/Users/jorda/OneDrive/Documents/Research/NIH/JLC-HMM/Data/JMHMMFitMix5Seed.rda")
re_prob5 <- to_save[[2]][[10]]
assignment5 <- apply(re_prob5,1,which.max)


load("c:/Users/jorda/OneDrive/Documents/Research/NIH/JLC-HMM/Data/JMHMMFitMix4Seed.rda")
re_prob4 <- to_save[[2]][[10]]
assignment4 <- apply(re_prob4,1,which.max)


load("c:/Users/jorda/OneDrive/Documents/Research/NIH/JLC-HMM/Data/JMHMMFitMix6Seed.rda")
re_prob6 <- to_save[[2]][[10]]
assignment6 <- apply(re_prob6,1,which.max)


table(assignment5,assignment4)

table(assignment6,assignment5)


soft54 <- t(re_prob5) %*% re_prob4
soft65 <- t(re_prob6) %*% re_prob5


# composition of each fit-5 class in terms of fit-4 classes
round(prop.table(soft54, margin = 1), 3)

# composition of each fit-6 class in terms of fit-5 classes
round(prop.table(soft65, margin = 1), 3)

# where each fit-4 class goes in the fit-5 solution
round(prop.table(soft54, margin = 2), 3)

# where each fit-5 class goes in the fit-6 solution
round(prop.table(soft65, margin = 2), 3)

