re_prob <- to_save[[2]][[10]]


#for application
#just 2 percent of inviduals have a max prob < .99 
sum(apply(re_prob,1,max) < .99)/nrow(re_prob)

mean_entropy <- ifelse(re_prob > 0, re_prob * log(re_prob), 0) |>
  rowSums() |> 
  mean() * -1

mean_entropy

load("c:/Users/jorda/OneDrive/Documents/Research/NIH/JLC-HMM/Data/JMHMMNoSurvMix5Seed.rda")
re_prob5_ns <- to_save[[2]][[10]]
assignment5_ns <- apply(re_prob5_ns,1,which.max)

load("c:/Users/jorda/OneDrive/Documents/Research/NIH/JLC-HMM/Data/JMHMMMix5Seed.rda")
re_prob5 <- to_save[[2]][[10]]
assignment5 <- apply(re_prob5,1,which.max)


load("c:/Users/jorda/OneDrive/Documents/Research/NIH/JLC-HMM/Data/JMHMMMix4Seed.rda")
re_prob4 <- to_save[[2]][[10]]
assignment4 <- apply(re_prob4,1,which.max)


load("c:/Users/jorda/OneDrive/Documents/Research/NIH/JLC-HMM/Data/JMHMMMix6Seed.rda")
re_prob6 <- to_save[[2]][[10]]
assignment6 <- apply(re_prob6,1,which.max)


table(assignment5,assignment4)
#Moving from 4 to 5 classes, the added class, fit-5 LC5, is largely made up of subsets from fit-4 LC2 and LC3 groups 
#The core of fit-4 LC2 and LC3 remains preserved as fit-5 LC2 and LC3, but a nontrivial fraction of each is separated into the new LC5.



table(assignment6,assignment5)
#Moving from 5 to 6 classes, there is not as clear of an added LC
#The additional class appears to primarily split the previous LC1 group. The 6-class LC2 is composed mostly of individuals who were assigned to LC1 in the 5-class model, with a secondary contribution from the previous LC4. 
#Capturing additional heterogeneity along the boundary between the previous LC1 and LC4 groups
#Fit5 LC3 and LC5 are largely preserved as Fit6 LC3 and LC5, respectively,



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

