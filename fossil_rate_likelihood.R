library(devtools)
library(TreeSim)
library(phytools)
library(phangorn)
library(ggplot2)
library(DirichletReg)
library(combn)

num_fossils_to_test <- 2 # number of fossils in the analysed dataset

node_height <- vector("list", num_fossils_to_test) # initial vector for the height of the nodes to which different fossils belong
diversity <- vector("list", num_fossils_to_test) # the diversity of each of the fossils at each of the time points
window_size <- vector("list", num_fossils_to_test) # the time interval within which the fossilisation probability is calculated
time_points <- vector("list", num_fossils_to_test) # the time slices to which diversity refers to for each of the fossils
oldest_fossil <- vector("list", num_fossils_to_test) # the age of eahc of the fossils we are analysing 

node_height[[1]] <- 1.726675
node_height[[2]] <- 1.847709
diversity[[1]] <- c(2, 2, 4, 5, 10, 11, 12, 14, 22, 32, 45, 61, 81, 96, 119, 146, 175, 223, 291, 374, 500)
diversity[[2]] <- c(2, 2, 2, 5, 6, 7, 11, 14, 17, 20, 27, 38, 50, 62, 92, 122, 150, 218, 293, 369, 500)
window_size[[1]] <- 0.05*node_height[[1]]
window_size[[2]] <- 0.05*node_height[[2]]
time_points[[1]] <- seq(node_height[[1]], 0, -window_size[[1]])
time_points[[2]] <- seq(node_height[[2]], 0, -window_size[[2]])
oldest_fossil[[1]] <- 0.8270752
oldest_fossil[[2]] <- 1.631107

rate_number_likelihood_scores <- vector("list", length(num_fossils_to_test)) #store likelihood scores for different numbers of rate regimes, with the number of regimes determined by how man fossils there are 
rate_number_aic_scores <- vector("list", length(num_fossils_to_test)) #same as above, but aic score with respect to best scoring model with fewer shifts
rate_estimates <- vector("list", length(num_fossils_to_test))

keep_testing <- "yes" # start search through models with increasing complexity
m=1

while (keep_testing == "yes"){
###get rate maps for a given complexity - so for a given number of rates we create a map of rates
mapper <- vector("list", 0)
storer <- vector(mode="numeric", length=0)
for (b in 1:50){
test <- sample(m, num_fossils_to_test, repl = TRUE)
if (length(unique(test)) == m){
if (paste(test, collapse="") %in% storer == FALSE){
mapper[[length(mapper)+1]] <- test
storer <- append(storer, paste(test, collapse=""))
}
}
}

###initiate the vectors to store the likelihood and aic scores for a given number of shifts
rate_number_likelihood_scores[[m]] <- vector(mode="numeric", length=0)
rate_number_aic_scores[[m]] <- vector(mode="numeric", length=0)
rate_estimates[[m]] <- vector(mode="numeric", length=0)

for (n in 1:length(mapper)){ # run through the map for a given number of rate shifts 

###generate_likelihood_equation###
##################################

### 1) probability that nothing happens = (1-pexp(window_size, fossil_rate)) ^ diversity
### 2) probability that something happens = 1 - (1-pexp(window_size, fossil_rate)) ^ diversity         ###prob of sampling earlier              ###prob of not sampling in time window if you havent sampled earlier
### 3) probability that nothing happens weighted by probability something may already have happened (i.e. more likely not to happen if already happened) = (1 - (1-pexp(window_size, fossil_rate)) ^ all_of_past_diversity) + (((1-pexp(window_size, fossil_rate)) ^ diversity) * (1-pexp(window_size, fossil_rate)) ^ all_of_past_diversity) # the two things are mutually exclusive
### 4) probability that something happend weighted by probability its not already happened = (1 - (1-pexp(window_size, fossil_rate)) ^ diversity) * ((1-pexp(window_size, fossil_rate)) ^ all_of_past_diversity)

##################################

rate_map_to_use <- mapper[[n]] # were are going to use one of the sets of mapped rates

probability_string <- "(" # start generating the parbability function
for (a in 1:(length(oldest_fossil))){ # we go through each one of the fossils
previous_diversity <- vector(mode="numeric", length=0)
for (b in 2:(length(time_points[[a]]))){ # and for each of the fossils we do each of the time points
if ((oldest_fossil[[a]] >= time_points[[a]][[b]]) & (oldest_fossil[[a]] < time_points[[a]][[b-1]])){ # if the fossil age is bigger than the current time point but less than the one before it
probability_string <- append(probability_string, paste("((1-(1-pexp(window_size[[a]], fossil_rate[", rate_map_to_use[[a]], "]))^", ((diversity[[a]][[b]]+diversity[[a]][[b-1]])/2), ") * ((1-pexp(window_size[[a]], fossil_rate[", rate_map_to_use[[a]], "]))^",sum(previous_diversity),"))*", sep="")) 
} else {
probability_string <- append(probability_string, paste("((1-(1-pexp(window_size[[a]], fossil_rate[", rate_map_to_use[[a]], "])) ^", sum(previous_diversity), ") +  (((1-pexp(window_size[[a]], fossil_rate[", rate_map_to_use[[a]], "])) ^", ((diversity[[a]][[b]]+diversity[[a]][[b-1]])/2), ") * ((1-pexp(window_size[[a]], fossil_rate[", rate_map_to_use[[a]], "])) ^", sum(previous_diversity), ")))*", sep=""))
}
previous_diversity <- append(previous_diversity, (diversity[[a]][[b]]+diversity[[a]][[b-1]])/2)
}
}

changer <- strsplit(probability_string[[length(probability_string)]], split="")
changer[[1]][[length(changer[[1]])]] <- ")"
probability_string[[length(probability_string)]] <- paste(changer[[1]], collapse="")

final_prob_calc <- paste(probability_string, collapse="")

likelihood <- function(fossil_rate){
eval(parse(text=final_prob_calc))*-1
}
likelihood_calculation <- nlm(likelihood,rep(0.2, length(unique(rate_map_to_use))))
rate_number_likelihood_scores[[m]] <- append(rate_number_likelihood_scores[[m]], likelihood_calculation$minimum)
rate_number_aic_scores[[m]] <- append(rate_number_aic_scores[[m]], (2*m - 2*likelihood_calculation$minimum) + ((2*m*(m+1))/(num_fossils_to_test-m-1)))
rate_estimates[[m]] <- append(rate_estimates[[m]], likelihood_calculation$estimate)
}
### conditions on whether or not to continue with the search 
if (m > 1){
if (m < num_fossils_to_test){
if (min(rate_number_aic_scores[[m]]) < min(rate_number_aic_scores[[m-1]])){
m=m+1
}
} else if (m == num_fossils_to_test){
if (min(rate_number_aic_scores[[m]]) < min(rate_number_aic_scores[[m-1]])){
message(paste("best model has ", m, " rates", sep=""))
keep_testing <- "no"
}
} else if (min(rate_number_aic_scores[[m]]) >= min(rate_number_aic_scores[[m-1]])){
message(paste("best model has ", m-1, " rates", sep=""))
keep_testing <- "no"
}
} else if (m == 1){
if (m != num_fossils_to_test){
m=m+1
} else if (m == num_fossils_to_test){
message(paste("best model has ", m, " rates", sep=""))
keep_testing <- "no"
}
}
}



