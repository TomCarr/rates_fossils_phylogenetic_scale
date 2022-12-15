library(devtools)
library(TreeSim)
library(phytools)
library(ggplot2)
###starting params###
#####################

fossilisation_rate <- 0.001
n_reps <- 50
n_extant_tips <-  100
speciation <- 0.5
extinction <- 0.45
bin <- 0.5

sim_fossils_probability_through_time <- function(fossilisation_rate, speciation, extinction, n_extant_tips, n_reps, bin){

###simulate the tree###
#######################
diversity_storage <- vector("list", n_reps)

for (a in 1:n_reps){
	tree <- sim.bd.taxa(n_extant_tips, 1, speciation, extinction, frac = 1, complete = TRUE, stochsampling = FALSE)[[1]]
	tree <- extract.clade(tree, findMRCA(tree, drop.tip(tree, getExtinct(tree))$tip.label, "node"))
	extant_tree <- drop.tip(tree, getExtinct(tree))
	tree_information <- cbind(tree[[1]], tree$edge.length) #get tree table with edge lengths
	tree_information <- tree_information[order(tree_information[,2]),] # order it by descendant node
	tree_information <- cbind(tree_information, c(max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[seq(1, length(tree$tip.label), 1)], (branching.times(tree)[-1] + (max(node.depth.edgelength(tree)) - branching.times(tree)[[1]])))) # add the branch end times
	tree_information <- cbind(tree_information, tree_information[,3]+tree_information[,4]) # add branch start times to the table 
	if (any(extant_tree[[1]][,2][which(extant_tree[[1]][,1] == length(extant_tree$tip.label) + 1)] <= length(extant_tree$tip.label)) == FALSE){
		extant_one <- extract.clade(tree, findMRCA(tree, extract.clade(extant_tree, extant_tree[[1]][,2][which(extant_tree[[1]][,1] == length(extant_tree$tip.label) + 1)][[1]])$tip.label, "node"))             
		extant_two <- extract.clade(tree, findMRCA(tree, extract.clade(extant_tree, extant_tree[[1]][,2][which(extant_tree[[1]][,1] == length(extant_tree$tip.label) + 1)][[2]])$tip.label, "node"))	
		extant_one_tree_information <- cbind(extant_one[[1]], extant_one$edge.length) #get tree table with edge lengths
		extant_one_tree_information <- extant_one_tree_information[order(extant_one_tree_information[,2]),] # order it by descendant node
		extant_one_tree_information <- cbind(extant_one_tree_information, c(max(node.depth.edgelength(extant_one)) - node.depth.edgelength(extant_one)[seq(1, length(extant_one$tip.label), 1)], (branching.times(extant_one)[-1] + (max(node.depth.edgelength(extant_one)) - branching.times(extant_one)[[1]])))) # add the branch end times
		extant_one_tree_information <- cbind(extant_one_tree_information, extant_one_tree_information[,3]+extant_one_tree_information[,4]) # add branch start times to the tree 
		extant_two_tree_information <- cbind(extant_two[[1]], extant_two$edge.length) #get tree table with edge lengths
		extant_two_tree_information <- extant_two_tree_information[order(extant_two_tree_information[,2]),] # order it by descendant node
		extant_two_tree_information <- cbind(extant_two_tree_information, c(max(node.depth.edgelength(extant_two)) - node.depth.edgelength(extant_two)[seq(1, length(extant_two$tip.label), 1)], (branching.times(extant_two)[-1] + (max(node.depth.edgelength(extant_two)) - branching.times(extant_two)[[1]])))) # add the branch end times
		extant_two_tree_information <- cbind(extant_two_tree_information, extant_two_tree_information[,3]+extant_two_tree_information[,4]) # add branch start times to the tree 
		extant_subtracter_information <- rbind(extant_one_tree_information, extant_two_tree_information) 
	} else {
		extant_one <- extract.clade(tree, findMRCA(tree, extract.clade(extant_tree, extant_tree[[1]][,2][which(extant_tree[[1]][,1] == length(extant_tree$tip.label) + 1)[which(extant_tree[[1]][,2][which(extant_tree[[1]][,1] == length(extant_tree$tip.label) + 1)] > length(extant_tree$tip.label))]])$tip.label, "node")) #!!!!!!
		extant_one_tree_information <- cbind(extant_one[[1]], extant_one$edge.length) #get tree table with edge lengths
		extant_one_tree_information <- extant_one_tree_information[order(extant_one_tree_information[,2]),] # order it by descendant node
		extant_one_tree_information <- cbind(extant_one_tree_information, c(max(node.depth.edgelength(extant_one)) - node.depth.edgelength(extant_one)[seq(1, length(extant_one$tip.label), 1)], (branching.times(extant_one)[-1] + (max(node.depth.edgelength(extant_one)) - branching.times(extant_one)[[1]])))) # add the branch end times
		extant_one_tree_information <- cbind(extant_one_tree_information, extant_one_tree_information[,3]+extant_one_tree_information[,4]) # add branch start times to the tree 
		extant_subtracter_information <- extant_one_tree_information
	}
	time_intervals <- seq(50, 0.5, -bin)
	counter <- vector(mode="numeric", length=0)
	subtracter <- vector(mode="numeric", length=0)
	for (b in 1:length(time_intervals)){
		less_thans <- which(tree_information[,4] <= time_intervals[[b]])
		more_thans <- which(tree_information[,5] > time_intervals[[b]])
		counter <- append(counter, length(which(more_thans %in% less_thans)))
	} 
	for (b in 1:length(time_intervals)){
		less_thans <- which(extant_subtracter_information[,4] <= time_intervals[[b]])
		more_thans <- which(extant_subtracter_information[,5] > time_intervals[[b]])
		subtracter <- append(subtracter, length(which(more_thans %in% less_thans)))	
	}		
	diversity_storage[[a]] <- append(diversity_storage[[a]], counter-subtracter)
}

###calculate fossilisation probabilities###
###########################################

probability <- vector("list", n_reps)
for (a in 1:length(diversity_storage)){
	probability_raw <- 1-(1-fossilisation_rate*bin)^diversity_storage[[a]]	
	probability_scaled <- vector(mode="numeric", length=0)
	probability_scaled[[1]] <- probability_raw[[1]]
		for (b in 2:length(diversity_storage[[a]])){
		probability_scaled <- append(probability_scaled, probability_raw[[b]]*(1-sum(probability_scaled[seq(1, b-1, 1)])))
		}
	probability[[a]] <- probability_scaled
}

mean_probability <- vector(mode="numeric", length=0)
for (a in 1:length(probability[[1]])){
	counter <- sapply(probability, "[[", a)
	mean_probability <- append(mean_probability, mean(counter))
}
}

##########
###PLOT###
##########

x_axis <- seq(50, 0.5, -bin)
table <- data.frame(x_axis, mean_probability)

plot <- ggplot(data=table, aes(x=table[,1], y=table[,2]))+
geom_point(aes(y=table[,2]), size=4.5) +
theme(panel.grid.major = element_line(colour = "grey", size=0.25), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(size=1,colour="black"), text=element_text(size=30,colour="black"), axis.ticks=element_line(size=1,colour="black"))






