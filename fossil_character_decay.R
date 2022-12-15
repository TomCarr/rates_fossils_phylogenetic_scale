library(devtools)
library(TreeSim)
library(phytools)
library(phangorn)
library(ggplot2)

###starting params###
#####################

fossilisation_rate <- 2
speciation <- 3
extinction <- 0
morpho_rate <- 0.5

n_extant_tips <- 25
seq_length <- 400
n_reps <- 1

sim_fossils_on_tree <- function(fossilisation_rate, speciation, extinction, seq_length, morpho_rate, n_extant_tips, n_reps){

###storage_variables###
#######################

absolute_fossil_ages_storage <<- vector(mode="numeric", length=0)
absolute_node_ages_storage <<- vector(mode="numeric", length=0)

###RUN THE SIMULATION###
########################

for (a in 1:n_reps){

	message(paste("rep", a, sep=" "))
	
###simulate tree
	tree <- sim.bd.taxa(n_extant_tips, 1, speciation, extinction, frac = 1, complete = TRUE, stochsampling = FALSE)[[1]] # the overall tree simulated from the origin
	crown_group_tree <- extract.clade(tree, findMRCA(tree, drop.tip(tree, getExtinct(tree))$tip.label, "node")) # the tree starting from the mrca of extant species
	fossil_addition_tree <- crown_group_tree # the tree on which fossil tips are sequentially added

###simulate fossil times
	table <- cbind(crown_group_tree[[1]], crown_group_tree$edge.length) # node information in crown group tree with edge lengths
	fossil_vector <- vector("list", nrow(table))
	fossil_names <- vector("list", nrow(table))
	for (b in 1:nrow(table)){ # fossils simulated for each branch in the tree
		continue = "yes" 
		fossil_vector[[b]] <- append(fossil_vector[[b]], 0) # initiate the vector
		while (continue == "yes"){
			fossil_vector[[b]] <- append(fossil_vector[[b]], fossil_vector[[b]][[length(fossil_vector[[b]])]] + rexp(1, fossilisation_rate)) # add next fossil sampling time according to a given rate
			if (fossil_vector[[b]][[length(fossil_vector[[b]])]] >= table[,3][[b]]){ # if the new fossil sampling time is greater than the length of the relevent branch you stop adding more fossils
				continue = "no"
			}
		}
		fossil_vector[[b]] <- fossil_vector[[b]][-which(fossil_vector[[b]] == 0)] # get rid of vector initiator
		fossil_vector[[b]] <- fossil_vector[[b]][-which(fossil_vector[[b]] > table[,3][[b]])] # get rid of over lengthed fossils
		if (length(fossil_vector[[b]]) > 0){ # if there are fossils on a branch, create names for them and then add them
			for (c in 1:length(fossil_vector[[b]])){
				fossil_names[[b]] <- append(fossil_names[[b]], paste("fossil_", b, "_", c, sep=""))
			} 

###add fossils to tree
			if (table[,2][[b]] > length(crown_group_tree$tip.label)){ # if its not a terminal branch you are adding the fossil tip to
				for (c in 1:length(fossil_vector[[b]])){
					fossil_addition_tree <- bind.tip(fossil_addition_tree, fossil_names[[b]][[c]], edge.length=0, where=findMRCA(fossil_addition_tree, extract.clade(crown_group_tree, table[,2][[b]])$tip.label, "node"), position=table[,3][[b]] - fossil_vector[[b]][[c]])
				}
			}
			if (table[,2][[b]] <= length(crown_group_tree$tip.label)){
				for (c in 1:length(fossil_vector[[b]])){
					fossil_addition_tree <- bind.tip(fossil_addition_tree, fossil_names[[b]][[c]], edge.length=0, where=which(fossil_addition_tree$tip.label == crown_group_tree$tip.label[[table[,2][[b]]]]), position=table[,3][[b]] - fossil_vector[[b]][[c]])
				}
			}
		}
	}
	message(paste("simulated a ", n_extant_tips, " tip tree with ", length(unlist(fossil_names)), " fossils", sep=""))
###simulate sequences###
########################
	message("now simulating sequences")
	sequences <- simSeq(fossil_addition_tree, seq_length, type = "USER", levels = c(0,1), rate = morpho_rate) # simulate morphological characters on the fossil addition tree

	extant_species <- drop.tip(crown_group_tree, getExtinct(crown_group_tree))$tip.label # get extant tips in the crown group tree
	extant_tree <- keep.tip(crown_group_tree, extant_species) # get the extant only crown group tree 

	extant_sequences <- subset(sequences, extant_tree$tip.label) # just the sequences of extant tips in crown group tree
	fossil_sequences <- subset(sequences, unlist(fossil_names)) # just the sequences of the fossils

###search for synapomorphies in extant_tree###
##############################################
	if (length(fossil_sequences) > 0){	
	
		message("now getting clade synapomorphies in extant tree")
	
		fossil_in_clade <- vector("list", length(fossil_sequences))
		for (b in which(extant_tree[[1]][,2] > length(extant_tree$tip.label))){ # going through all nodes, excluding terminals, and the root node
			clade_tips <- extract.clade(extant_tree, extant_tree[[1]][,2][[b]])$tip.label # get tips of the relevent clade
			non_clade_tips <- extant_tree$tip.label[-which(extant_tree$tip.label %in% clade_tips)] # get tips of everything else
			clade_sequences <- subset(extant_sequences, clade_tips) # get sequences of the clade
			non_clade_sequences <- subset(extant_sequences, non_clade_tips) # get sequences of everything else
			clade_synapomorphies <- vector("list", length=0)
			for (c in 1:seq_length){ # loop progressing along sites, going to search for synapomorphies for clade at each site 
				if (length(unique(as.numeric(unlist(subset(extant_sequences, select = c, site.pattern=FALSE))))) != 1){ #if its not an invariant site with respect to all extant sequences
					clade_site_states <- as.numeric(unlist(subset(clade_sequences, select = c, site.pattern=FALSE))) # get states at that site within the clade of interest
					most_common_at_site <- sort(table(clade_site_states),decreasing=TRUE)[[1]] # get the frequency of the most common state at that site within the clade of interest
					if (most_common_at_site > 0.9*length(clade_tips)){ #if the highest frequency state is in 90% of the tips...
						pre_synapomorphy <- as.numeric(names(sort(table(clade_site_states),decreasing=TRUE)))[[1]] # assign the state to pre-synapomorphy status
						non_clade_site_states <- as.numeric(unlist(subset(non_clade_sequences, select = c, site.pattern=FALSE))) # get states at the site outside the clade of interest 
						if ((pre_synapomorphy %in% non_clade_site_states) == FALSE){ # see if its outside the clade too, if not we call it a synapomorphy
							clade_synapomorphies[[length(clade_synapomorphies)+1]] <- c(c, pre_synapomorphy) # get a lits of all the synapomorphies for the relevent clade
						}
					}
				}
			}


###find which fossils can go to the clade###
############################################
		
			if (length(clade_synapomorphies) > 0){
				
				message(paste("assigining fossils to clade", extant_tree[[1]][,2][[b]], sep=""))
				
				for (c in 1:length(unlist(fossil_names))){
					test_sequence <- subset(sequences, unlist(fossil_names)[[c]]) # get the sequence for the fossil we are testing
					fossil_site_vector <- vector(mode="numeric", length=0)
					for (d in 1:length(clade_synapomorphies)){# have to do with this loop becaus of site pattern problems
						fossil_site_vector <- append(fossil_site_vector, as.numeric(unlist(subset(test_sequence, select = clade_synapomorphies[[d]][[1]], site.pattern=FALSE))))
					}
					if (sum(fossil_site_vector == sapply(clade_synapomorphies, "[[", 2)) > 0.9*length(clade_synapomorphies)){ # if the fossil has 90% of the synapomorphies
						fossil_in_clade[[c]] <- append(fossil_in_clade[[c]], extant_tree[[1]][,2][[b]]) #then the fossil belongs to the clade we are investigating
					}
				}
			}
		}
		if (length(unlist(fossil_in_clade)) > 0){

###now get oldest fossil in each clade###
#########################################

			message("getting oldest fossil in each clade")
			
			removal <- which(lengths(fossil_in_clade) == 0)
			if (length(removal) > 0){
				fossil_in_clade <- sapply(fossil_in_clade[-removal], max)
				final_fossil_names <- unlist(fossil_names)[-removal]
			} else {			
				fossil_in_clade <- sapply(fossil_in_clade, max)
				final_fossil_names <- unlist(fossil_names)
			}
			
			absolute_fossil_ages <- vector(mode="numeric", length=0) # get the absolute age of each fossil that has been assigned to a clade in the tree
			for (b in 1:length(final_fossil_names)){
				absolute_fossil_ages <- append(absolute_fossil_ages, max(node.depth.edgelength(fossil_addition_tree)) - node.depth.edgelength(fossil_addition_tree)[which(fossil_addition_tree$tip.label == final_fossil_names[[b]])]) #get the absolute age of each of the fossils
			}

				fossil_in_clade_removal <- vector(mode="numeric", length=0) # this removes fossils if there is an older fossil at the same node or in descendants
			for (b in 1:length(fossil_in_clade)){
				for (c in 1:length(fossil_in_clade[-b])){
					if ((fossil_in_clade[[c]] == fossil_in_clade[[b]]) & (absolute_fossil_ages[[c]] > absolute_fossil_ages[[b]])){
						fossil_in_clade_removal <- append(fossil_in_clade_removal, b)
					} else if ((fossil_in_clade[[c]] %in% getDescendants(extant_tree, fossil_in_clade[[b]])) & (absolute_fossil_ages[[c]] >= absolute_fossil_ages[[b]])){
						fossil_in_clade_removal <- append(fossil_in_clade_removal, b)
					}
				}
			}

			final_fossil_names <- final_fossil_names[-unique(fossil_in_clade_removal)]
			fossil_in_clade <- fossil_in_clade[-unique(fossil_in_clade_removal)]
			absolute_fossil_ages <- absolute_fossil_ages[-unique(fossil_in_clade_removal)]

			node_ages <- vector(mode="numeric", length=0)
			for (b in 1:length(absolute_fossil_ages)){
				node_ages <- append(node_ages, max(node.depth.edgelength(extract.clade(extant_tree, Ancestors(extant_tree, fossil_in_clade[[b]])[[1]])))) 
			}
			absolute_fossil_ages_storage <<- append(absolute_fossil_ages_storage, absolute_fossil_ages)
			absolute_node_ages_storage <<- append(absolute_node_ages_storage, node_ages)
		}	
	}
}

}

percentage_difference <- (absolute_fossil_ages_storage-absolute_node_ages_storage)/absolute_node_ages_storage

table <- data.frame(absolute_node_ages_storage, percentage_difference)

plot <- ggplot(data=table, aes(x=table[,1], y=table[,2]))+
geom_point(aes(y=table[,2]), size=4.5) +
theme(panel.grid.major = element_line(colour = "grey", size=0.25), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(size=1,colour="black"), text=element_text(size=30,colour="black"), axis.ticks=element_line(size=1,colour="black"))







