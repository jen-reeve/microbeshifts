## Full New Cyanos Analysis over trees
## # Evolution of oxygen metabolism in freshwater and marine habitats in the Archaea
## Load packages and prepare environment
require(phytools)
require(diversitree)
require(geiger)
require(devtools)
require(testthat)
require(roxygen2)
require(rncl)
require(foreach)
require(doParallel)
require(aRbor)
require(optimx)
#setwd("~/repos/microbeshifts/analysis")
source("./shiftfunctions.R")
source("./cyanofunctions.R")

## Read in posterior distribution of trees
tree1 <- read.tree("../data/Tree1MrB.tre")

lbaFiles <- list.files("../data/Posterior Probability Sampling of Trees/Trees Exhibiting Long Branch Attraction/Posterior probability distribution of time trees/")
lbaTrees <- lapply(lbaFiles,function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Exhibiting Long Branch Attraction/Posterior probability distribution of time trees/", x, sep="")))

goodFiles <- list.files("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/")
goodTrees <- lapply(goodFiles[grep("good", goodFiles)], function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/", x, sep="")))
cyanosTrees <- lapply(goodFiles[grep("Cyanos", goodFiles)], function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/", x, sep="")))

## For Josef: Figure out how these trees all differ
lba.ntips <- lapply(lbaTrees, function(x) unname(sapply(x, function(y) length(y$tip.label))))
good.ntips <- lapply(goodTrees, function(x) unname(sapply(x, function(y) length(y$tip.label))))
cyanos.ntips <- lapply(cyanosTrees, function(x) unname(sapply(x, function(y) length(y$tip.label))))

lba.species <- unique(unlist(sapply(lbaTrees, function(x) x[[1]]$tip.label)))
good.species <- unique(unlist(sapply(goodTrees, function(x) x[[1]]$tip.label)))
cyanos.species <- unique(unlist(sapply(cyanosTrees, function(x) x[[1]]$tip.label)))
all.species <- unique(c(lba.species, good.species, cyanos.species))


## # Tasks:
## 1. Create new profile plots over the posterior of trees, so how much uncertainty this introduces
## 2. Simulate no shift models on on set of trees, estimate shifts on another. See if artifactual shifts are produced.
## 3. Manually shuffle trees and see if artifactual shifts are produced. 

## Create profile plots over the posterior of trees
## First create a function that matches a given tree to the data and creates the set of bisse functions we want
dat <- clean.data()
tree <- tree1
#trees <- list(goodTrees[[1]][1:2], goodTrees[[2]][1:2], goodTrees[[3]][1:2], goodTrees[[4]][1:2])


