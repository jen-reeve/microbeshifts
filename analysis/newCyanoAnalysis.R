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
setwd("~/repos/microbeshifts/analysis")
source("./shiftfunctions.R")
source("./cyanofunctions.R")

## Read in posterior distribution of trees

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

tree <- lbaTrees[[1]][[1]]
fns <- make.bisse.fns(tree, dat)
registerDoParallel(cores=8)
fits.ARDnotime <- fitFns(5, fns$ARD.notime, fns$tdList, res=NULL)


#seq <- seq(0.2, 3.6, 0.5); Fns <- fns; fns <- Fns$ARD.R2; starts<-fits.ARDnotime; res=NULL; tds <- fns$tdList
prof.ARDR2 <- profiles(1, fns$ARD.R2, fns$tdList, fits.ARDnotime, res=NULL, seq=seq(0.1, 3.6, length.out=32), cores=8)

lnL.ARDR2 <- lapply(prof.ARDR2, function(x) sapply(x, function(y) max(y$value)))
lnL.ARDR2 <- do.call(rbind, lnL.ARDR2)
cslnL.ARDR2 <- apply(lnL.ARDR2, 2, cumsum)
plot(seq(3.6, 0.1, length.out=32), rev(cslnL.ARDR2[25,]))



