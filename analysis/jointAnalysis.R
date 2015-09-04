## # Single multiplier & time-shift point
## Full New Cyanos Analysis over trees
## # Evolution of oxygen metabolism in freshwater and marine habitats in the Archaea
## Load packages and prepare environment
require(phytools)
require(diversitree)
require(geiger)
require(devtools)
require(rncl)
require(foreach)
require(doParallel)
require(aRbor)
require(optimx)
setwd("~/repos/microbeshifts/analysis")
source("./shiftfunctions.R")
source("./cyanofunctions.R")

## Read in posterior distribution of trees

youngFiles <- list.files("../data/Posterior Probability Sampling of Trees/Younger Root MaxAge/")
youngTrees <- lapply(youngFiles,function(x) read.trees.from.mesquite(paste("../data/Posterior Probability Sampling of Trees/Younger Root MaxAge/", x, sep="")))

#lbaFiles <- list.files("../data/Posterior Probability Sampling of Trees/Trees Exhibiting Long Branch Attraction/Posterior probability distribution of time trees/")
#lbaTrees <- lapply(lbaFiles,function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Exhibiting Long Branch Attraction/Posterior probability distribution of time trees/", x, sep="")))

#goodFiles <- list.files("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/")
#goodTrees <- lapply(goodFiles[grep("good", goodFiles)], function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/", x, sep="")))
#cyanosTrees <- lapply(goodFiles[grep("Cyanos", goodFiles)], function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/", x, sep="")))

## For Josef: Figure out how these trees all differ
#lba.ntips <- lapply(lbaTrees, function(x) unname(sapply(x, function(y) length(y$tip.label))))
#good.ntips <- lapply(goodTrees, function(x) unname(sapply(x, function(y) length(y$tip.label))))
#cyanos.ntips <- lapply(cyanosTrees, function(x) unname(sapply(x, function(y) length(y$tip.label))))

#lba.species <- unique(unlist(sapply(lbaTrees, function(x) x[[1]]$tip.label)))
#good.species <- unique(unlist(sapply(goodTrees, function(x) x[[1]]$tip.label)))
#cyanos.species <- unique(unlist(sapply(cyanosTrees, function(x) x[[1]]$tip.label)))
#all.species <- unique(c(lba.species, good.species, cyanos.species))


## # Tasks:
## 1. Create new profile plots over the posterior of trees, so how much uncertainty this introduces
## 2. Simulate no shift models on on set of trees, estimate shifts on another. See if artifactual shifts are produced.
## 3. Manually shuffle trees and see if artifactual shifts are produced. 

## Create profile plots over the posterior of trees
## First create a function that matches a given tree to the data and creates the set of bisse functions we want
dat <- clean.data()
dat1 <- select(dat, Thermophilic, Akinetes, Heterocysts, Nonfreshwater_habitat, Nitrogen_fixation, Morphology, Habit, Freeliving, Mats, Epi_Endolithic, Epiphytic,
              Periphytic, Motility, Hormogonia, Gas_vesicles, False_Branching, True_Branching, Fission_in_multiple_planes, Uniseriate_trichome, 
              Multiseriate_trichomes, Baeocytes, Extracellular_sheath, Mucilage, celldiam_mean)
trees <- c(youngTrees[[1]][seq(1,100, length.out=10)], youngTrees[[2]][seq(1,100, length.out=10)])
trees <- lapply(trees, function(x){x$edge.length <- x$edge.length/1000; x})
TLs <- sapply(1:length(trees), function(x) max(branching.times(trees[[x]])))

Fns <- lapply(trees, function(x) make.bisse.fns(x, dat1))
registerDoParallel(cores=5)
#fits.ARDnotime <- lapply(Fns, function(fns) fitFns(5, fns$ARD.notime, fns$tdList, res=NULL))
#saveRDS(fits.ARDnotime, file="../output/fits.ARDnotime.joint.rds")
fits.ARDnotime <- readRDS("../output/fits.ARDnotime.joint.rds")
starts <- lapply(fits.ARDnotime, function(x) x[,1:2])

fns <- lapply(1:length(trees), function(x) Fns[[x]]$ARD.R1)
rm(Fns)

## Fit model with a single, estimated time shift and a single rate multiplier
## Define joint likelihood function
#Parameter order c(t, r, q01.1, q1.10.1, ..., q01.22, q10.22); total of 46 parameters
ntrait <- ncol(dat1)
startx <-lapply(1:length(fns), function(x) c(runif(1,0.3,2.3), 0, unlist(do.call(c, lapply(1:nrow(starts[[x]]), function(y) starts[[x]][y,])))))
liks <- list()
make.liks <- function(i) {
  fx <- fns[[i]]
  lik <- function(pars){
    time <- pars[1]
    r <- exp(pars[2])
    qs <- (matrix(pars[3:length(pars)], ncol=2, byrow=TRUE))
    pars2 <- cbind(r, qs, time)
    lnLs <- sapply(1:length(fx), function(x) fx[[x]](pars2[x,]))
    res <- -1*sum(lnLs)
    return(res)
    }
  return(lik)
}

liks <- lapply(1:length(fns), function(x) make.liks(x))

## Set function attributes
#liks <- lapply(liks, function(x){attributes(x)$varnames <- c("t", "r", paste(c("q01", "q10"), sort(rep(1:ntrait,2)), sep="."));
#  attributes(x)$converter <- function(pars){
#    pars[2:length(pars)] <- exp(pars[2:length(pars)])
#    names(pars) <- attributes(x)$varnames
#    pars
#  }; attributes(x)$starter <- function(){
#      c(runif(1,0.5, min(TLs)-0.5), rnorm(2*ntrait+1, 0, 0.5))
#    }; attributes(x)$lower <- c(0, log(0.001), log(rep(0.01/100, ntrait*2)));
#    attributes(x)$upper <- c(max(TLs), log(1000000), log(rep(0.5*100, ntrait*2))); x})

## Optimize joint likelihood functions
for(i in 1:length(liks)){
  tmp <- optimx(startx[[i]], liks[[i]], lower=c(0, log(0.001), rep(0, ntrait*2)), upper=c(max(TLs), log(1000000),  rep(0.5*100, ntrait*2)), control=list(all.methods=TRUE))
  saveRDS(tmp, file=paste("../output/jointoptim",i,".rds",sep=""))
} 

#fits.ARDnotime[[i]]$lnL-lnLs

#tmp <- find.mle(liks[[i]], x.init=attributes(liks[[i]])$starter(), method="optim",  lower=attributes(liks[[i]])$lower, upper=attributes(liks[[i]])$upper)

#res2 <- foreach(i=1:5) %dopar% find.mle(fn, x.init=attributes(fn)$smartstart(0.5), method="optim", lower=attributes(fn)$lower, upper=attributes(fn)$upper)
#res.liks <- sapply(res, function(x) x$lnLik)
#best.id <- which(res.liks == max(res.liks))
#best.res <- res[[best.id]]

