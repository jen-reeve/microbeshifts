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
require(plotrix)
require(grid)
require(bayou)
#setwd("~/repos/microbeshifts/analysis")
source("./shiftfunctions.R")
source("./cyanofunctions.R")

## Read in posterior distribution of trees
goodFiles <- list.files("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/")
goodTrees <- lapply(goodFiles[grep("good", goodFiles)], function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/", x, sep="")))

youngFiles <- list.files("../data/Posterior Probability Sampling of Trees/Younger Root MaxAge/")
youngTrees <- lapply(youngFiles,function(x) read.trees.from.mesquite(paste("../data/Posterior Probability Sampling of Trees/Younger Root MaxAge/", x, sep="")))

dat <- clean.data()

trees <- c(goodTrees[[1]][seq(1,100,length.out=5)],goodTrees[[2]][seq(1,100,length.out=5)],goodTrees[[3]][seq(1,100,length.out=5)],goodTrees[[4]][seq(1,100,length.out=5)], youngTrees[[1]][seq(1,100, length.out=10)], youngTrees[[2]][seq(1,100, length.out=10)])
trees <- lapply(trees, function(x){x$edge.length <- x$edge.length/1000; x})
TLs <- sapply(1:length(trees), function(x) max(branching.times(trees[[x]])))
rm(youngFiles)
rm(goodFiles)
rm(goodTrees)
rm(youngTrees)
tree <- trees[[1]]
fns <- make.bisse.fns(tree, dat)


Fns <- lapply(trees, function(x) make.bisse.fns(x, dat))
registerDoParallel(cores=5)
#fits.ARDnotime <- lapply(Fns, function(fns) fitFns(20, fns$ARD.notime, fns$tdList, res=NULL))
#saveRDS(fits.ARDnotime, file="../output/fits.ARDnotime.rds")
fits.ARDnotime <- readRDS("../output/fits.ARDnotime.rds")
allfits <- do.call(rbind, fits.ARDnotime)
nsims <- 100 

Qs <- lapply(sample(1:nrow(allfits), nsims, replace=FALSE), function(x){rates <- allfits[x, 1:2]})#; Q <- matrix(c(-1*rates[1], rates[1], rates[2], -1*rates[2]), ncol=2, nrow=2, byrow=TRUE); Q})
Qs <- lapply(Qs, function(x) do.call(c, ifelse(x<0.01/100, 0.01/100, x)))
Qs <- lapply(Qs, function(x) ifelse(x>0.5*100, 0.5*100, x))

lbaFiles <- list.files("../data/Posterior Probability Sampling of Trees/Trees Exhibiting Long Branch Attraction/Posterior probability distribution of time trees/")
lbaTrees <- lapply(lbaFiles,function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Exhibiting Long Branch Attraction/Posterior probability distribution of time trees/", x, sep="")))
goodFiles <- list.files("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/")
goodTrees <- lapply(goodFiles[grep("good", goodFiles)], function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/", x, sep="")))
youngFiles <- list.files("../data/Posterior Probability Sampling of Trees/Younger Root MaxAge/")
youngTrees <- lapply(youngFiles,function(x) read.trees.from.mesquite(paste("../data/Posterior Probability Sampling of Trees/Younger Root MaxAge/", x, sep="")))

trees <- do.call(c, list(do.call(c, lbaTrees), do.call(c, goodTrees), do.call(c, youngTrees)))
treesources <- do.call(c, list(
  do.call(c, lapply(1:length(lbaTrees), function(x) rep(paste("lba", x, sep=""), length(lbaTrees[[x]])))), 
  do.call(c, lapply(1:length(goodTrees), function(x) rep(paste("good", x, sep=""), length(goodTrees[[x]])))), 
  do.call(c, lapply(1:length(youngTrees), function(x) rep(paste("young", x, sep=""), length(youngTrees[[x]]))))
))
rm(Fns, lbaTrees, lbaFiles, goodFiles, goodTrees, youngFiles, youngTrees)
gc()
simtreeids <- lapply(1:nsims, function(x) sample(1:length(trees), 2, replace=FALSE))
simtreeids <- do.call(rbind, simtreeids)
simtreesources <- apply(simtreeids, 2, function(x) treesources[x])
simtreeTLs <- cbind(sapply(trees[simtreeids[,1]], function(x) max(branching.times(x))),sapply(trees[simtreeids[,2]], function(x) max(branching.times(x))))

## Simulations on one posterior tree, estimation on another
trees <- lapply(trees, function(x){x$edge.length <- x$edge.length/1000; x})
TLs <- lapply(1:nrow(simtreeids), function(x) max(branching.times(trees[[simtreeids[x,2]]])))

fits <- list()
true.lnL <- list()
true.pars <- list()
treeids <- list()
traits <- list()
for(i in 1:nsims){
  treeids[[i]] <- simtreeids[i,]
  trait <- 1
  while(length(unique(trait))<=1){
    trait <- diversitree::sim.character(trees[[treeids[[i]][1]]], pars=c(Qs[[i]][1], Qs[[i]][2]), x0=sample(c(0,1), 1), model = "mk2")
  }
  traits[[i]] <- trait
  strait <- trait
  shuffled <- sample(1:length(trait), 10, replace=FALSE)
  strait[shuffled] <- trait[shuffled][sample(shuffled, 10, replace=FALSE)]
  fns <- make.bisse.fns(trees[[treeids[[i]][1]]], data.frame("trait"=strait))
  fits.ARDnotime <- fitFns(10, fns$ARD.notime, fns$tdList, res=NULL)
  true.lnL[[i]] <- fits.ARDnotime$lnL
  #fits[[i]] <- profiles(1, fns$ARD.R2, fns$tdList, starts=cbind(1, 1, Qs[[i]][[1,2]], Qs[[i]][[2,1]]), res=NULL, seq=seq(0.1, TLs[[i]]-0.1, length.out=37), cores=8)[[1]]
  startx <- do.call(rbind, lapply(1:37, function(x) cbind(1, 1, Qs[[i]][1], Qs[[i]][2])))
  fits[[i]] <- profiles(1, fns$ARD.R2, fns$tdList, starts=list(startx), res=NULL, seq=seq(0.1, TLs[[i]]-0.1, length.out=37), cores=8, start4=TRUE) 
}
lnLs <- unlist(sapply(lnLs,function(x) as.numeric(as.character(x))))
profs <- lapply(fits, function(x) do.call(rbind, lapply(x[[1]], function(z) z[which(z$value==max(z$value)), c(5, 1:4)])))
profs2 <- lapply(1:length(profs), function(x) profs[[x]][,1]-lnLs[x])
