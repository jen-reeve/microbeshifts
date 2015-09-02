## Simulation of Mk2 model cyanobacteria
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
goodFiles <- list.files("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/")
goodTrees <- lapply(goodFiles[grep("good", goodFiles)], function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/", x, sep="")))

dat <- clean.data()

tree <- goodTrees[[1]][[1]]
fns <- make.bisse.fns(tree, dat)
registerDoParallel(cores=8)

fns$ARD.notime[[1]]
Q <- matrix(c(-0.03, 0.03, 
        0.23, -0.23), nrow=2, ncol=2, byrow=TRUE)

## Simulations on one posterior tree, estimation on another
trees <- unlist(goodTrees, recursive=FALSE, use.names = FALSE)
#trees <- lapply(trees, function(x){x$edge.length <- x$edge.length/1000; x})
nsims <- 100
fits <- list()
true.lnL <- list()
true.pars <- list()
treeids <- list()
for(i in 1:nsims){
  treeids[[i]] <- sample(1:length(trees), 2, replace=FALSE)
  trait <- diversitree::sim.character(trees[[treeids[[i]][1]]], pars=true.pars[[i]], x0=sample(c(0,1), 1), model = "mk2")
  fns <- make.bisse.fns(trees[[treeids[[i]][2]]], data.frame("trait"=trait))
  fits.ARDnotime <- fitFns(5, fns$ARD.notime, fns$tdList, res=NULL)
  true.lnL[[i]] <- fits.ARDnotime$lnL
  fits[[i]] <- profiles(1, fns$ARD.R2, fns$tdList, cbind(0.03, 0.23), res=NULL, seq=seq(0.1, 3.6, 0.1), cores=8)[[1]]
}

results <- list(fits=fits, lnL=true.lnL, treeids=treeids, true.pars=true.pars)
saveRDS(results, file="../output/simulation_study1.rds")

res <- lapply(fits, function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
res <- do.call(cbind, res)
res <- apply(res, 2, function(x) x-min(x))
res1 <- res


## Tip shuffling on the phylogeny
nsims <- 100
fits <- list()
true.lnL <- list()
true.pars <- list()
TL <- list()
treeid <- list()
for(i in 1:nsims){
  treeid[[i]] <- sample(1:length(trees), 1)
  true.pars[[i]] <- runif(2, 0.1, 1)
  trait <- diversitree::sim.character(trees[[treeid[[i]]]], pars=true.pars[[i]], x0=sample(c(0,1), 1), model = "mk2")
  strait <- trait
  shuffled <- sample(1:length(trait), 10, replace=FALSE)
  strait[shuffled] <- trait[shuffled][sample(shuffled, 10, replace=FALSE)]
  fns <- make.bisse.fns(trees[[treeid[[i]]]], data.frame("trait"=strait))
  fits.ARDnotime <- fitFns(5, fns$ARD.notime, fns$tdList, res=NULL)
  true.lnL[[i]] <- fits.ARDnotime$lnL
  TL[[i]] <- max(nodeHeights(fns$tdList[[1]]$phy))
  fits[[i]] <- profiles(1, fns$ARD.R2, fns$tdList, rbind(true.pars[[i]]), res=NULL, seq=seq(0.1, TL[[i]]-0.1, 0.1), cores=8)[[1]]
}

treeid <- 392
true.pars <- c(0.1523248, 0.4984838)
fits <- list()
true.lnL <- list()
#true.pars <- list()
TL <- list()
#treeid <- list()
switched <- list(shuffled=list(), replaced=list())
for(i in 1:nsims){
  #treeid[[i]] <- sample(1:length(trees), 1)
  #true.pars[[i]] <- runif(2, 0.1, 1)
  trait <- diversitree::sim.character(trees[[treeid]], pars=true.pars, x0=sample(c(0,1), 1), model = "mk2")
  strait <- trait
  switched$shuffled[[i]] <- sample(1:length(trait), 10, replace=FALSE)
  switched$replaced[[i]] <- sample(switched$shuffled[[i]], 10, replace=FALSE)
  strait[switched$shuffled[[i]]] <- trait[switched$shuffled[[i]]][switched$replaced[[i]]]
  fns <- make.bisse.fns(trees[[treeid]], data.frame("trait"=strait))
  fits.ARDnotime <- fitFns(5, fns$ARD.notime, fns$tdList, res=NULL)
  true.lnL[[i]] <- fits.ARDnotime$lnL
  TL[[i]] <- max(nodeHeights(fns$tdList[[1]]$phy))
  fits[[i]] <- profiles(1, fns$ARD.R2, fns$tdList, starts=rbind(true.pars), res=NULL, seq=seq(0.1, TL[[i]]-0.1, 0.1), cores=8, start4=FALSE)[[1]]
}

results <- list(fits=fits, TL=TL, lnL=true.lnL, pars=true.pars, treeid=treeid, switched=switched)
res <- lapply(fits, function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
res <- do.call(cbind, res)
res <- apply(res, 2, function(x) x-min(x))
saveRDS(results, "../output/simulationstudy2.rds")

pdf("../output/simulation_studies.pdf")
  plot(c(0, 3.8), c(0, 2), type="n")
  lapply(1:ncol(res1),function(x) lines(seq(0.1, 3.6, 0.1), res1[,x], col=x+1))
  res <- lapply(fits, function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
  res <- do.call(cbind, res)
  res <- apply(res, 2, function(x) x-min(x))
  plot(c(3.8, 0), c(0, 2), type="n", xlim=c(3.8,0))
  lapply(1:ncol(res),function(x) lines(seq(0.1, TL[[x]]-0.1, length.out=length(res[,x])), res[,x], col=x+1))
  abline(v=0.6)
dev.off()

