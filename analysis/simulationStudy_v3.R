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
    fns <- make.bisse.fns(trees[[treeids[[i]][2]]], data.frame("trait"=trait))
    fits.ARDnotime <- fitFns(10, fns$ARD.notime, fns$tdList, res=NULL)
    true.lnL[[i]] <- fits.ARDnotime$lnL
    #fits[[i]] <- profiles(1, fns$ARD.R2, fns$tdList, starts=cbind(1, 1, Qs[[i]][[1,2]], Qs[[i]][[2,1]]), res=NULL, seq=seq(0.1, TLs[[i]]-0.1, length.out=37), cores=8)[[1]]
    startx <- do.call(rbind, lapply(1:37, function(x) cbind(1, 1, Qs[[i]][1], Qs[[i]][2])))
    fits[[i]] <- profiles(1, fns$ARD.R2, fns$tdList, starts=list(startx), res=NULL, seq=seq(0.1, TLs[[i]]-0.1, length.out=37), cores=8, start4=TRUE) 
}
lnLs <- unlist(sapply(lnLs,function(x) as.numeric(as.character(x))))
profs <- lapply(fits, function(x) do.call(rbind, lapply(x[[1]], function(z) z[which(z$value==max(z$value)), c(5, 1:4)])))
profs2 <- lapply(1:length(profs), function(x) profs[[x]][,1]-lnLs[x])

plot(c(0, 4), c(0, 6), type="n", xlim=c(3.8, 0))
lapply(1:length(profs2), function(x) lines(seq(0.1, TLs[notnas][[x]]-0.1, length.out=37), profs2[[x]], col="gray20"))
cs <- apply(do.call(rbind, profs2),2, mean)
plot(seq(0, 3.5, length.out=37), cs, xlim=c(3.8, 0))
abline(v=c(0.6, 1.7))

maxInd <- lapply(profs2, function(x) which(x==max(x)))
maxT <- sapply(1:length(maxInd), function(x) seq(0.1, TLs[notnas][[x]]-0.1, length.out=37)[maxInd[[x]]])
maxD <- sapply(1:length(maxInd), function(x) profs2[[x]][maxInd[[x]]])
o <- order(maxD, decreasing=TRUE)

par(ask=TRUE, mfrow=c(3,3))
for(i in 1:length(maxInd)){
  plot(seq(0.1, TLs[notnas][o][[i]]-0.1, length.out=37), profs2[o][[i]], xlim=c(3.8, 0))
  abline(v=maxT[o][i], lty=2)
  abline(v=c(0.6, 1.7), lty=2, col="red")
}

pdf("../output/simulationstudytrees.pdf", height=20, width=10)
i=10
par(mfrow=c(1,2))
plot(trees[[simtreeids[i,1]]], cex=0.5)
abline(v=max(branching.times(trees[[simtreeids[i,1]]]))-c(0.6, 1.7), lty=2, col="red")
plot(trees[[simtreeids[i,2]]], cex=0.5)
abline(v=max(branching.times(trees[[simtreeids[i,2]]]))-c(0.6, 1.7), lty=2, col="red")
par(mfrow=c(1,1))
plot(c(0, 3.8), c(0, 5), type="n")
tmp1 <- ltt(trees[[simtreeids[i,1]]], plot=FALSE)
tmp2 <- ltt(trees[[simtreeids[i,2]]], plot=FALSE)
lines(tmp2$times, log(tmp2$ltt), col="red")
lines(tmp1$times, log(tmp1$ltt), col="black")
abline(v=c(3.8-1.7, 3.8-0.6))
dev.off()

profs <- lapply(fits[[1]], function(x) lapply(x[[1]], function(y) do.call(rbind, lapply(y, function(z) z[which(z$value==max(z$value)), c(5, 1:4)]))))

results <- list(fits=fits, lnL=true.lnL, treeids=treeids, true.pars=Qs, traits=traits)
saveRDS(results, file="../output/simulation_study_r2.rds")
results2 <- readRDS("../output/simulation_study_r2.rds")
results1 <- readRDS("../output/simulation_study1.rds")
fits <- results$fits
lnL <- results$true.lnL
treeids <- results$treeids
true.pars <- results$treeids

getfits <- function(fits){
  notnas <- which(is.na(fits)==FALSE)
  res <- lapply(fits[notnas], function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
  res <- do.call(cbind, res)
  res <- apply(res, 2, function(x) x-min(x))
  res1 <- res
  res1
}

r1 <- getfits(results1$fits)
r2 <- getfits(lapply(results2$fits, function(x) x[[1]]))
plot(apply(r2, 1, sum))
plot(apply(r1, 1, sum))
plot(seq(0, 3.7, length.out=36), apply(r1, 1, sum)+apply(r2, 1, sum)[1:36], )
for(i in 1:100){
  plot(seq(0.1, 3.7, length.out=37), r2[,i],xlim=c(3.7,0))
  lines(seq(0.1, 3.7, length.out=37), r2[,i])
  abline(v=c(0.6, 1.8))
}
maxes <- apply(r2, 2, function(x) which(x==max(x)))
par(mfrow=c(3,3), ask=TRUE)
falsepos <- which(maxes %in% 16:20)
for(i in which(maxes %in% 16:20)){
  plot(seq(0.1, 3.7, length.out=37), r2[,i],xlim=c(3.7,0))
  lines(seq(0.1, 3.7, length.out=37), r2[,i])
  abline(v=c(0.6, 1.8))
}

treeids <- results2$treeids[!sapply(results2$fits, is.na)]
sum.falsepos <- data.frame(t(sapply(treeids[falsepos], function(x) treesources[x])), maxlnL=apply(r2[,falsepos], 2, max))
sum.falsepos


tmp <- trees[treeids[[falsepos[2]]]]
par(mfrow=c(1,2))
plot(tmp[[1]],show.tip.label =FALSE)
plot(tmp[[2]], show.tip.label =FALSE)
plot(tmp[[1]]$edge.length, tmp[[2]]$edge.length)
abline(0,1)

## Identical trees (essentially gives two peaks? WTH?)
fits <- list()
true.lnL <- list()
true.pars <- list()
treeids <- list()
traits <- list()
nsims <- 10
for(i in 1:nsims){
  treeids[[i]] <- c(910, 857)
  trait <- 1
  while(length(unique(trait))<=1){
    trait <- diversitree::sim.character(trees[[treeids[[i]][1]]], pars=c(Qs[[i]][1], Qs[[i]][2]), x0=sample(c(0,1), 1), model = "mk2")
  }
  traits[[i]] <- trait
  fns <- make.bisse.fns(trees[[treeids[[i]][2]]], data.frame("trait"=trait))
  fits.ARDnotime <- fitFns(10, fns$ARD.notime, fns$tdList, res=NULL)
  true.lnL[[i]] <- fits.ARDnotime$lnL
  #fits[[i]] <- profiles(1, fns$ARD.R2, fns$tdList, starts=cbind(1, 1, Qs[[i]][[1,2]], Qs[[i]][[2,1]]), res=NULL, seq=seq(0.1, TLs[[i]]-0.1, length.out=37), cores=8)[[1]]
  startx <- do.call(rbind, lapply(1:37, function(x) cbind(1, 1, Qs[[i]][1], Qs[[i]][2])))
  fits[[i]] <- profiles(1, fns$ARD.R2, fns$tdList, starts=list(startx), res=NULL, seq=seq(0.1, TLs[[i]]-0.1, length.out=37), cores=8, start4=TRUE) 
}
fits <- lapply(fits, function(x) x[[1]])

r3 <- getfits(fits)
r4 <- getfits(fits)
r5 <- do.call(cbind, lapply(1:ncol(r4), function(x) ifelse(r4[,x]<as.numeric(as.character(true.lnL[[x]])), as.numeric(as.character(true.lnL[[x]])), as.numeric(as.character(r4[,x])))))

model.fits <- lapply(fits, function(x) do.call(rbind, lapply(x, function(y) y[which(y$value==max(y$value))[1],])))
plot(c(0, 3.7), c(0, max(r3)+0.5), type="n", xlim=c(3.7, 0))
lapply(1:10, function(x) lines(seq(0, 3.7, length.out=37), r3[,x]))
lapply(1:10, function(x) lines(seq(0, 3.7, length.out=37), r4[,x]))
abline(v=c(1.8, 0.6))

plot(seq(0, 3.7, length.out=37), apply(r3, 1, mean), xlim=c(3.7, 0))
lines(seq(0, 3.7, length.out=37), apply(r3, 1, mean))
plot(seq(0, 3.7, length.out=37), apply(r4, 1, mean), xlim=c(3.7, 0))
lines(seq(0, 3.7, length.out=37), apply(r4, 1, mean))
abline(v=c(1.8, 0.6))

par(mfrow=c(3,3), ask=TRUE)
for(i in 1:10){
  plot(seq(0.1, 2.7, length.out=37), r3[,i], type="n", xlim=c(2.7,0), ylim=c(0,2))
  lines(seq(0.1, 2.7, length.out=37), r3[,i])
  lines(seq(0.1, 2.7, length.out=37), r4[,i], lty=2)
  lines(seq(0.1, 2.7, length.out=37), r5[,i], lty=4, col="red")
  abline(v=c(0.6, 1.8))
}
Qs[1]
par(mfrow=c(1,2), ask=FALSE)
plot(trees[[treeids[[1]][1]]], show.tip.label=FALSE)
tiplabels(pch=21, bg=traits[[1]])
plot(trees[[treeids[[1]][2]]], show.tip.label=FALSE)
tiplabels(pch=21, bg=traits[[1]])


for(i in 1:10){
  par(mfrow=c(1,2), ask=TRUE)
  parsi <- c(unlist(model.fits[[i]][which(model.fits[[i]]$value==max(model.fits[[i]]$value))[1],1:4]), seq(0.1, max(branching.times(trees[[treeids[[i]][1]]])), length.out=37)[which(model.fits[[i]]$value==max(model.fits[[i]]$value))[1]])
  plot.asrR2(parsi, trees[[treeids[[i]][1]]], traits[[i]], show.tip.label=FALSE)
  tiplabels(pch=21, bg=c("gray10", "gray90")[traits[[i]]+1])
  #plot(trees[[treeids[[i]][1]]], type="fan", show.tip.label=FALSE)
  #tiplabels(pch=21, bg=c("gray10", "gray90")[traits[[i]]+1])
  #nodelabels(pch=21,bg=c("gray10", "gray90")[attr(traits[[i]], "node.state")+1] )
  
  
  parsi <- c(unlist(model.fits[[i]][which(model.fits[[i]]$value==max(model.fits[[i]]$value))[1],1:4]), seq(0.1, max(branching.times(trees[[treeids[[i]][2]]])), length.out=37)[which(model.fits[[i]]$value==max(model.fits[[i]]$value))[1]])
  plot.asrR2(parsi, trees[[treeids[[i]][2]]], traits[[i]], show.tip.label=FALSE)
  tiplabels(pch=21, bg=c("gray10", "gray90")[traits[[i]]+1])
  #plot(trees[[treeids[[i]][2]]], type="fan", show.tip.label=FALSE)
  #tiplabels(pch=21, bg=c("gray10", "gray90")[traits[[i]]+1])
  #nodelabels(pch=21,bg=c("gray10", "gray90")[attr(traits[[i]], "node.state")+1] )
}

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
results <- readRDS("../output/simulationstudy2.rds")
fits <- results$fits
TL <- results$TL
res <- lapply(fits, function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
res <- do.call(cbind, res)
res <- apply(res, 2, function(x) x-x[length(x)])
res[res < 0] <- 0

pdf("../output/simulation_study2.pdf")
  #plot(c(0, 3.8), c(0, 2), type="n")
  #lapply(1:ncol(res1),function(x) lines(seq(0.1, 3.6, 0.1), res1[,x], col=x+1))
  #res <- lapply(fits, function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
  #res <- do.call(cbind, res)
  #res <- apply(res, 2, function(x) x-min(x))
  plot(c(3.8, 0), c(0, 1.5), type="n", xlim=c(3.8,0), xlab="Time", ylab="")
  lapply(1:ncol(res),function(x) lines(seq(0.1, TL[[x]]-0.1, length.out=length(res[,x])), res[,x], col="gray70"))
  lines(seq(0.1, mean(unlist(TL))-0.1, length.out=nrow(res)), apply(res, 1, mean), lwd=3, col="#e74c3c")
  abline(v=0.6)
dev.off()

