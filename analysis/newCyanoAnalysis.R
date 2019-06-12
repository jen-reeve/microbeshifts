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
require(treeplyr)
setwd("~/GitHub/microbeshifts/analysis")
file_path<-"shiftfunctions.R"
source(file_path)
file_path<-"cyanofunctions.R"
source(file_path)

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
tree <- lbaTrees[[2]][[1]]
trees <- list(goodTrees[[1]][1:2], goodTrees[[2]][1:2], goodTrees[[3]][1:2], goodTrees[[4]][1:2])

fns <- make.bisse.fns(tree, dat)
registerDoParallel(cores=8)
fits.ARDnotime <- fitFns(5, fns$ARD.notime, fns$tdList, res=NULL)


#seq <- seq(0.2, 3.6, 0.5); Fns <- fns; fns <- Fns$ARD.R2; starts<-fits.ARDnotime; res=NULL; tds <- fns$tdList
prof.ARDR2 <- profiles(1, fns$ARD.R2, fns$tdList, fits.ARDnotime, res=NULL, seq=seq(0.1, 3.6, length.out=32), cores=8)


saveRDS(prof.ARDR2, "../output/prof_lba2.1_r1.rds")
sprof.ARDR2 <- smooth.profiles(fns$ARD.R2,fns$tdList, prof.ARDR2, seq=seq(0.1, 3.6, length.out=32, cores=8))
tryerr <- (sapply(sapply(prof.ARDR2, function(x) class(x[[1]])=="try-error"), function(x) x[1]))
lnL.ARDR2 <- lapply(prof.ARDR2[-which(tryerr)], function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
lnL.ARDR2 <- do.call(rbind, lnL.ARDR2)
cslnL.ARDR2 <- apply(lnL.ARDR2, 2, cumsum)
plot(seq(3.6, 0.1, length.out=32), rev(cslnL.ARDR2[nrow(cslnL.ARDR2),]))

pdf("../output/goodTreesPosterior.pdf")
plot(seq(3.6, 0.1, length.out=32), rev(cslnL.ARDR2[nrow(cslnL.ARDR2),]), type="n", ylim=c(0, 40))
profs <- list()
load("./new.researchFull.rds")
oldresults <- new.researchFull[names(dat)]
oldres <- lapply(oldresults, function(trait) do.call(rbind, lapply(1:length(trait), function(time) c(trait[[time]]$lnLik, trait[[time]]$par))))
oldcumsum <- apply(sapply(oldres, function(x) x[,1]), 1, sum)
oldcumsum <- oldcumsum-min(oldcumsum)

for(i in 1:8){
  profs[[i]] <- readRDS(paste("../output/prof_good_",i,".rds", sep=""))
  #tryerr <- (sapply(sapply(prof.ARDR2, function(x) any(sapply(x, class)=="try-error")), function(x) x[1]))
  lnL.ARDR2 <- lapply(profs[[i]], function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
  lnL.ARDR2 <- do.call(rbind, lnL.ARDR2)
  cslnL.ARDR2 <- apply(lnL.ARDR2, 2, cumsum)
  lines(seq(3.6, 0.1, length.out=32), rev(cslnL.ARDR2[nrow(cslnL.ARDR2),]-min(cslnL.ARDR2)))
  lines(rev(seq(0.1, 3.7, length.out=37)), oldcumsum, xlim=c(3.8, 0), col="red")
}

tmp <- lapply(profs, function(reps) lapply(reps, function(traits) do.call(rbind, lapply(traits, function(times) times[which(times$value==max(times$value)),]))))
tmp <- lapply(1:25, function(x) lapply(tmp, function(y) y[[x]]))
allprofs <- lapply(1:32, function(x) lapply(tmp, function(y) do.call(rbind, lapply(y, function(z) z[x,]))))
allprofs <- lapply(1:25, function(x) lapply(allprofs, function(y) y[[x]]))
allprofs <- lapply(allprofs, function(x) do.call(rbind, lapply(1:32, function(y) data.frame(time=seq(0.1,3.6, length.out=32)[y], x[[y]]))))



par(mfrow=c(2,2), ask=TRUE)
for(i in 1:25){
  for(j in 2:5){
    plot(c(allprofs[[i]][,1],rev(seq(0.1,3.7, 0.1))), c(log(allprofs[[i]][,j]+0.0001),log(oldres[[i]][,j]+0.0001)),type="n", pch=21, bg="gray20", cex=0.5)
    points(allprofs[[i]][,1], log(allprofs[[i]][,j]+0.0001), pch=21, cex=(allprofs[[i]][,'value']-min(allprofs[[i]][,'value'])+0.001)/10, bg="gray20")
    lines(rev(seq(0.1,3.7, 0.1)), log(oldres[[i]][,j]+0.0001))
  }
}



for(i in c("1_r1", "1_r2", "2.1_r1")){
  prof.ARDR2 <- readRDS(paste("../output/prof_lba",i,".rds", sep=""))
  tryerr <- (sapply(sapply(prof.ARDR2, function(x) any(sapply(x, class)=="try-error")), function(x) x[1]))
  if(sum(tryerr)>0){
    lnL.ARDR2 <- lapply(prof.ARDR2[-which(tryerr)], function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
  } else {
    lnL.ARDR2 <- lapply(prof.ARDR2, function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
  }
  lnL.ARDR2 <- do.call(rbind, lnL.ARDR2)
  cslnL.ARDR2 <- apply(lnL.ARDR2, 2, cumsum)
  lines(seq(3.6, 0.1, length.out=32), rev(cslnL.ARDR2[nrow(cslnL.ARDR2),]-min(cslnL.ARDR2)), col="red")
}
dev.off()

