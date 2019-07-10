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

lbaFiles <- list.files("../data/Posterior Probability Sampling of Trees/Trees Exhibiting Long Branch Attraction/Posterior probability distribution of time trees/")
lbaTrees <- lapply(lbaFiles,function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Exhibiting Long Branch Attraction/Posterior probability distribution of time trees/", x, sep="")))

goodFiles <- list.files("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/")
goodTrees <- lapply(goodFiles[grep("good", goodFiles)], function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/", x, sep="")))
cyanosTrees <- lapply(goodFiles[grep("Cyanos", goodFiles)], function(x) read.tree(paste("../data/Posterior Probability Sampling of Trees/Trees Without Known Long Branch Attraction Artifacts/Posterior distribution of trees where branch lengths are scaled to age/", x, sep="")))

lba.ntips <- lapply(lbaTrees, function(x) unname(sapply(x, function(y) length(y$tip.label))))
good.ntips <- lapply(goodTrees, function(x) unname(sapply(x, function(y) length(y$tip.label))))
cyanos.ntips <- lapply(cyanosTrees, function(x) unname(sapply(x, function(y) length(y$tip.label))))

lba.species <- unique(unlist(sapply(lbaTrees, function(x) x[[1]]$tip.label)))
good.species <- unique(unlist(sapply(goodTrees, function(x) x[[1]]$tip.label)))
cyanos.species <- unique(unlist(sapply(cyanosTrees, function(x) x[[1]]$tip.label)))
all.species <- unique(c(lba.species, good.species, cyanos.species))

dat <- clean.data()
tree <- lbaTrees[[2]][[1]]
trees <- list(goodTrees[[1]][1:2], goodTrees[[2]][1:2], goodTrees[[3]][1:2], goodTrees[[4]][1:2])

make.treedat <- function(tree, dat){
  ## set branch lengths to by
  #tree$edge.length <- tree$edge.length/1000
  
  ## Combine tree 
  tdcy <- make.treedata(tree, dat, name_column=0)
  rownames(tdcy$dat) <- tdcy$phy$tip.label
  colnames(tdcy$dat) <- gsub("/", "_", colnames(tdcy$dat), fixed=TRUE)
  nc <- ncol(tdcy$dat)
  tdcyList <- lapply(1:nc, function(x) select(tdcy, x))
  tdcyList <- lapply(1:nc, function(x) filter_(tdcyList[[x]], paste("!is.na(",names(tdcyList[[x]]$dat),")", sep="")))
  treedat<-tdcy
  treedat <- list(td=tdcy,tdList=tdcyList)
  return(treedat)
  
}

tree.dat <- make.treedat(tree, dat)

plot(tree.dat$td$phy)

## @knitr simplest salinity as continuous

dat_simp <- dat %>% 
  mutate(sal_simp = case_when(Freshwater_habitat == 1  ~ 0, 
                        Freshwater_habitat == 0 ~ 35,
                        TRUE ~ 0))
row.names(dat_simp) <- row.names(dat)

tree.dat <- make.treedat(tree, dat_simp)

obj<-contMap(tree.dat$td$phy,tree.dat$td[["sal_simp"]],plot=FALSE)
plot(obj,type="fan",legend=0.7*max(nodeHeights(tree)))

## @knitr drawing from distributions
require(truncdist)

dat_dist <- dat # replicate the dataset
dat_dist$sal_dist <- c(rep(NA,length(dat_dist$Freshwater_habitat))) # make an empty salinity distribution column
for (i in 1:length(dat_dist$Freshwater_habitat)){
  if (dat_dist$Freshwater_habitat[i] == 1){
    dat_dist$sal_dist[i] <- rtrunc(1,a=0,b=Inf,spec="norm",mean=0,sd=5) # truncated normal distribution so no negative salinities
  }
  if (dat_dist$Freshwater_habitat[i] == 0){
    dat_dist$sal_dist[i] <- rtrunc(1,a=-Inf,b=Inf,spec="norm",mean=35,sd=4)
  }
} # loops over the dataframe and draws from one of two distributions to assign salinity optima based on habitat

tree.dat <- make.treedat(tree,dat_dist)

obj2<-contMap(tree.dat$td$phy,tree.dat$td[["sal_dist"]],plot=FALSE)
plot(obj2,legend=0.7*max(nodeHeights(tree)))


# cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))
# tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
# add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
#                   y=-max(nodeHeights(tree)),fsize=0.8)
# 
