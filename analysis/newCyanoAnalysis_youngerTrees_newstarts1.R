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
trees <- c(youngTrees[[1]][seq(1,100, length.out=10)], youngTrees[[2]][seq(1,100, length.out=10)])
trees <- lapply(trees, function(x){x$edge.length <- x$edge.length/1000; x})
TLs <- sapply(1:length(trees), function(x) max(branching.times(trees[[x]])))

Fns <- lapply(trees, function(x) make.bisse.fns(x, dat))
#registerDoParallel(cores=5)
#fits.ARDnotime <- lapply(Fns, function(fns) fitFns(5, fns$ARD.notime, fns$tdList, res=NULL))

load("./new.researchFull.rds")
oldresults <- new.researchFull[names(new.researchFull)]
oldres <- lapply(oldresults, function(trait) do.call(rbind, lapply(1:length(trait), function(time) c(trait[[time]]$lnLik, trait[[time]]$par))))
oldstarts <- lapply(oldres, function(x) x[,-1])
oldnames <- names(oldstarts)
newnames <- sapply(Fns[[1]]$tdList, function(x) colnames(x$dat)[1])
newnames[newnames=="temp"] <- "Nonfreshwater_habitat"
oldnames <- gsub("/", "_", oldnames, fixed=TRUE)
starts <- oldstarts[match(newnames, oldnames)]
starts[[which(sapply(starts, length)==0)]] <- starts[[3]]
starts <- lapply(starts, function(x) x[nrow(x):1,])
starts1 <- readRDS("../output/newstarts1.rds")
starts1 <- lapply(starts1, as.data.frame)
names(starts1) <- names(starts)


#seq <- seq(0.2, 3.6, 0.5); Fns <- fns; fns <- Fns$ARD.R2; starts<-fits.ARDnotime; res=NULL; tds <- fns$tdList
seqs <- lapply(1:length(TLs), function(x) seq(0.1, TLs[x], length.out=37))
profs <- list()
for(i in 1:length(Fns)){
  fns <- Fns[[i]]
  profs[[i]] <- profiles(1, fns$ARD.R2, fns$tdList, starts1, res=NULL, seq=seq(0.1, TLs[[i]], length.out=37), cores=10, start4=TRUE) 
}

#fns <- Fns[[1]]
#prof.tree1 <- profiles(1, fns$ARD.R2, fns$tdList, starts, res=NULL, seq=seq(0.1, 3.7, length.out=37), cores=8, start4=TRUE)
#fns <- Fns[[2]]
#prof.tree2 <- profiles(1, fns$ARD.R2, fns$tdList, starts, res=NULL, seq=seq(0.1, 3.7, length.out=37), cores=8, start4=TRUE)
#saveRDS(prof.tree1, "../output/proftree1.rds")
#saveRDS(prof.tree2, "../output/proftree2.rds")
#prof.tree1 <- readRDS("../output/proftree1.rds")
#prof.tree2 <- readRDS("../output/proftree2.rds")
saveRDS(profs, "../output/profs_20youngTrees_starts1.rds")
profs <- readRDS("../output/profs_20youngTrees_starts1.rds")

profiles <- lapply(profs, function(x) lapply(x, function(y) do.call(rbind, lapply(y, function(z) z[which(z$value==max(z$value)), c(5, 1:4)]))))
sumprof <- lapply(profiles, function(x) lapply(lapply(x, function(y) cbind(y$value)), summarizeProfile))
cumsums <- lapply(sumprof, function(x) sapply(x, function(y) y$slnL[1:37]))
cumsums <- lapply(cumsums, function(x){x[x<0] <- 0; x})
cumsums <- lapply(cumsums, function(x) apply(x, 1, sum))
allcs <- do.call(rbind, cumsums)
allcs <- apply(allcs, 2, median)
allcstime <- apply(do.call(rbind, seqs), 2, median)

plot(c(0,2.7), c(0,33), type="n", xlim=c(2.7, 0))
lapply(1:length(cumsums), function(x) lines(seqs[[x]], cumsums[[x]], col="gray60"))
lines(allcstime, allcs, col="black", lwd=5, lty=2)
abline(v=c(0.53, 1.68))

#tmp <- lapply(1:length(profiles), function(tr) lapply(1:length(profiles[[1]]), function(trait) data.frame(profiles[[tr]][[trait]][1:37,],tree=tr, trait=trait, seq=1:37)))
#tmp <- lapply(1:length(tmp), function(x) do.call(rbind, tmp[[x]]))
#tmp <- do.call(rbind, tmp)
#tmp <- group_by(tmp, trait, seq)
#newstarts <- summarize(tmp, max(value), r1[value==max(value)], r2[value==max(value)], q01.y0[value==max(value)], q10.y0[value==max(value)])
#colnames(newstarts) <- c("trait", "seq", "value", "r1", "r2", "q01", "q10")
#starts.new <- lapply(1:26, function(x) newstarts[newstarts$trait==x,4:7])
#starts1 <- lapply(1:26, function(x) rbind(starts.new[[x]][2:37,], starts.new[[x]][37,]))
#starts2 <- lapply(1:26, function(x) rbind(starts.new[[x]][1,], starts.new[[x]][1:36,]))
#saveRDS(starts1, "../output/newstarts1.rds")
#saveRDS(starts2, "../output/newstarts2.rds")


#tmp1 <- sapply(prof.tree1, function(x) sapply(x, function(y) y[which(y[,5]==max(y[,5])),5][1]))
#tmp2 <- sapply(prof.tree2, function(x) sapply(x, function(y) y[which(y[,5]==max(y[,5])),5][1]))
#tmp1 <- apply(tmp1, 2, function(x) x-x[37])
#tmp2 <- apply(tmp2, 2, function(x) x-x[37])
#tmp1[tmp1 < 0] <- 0
#tmp2[tmp2 < 0] <- 0
#cs1 <- apply(tmp1, 1, cumsum)
#cs2 <- apply(tmp2, 1, cumsum)
#plot(seq,cs1[26,], xlim=c(2.7,0), type="n", ylim=c(0,35))
#lines(seq, cs1[26,], col="red")
#lines(seq, cs2[26,], col="blue")
#abline(v=c(0.7,1.7))

#habitat1 <- prof.tree1[which(newnames %in% c("Freshwater_habitat", "Nonfreshwater_habitat"))]
#habitat2 <- prof.tree2[which(newnames %in% c("Freshwater_habitat", "Nonfreshwater_habitat"))]
#habitat2 <- lapply(habitat2, function(x) sapply(x, function(y) y[which(y[,5]==max(y[,5])),5][1]))
#habitat2 <- lapply(habitat2, function(x) x-x[length(x)])
#habitat2 <- lapply(habitat2, function(x){x[x<0] <- 0; x})
#habitat1 <- lapply(habitat1, function(x) sapply(x, function(y) y[which(y[,5]==max(y[,5])),5][1]))
#habitat1 <- lapply(habitat1, function(x) x-x[length(x)])
#habitat1 <- lapply(habitat1, function(x){x[x<0] <- 0; x})

#pdf("../output/habitatShifts.pdf")
#par(mfrow=c(1,2))
#plot(seq(0.1, 3.7, length.out=37), habitat1[[1]],type="n", ylim=c(0,5), xlim=c(2.7,0), main="Freshwater")
#lines(seq(0.1, 3.7, length.out=37), habitat1[[1]], lty=2, col="blue")
#lines(seq(0.1, 3.7, length.out=37), habitat2[[1]], lty=1, col="blue")
#lines(seq(0.1, 3.7, length.out=37), rev(os))
#plot(seq(0.1, 3.7, length.out=37), habitat1[[1]],type="n", ylim=c(0,5), xlim=c(2.7,0), main="Nonfreshwater (temp)")
#lines(seq(0.1, 3.7, length.out=37), habitat2[[2]], lty=1, col="red")
#lines(seq(0.1, 3.7, length.out=37), habitat1[[2]], lty=2, col="red")
#lines(seq(0.1, 3.7, length.out=37), rev(os))
#dev.off()

#tmp <- lapply(prof.tree1, function(x) t(sapply(x, function(y) y[which(y[,5]==max(y[,5])),1:5])))
#plot(rev(seq), unlist(tmp[[21]][,'value']))
#lines(seq, tmp[[21]][,'value'], xlim=c(0,3.8))


#saveRDS(prof.ARDR2, "../output/prof_lba2.1_r1.rds")
#sprof.ARDR2 <- smooth.profiles(fns$ARD.R2,fns$tdList, prof.ARDR2, seq=seq(0.1, 3.6, length.out=32, cores=8))
#tryerr <- (sapply(sapply(prof.ARDR2, function(x) class(x[[1]])=="try-error"), function(x) x[1]))
#lnL.ARDR2 <- lapply(prof.ARDR2[-which(tryerr)], function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
#lnL.ARDR2 <- do.call(rbind, lnL.ARDR2)
#cslnL.ARDR2 <- apply(lnL.ARDR2, 2, cumsum)
#plot(seq(3.6, 0.1, length.out=32), rev(cslnL.ARDR2[nrow(cslnL.ARDR2),]))

#pdf("../output/goodTreesPosterior.pdf")
#plot(seq(3.6, 0.1, length.out=32), rev(cslnL.ARDR2[nrow(cslnL.ARDR2),]), type="n", ylim=c(0, 40))
#profs <- list()
#load("./new.researchFull.rds")
#oldresults <- new.researchFull[names(dat)]
#oldres <- lapply(oldresults, function(trait) do.call(rbind, lapply(1:length(trait), function(time) c(trait[[time]]$lnLik, trait[[time]]$par))))
#oldcumsum <- apply(sapply(oldres, function(x) x[,1]), 1, sum)
#oldcumsum <- oldcumsum-min(oldcumsum)

#for(i in 1:8){
#  profs[[i]] <- readRDS(paste("../output/prof_good_",i,".rds", sep=""))
#  #tryerr <- (sapply(sapply(prof.ARDR2, function(x) any(sapply(x, class)=="try-error")), function(x) x[1]))
#  lnL.ARDR2 <- lapply(profs[[i]], function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
#  lnL.ARDR2 <- do.call(rbind, lnL.ARDR2)
#  cslnL.ARDR2 <- apply(lnL.ARDR2, 2, cumsum)
#  lines(seq(3.6, 0.1, length.out=32), rev(cslnL.ARDR2[nrow(cslnL.ARDR2),]-min(cslnL.ARDR2)))
#  lines(rev(seq(0.1, 3.7, length.out=37)), oldcumsum, xlim=c(3.8, 0), col="red")
#}

#tmp <- lapply(profs, function(reps) lapply(reps, function(traits) do.call(rbind, lapply(traits, function(times) times[which(times$value==max(times$value)),]))))
#tmp <- lapply(1:25, function(x) lapply(tmp, function(y) y[[x]]))
#allprofs <- lapply(1:32, function(x) lapply(tmp, function(y) do.call(rbind, lapply(y, function(z) z[x,]))))
#allprofs <- lapply(1:25, function(x) lapply(allprofs, function(y) y[[x]]))
#allprofs <- lapply(allprofs, function(x) do.call(rbind, lapply(1:32, function(y) data.frame(time=seq(0.1,3.6, length.out=32)[y], x[[y]]))))

#par(mfrow=c(2,2), ask=TRUE)
#for(i in 1:25){
#  for(j in 2:5){
#    plot(c(allprofs[[i]][,1],rev(seq(0.1,3.7, 0.1))), c(log(allprofs[[i]][,j]+0.0001),log(oldres[[i]][,j]+0.0001)),type="n", pch=21, bg="gray20", cex=0.5)
#    points(allprofs[[i]][,1], log(allprofs[[i]][,j]+0.0001), pch=21, cex=(allprofs[[i]][,'value']-min(allprofs[[i]][,'value'])+0.001)/10, bg="gray20")
#    lines(rev(seq(0.1,3.7, 0.1)), log(oldres[[i]][,j]+0.0001))
#  }
#}

#for(i in c("1_r1", "1_r2", "2.1_r1")){
#  prof.ARDR2 <- readRDS(paste("../output/prof_lba",i,".rds", sep=""))
#  tryerr <- (sapply(sapply(prof.ARDR2, function(x) any(sapply(x, class)=="try-error")), function(x) x[1]))
#  if(sum(tryerr)>0){
#    lnL.ARDR2 <- lapply(prof.ARDR2[-which(tryerr)], function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
#  } else {
#    lnL.ARDR2 <- lapply(prof.ARDR2, function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
#  }
#  lnL.ARDR2 <- do.call(rbind, lnL.ARDR2)
#  cslnL.ARDR2 <- apply(lnL.ARDR2, 2, cumsum)
#  lines(seq(3.6, 0.1, length.out=32), rev(cslnL.ARDR2[nrow(cslnL.ARDR2),]-min(cslnL.ARDR2)), col="red")
#}
#dev.off()

## Figures
#sumRes <- lapply(oldres, summarizeProfile)
#o <- order(sapply(sumRes, function(x) x$maxT))
#sumRes <- sumRes[o]

### Shift Support Figure
#pdf("../output/shiftsupportFigure1.pdf")
#shiftFigure <- function(sumRes){
#  par(mar=c(3,5,12,0.5))
#  plot(c(1, length(sumRes)+2), c(3.8, 0), ylim=c(3.8, 0), type="n", ylab="Time (billion years before present)",  xaxt="n")
#  axis(3, at=c(1:length(names(sumRes)),27), labels=gsub("_", " ", c(names(sumRes), "Cumulative")), las=3, cex=0.5)
#  polygon(c(0,0,28,28), c(0.7, 0.5, 0.5 ,0.7), col=makeTransparent("gray80", 100), border="white")
#  polygon(c(0,0,28,28), c(2.2, 2.4, 2.4 ,2.2), col=makeTransparent("gray80", 100), border="white")
#  lapply(1:length(sumRes), function(x) polygon(c(x, rep(x,37)-rev(sumRes[[x]]$slnL)/12,x, rep(x,37)+(sumRes[[x]]$slnL)/12), c(seq(0, 3.8,0.1), seq(3.7,0.1,-0.1)), border = "gray80",col=dlnLColors(rev(sumRes[[x]]$dlnL), "red", max=7)))
#  #lapply(1:length(sumRes), function(x) points(rep(x,37), seq(0.1,3.7,0.1), pch="|", col=dlnLColors(rev(sumRes[[x]]$slnL), "red")))
#  lapply(1:length(sumRes), function(x) points(x, 3.8-sumRes[[x]]$maxT, pch=21, col=dlnLColors(sumRes[[x]]$dlnL,"black", max=7), bg=dlnLColors(sumRes[[x]]$dlnL,"red", max=7)))
#  oldcumsum <- apply(sapply(sumRes, function(x) x$slnL),1,sum)
#  cs <-rev(oldcumsum)-20
#  cs[cs<0] <- 0
#  polygon(c(27, rep(27,37)-cs/40,27, rep(27,37)+rev(cs)/40), c(seq(0, 3.8,0.1), seq(3.7,0.1,-0.1)), border = "gray80",col=dlnLColors(rev(max(oldcumsum)), "red", max=max(oldcumsum)))
#  points(rep(27,37), seq(0.1,3.7,0.1),cex=1, pch="|", col=dlnLColors(cs,"red", max=max(cs)))
#  r <- rle(cs)
#  peaks <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, times = r$lengths))
#  maxpeak <- which(cs[peaks]==max(cs[peaks]))
#  peak2 <- which(cs[peaks[-maxpeak]]==max(cs[peaks[-maxpeak]]))
#  peaks <- peaks[c(peak2, maxpeak)]
#  points(rep(27, length(peaks)), seq(0.1,3.7,0.1)[peaks], cex=1, pch=21, col="black", bg="red")
#  abline(h=seq(0.1,3.7,0.1)[peaks], lty=2)
#  #points(rep(27,2), seq(0.1,3.7,0.1)[which(cs>18)],cex=1, pch=21, col="black",bg="red")
#}
#shiftFigure(sumRes)
#dev.off()

### Write a table of results
#o <- order(sapply(sumRes, function(x) x$dlnL))
#sumRes <- sumRes[rev(o)]
#timeslice <- sapply(sumRes, function(x) 3.8-x$maxT)
#shiftcluster <- rep("", length(timeslice))
#shiftcluster[abs(timeslice-0.6) <=0.2] <- "NeoProterozoic"
#shiftcluster[abs(timeslice-1.8) <=0.2] <- "PaleoProterozoic"
#sumTable <- data.frame("Trait"=names(sumRes), "Max Support"=sapply(sumRes, function(x) x$dlnL), "Time Estimate"=timeslice, "Shift Cluster"=shiftcluster)
#o <- order(sumTable$Shift.Cluster, sumTable[,2], decreasing = TRUE)
#write.csv(sumTable[o,], "../output/ResultsTable.csv",row.names=FALSE)
