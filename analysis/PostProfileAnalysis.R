## Putting together different runs across posterior of trees:
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


youngFiles <- list.files("../data/Posterior Probability Sampling of Trees/Younger Root MaxAge/")
youngTrees <- lapply(youngFiles,function(x) read.trees.from.mesquite(paste("../data/Posterior Probability Sampling of Trees/Younger Root MaxAge/", x, sep="")))

dat <- clean.data()
trees <- c(youngTrees[[1]][seq(1,100, length.out=10)], youngTrees[[2]][seq(1,100, length.out=10)])
trees <- lapply(trees, function(x){x$edge.length <- x$edge.length/1000; x})
TLs <- sapply(1:length(trees), function(x) max(branching.times(trees[[x]])))

Fns <- lapply(trees, function(x) make.bisse.fns(x, dat))
#registerDoParallel(cores=5)
#fits.ARDnotime <- lapply(Fns, function(fns) fitFns(5, fns$ARD.notime, fns$tdList, res=NULL))
seqs <- lapply(1:length(TLs), function(x) seq(0.1, TLs[x], length.out=37))

load("./new.researchFull.rds")
oldresults <- new.researchFull[names(new.researchFull)]
oldres <- lapply(oldresults, function(trait) do.call(rbind, lapply(1:length(trait), function(time) c(trait[[time]]$lnLik, trait[[time]]$par))))
oldstarts <- lapply(oldres, function(x) x[,-1])
oldnames <- names(oldstarts)
newnames <- sapply(Fns[[1]]$tdList, function(x) colnames(x$dat)[1])
newnames[newnames=="temp"] <- "Nonfreshwater_habitat"
oldnames <- gsub("/", "_", oldnames, fixed=TRUE)

prof1 <- readRDS("../output/profs_20youngTrees.rds")
prof2 <- readRDS("../output/profs_20youngTrees_starts1.rds")
prof3 <- readRDS("../output/profs_20youngTrees_starts2.rds")
prof4 <- readRDS ("../output/profs_20youngTrees_r2.rds")
profs <- list(prof1, prof2, prof3, prof4)

profiles <- lapply(profs, function(j) lapply(j, function(x) lapply(x, function(y) do.call(rbind, lapply(y, function(z) z[which(z$value==max(z$value)), c(5, 1:4)])))))
sumprof <- lapply(1:length(profiles[[1]]), function(x) list())
for(i in 1:length(profiles[[1]][[1]])){
  for(j in 1:length(profiles[[1]])){
    prfs <- lapply(profiles, function(x) x[[j]][[i]])
    byseq <- lapply(1:nrow(prfs[[1]]), function(x) do.call(rbind, lapply(1:length(prfs), function(z) setNames(prfs[[z]][x,], colnames(prfs[[1]])))))
    collapseback <- do.call(rbind, lapply(byseq, function(x) x[which(x[,1]==max(x[,1]))[1],]))
    sumprof[[j]][[i]] <- collapseback
  }
}

sumprof2 <- lapply(sumprof, function(x) lapply(lapply(x, function(y) cbind(y$value)), summarizeProfile))
cumsums <- lapply(sumprof2, function(x) sapply(x, function(y) y$slnL[1:37]))
cumsums <- lapply(cumsums, function(x){x[x<0] <- 0; x})
#tmp <- lapply(1:ncol(cumsums[[1]][[1]]),function(x) lapply(1:length(cumsums[[1]]), function(z)  sapply(1:length(cumsums), function(y) cumsums[[y]][[z]][,x])))
allRes <- lapply(1:ncol(cumsums[[1]]), function(x) do.call(cbind, lapply(1:length(cumsums), function(y) cumsums[[y]][,x])))
names(allRes) <- newnames
keep <- 1:26[-2]
medAllRes <- lapply(allRes, function(x) apply(x,1, median))
cumsumAll <- apply(do.call(rbind, medAllRes[keep]), 2, sum)
cumsums <- lapply(cumsums, function(x) apply(x, 1, sum))
plot(cumsumAll)

#par(ask=TRUE)
#for(i in 1:26){
#  plot(c(0, 3), y=c(-1, max(unlist(tmp[[i]]))+0.2), type="n", main=newnames[i])
#  lapply(1:20, function(x) lines(seqs[[x]], tmp[[i]][[x]], col=x+1))
#  abline(v=c(0.6, 1.8), lty=2, lwd=3, col="red")
#}

allcstime <- apply(do.call(rbind, seqs), 2, median)

#load("./new.researchFull.rds")
#oldresults <- new.researchFull[names(dat)]
#oldres <- lapply(oldresults, function(trait) do.call(rbind, lapply(1:length(trait), function(time) c(trait[[time]]$lnLik, trait[[time]]$par))))
#oldres <- lapply(lapply(oldres, function(x) x[37:1,]), summarizeProfile)
#oldres <- do.call(cbind, sapply(1:length(oldres), function(x) oldres[[x]]$slnL))
#oldcumsum <- apply(oldres, 1, sum)
#oldcumsum <- oldcumsum-min(oldcumsum)


pdf("../output/LikelihoodSurfaceFigure.pdf")
allcs <- cumsumAll
plot(allcstime, cumsumAll, ylim=c(0, 35), xlim=c(2.7,0), type="n", xlab="Billion years before present", ylab="Likelihood Support")
abline(v=c(0.6, 1.7), lty=2, lwd=3, col="red")
lapply(1:length(cumsums), function(x) lines(seqs[[x]], cumsums[[x]], col="gray80"))
lines(allcstime, cumsumAll, lty=1, lwd=3)
polygon(c(-20, 20, 20, -20), c(max(allcs), max(allcs), max(allcs)-4, max(allcs)-4), border=NA, col=makeTransparent("gray60", 50))
#lines(seq(0.1,3.7, 0.1), oldcumsum)
dev.off()


pdf("../output/shiftsupportFigure1.pdf")
sumRes <- lapply(allRes, function(x) rev(apply(x, 1, median)))[keep]
o <- order(sapply(sumRes, function(x) which(x==max(x))[1]))
sumRes <- sumRes[o]
shiftFigure <- function(sumRes, col="#e74c3c", width=12){
  maxTimes <- sapply(sumRes, function(x) rev(allcstime)[which(x==max(x))])
  dlnLs <- sapply(sumRes, function(x) max(x))
  transramp <- round(aRbor:::.scale(dlnLs, range = c(100,255), qq=c(0,0.25, 0.5, 0.75, 1), n = 100, log=FALSE),0)
  cols <- sapply(transramp, function(x) makeTransparent(col, x))
  par(mar=c(3,5,12,0.5))
  plot(c(1, 29), c(2.7, 0), ylim=c(2.7, 0), type="n", ylab="Time (billion years before present)",  xaxt="n")
  axis(3, at=c(1:length(names(sumRes)),28), labels=gsub("_", " ", c(gsub("_", " ", names(sumRes)), "Cumulative")), las=3, cex=0.5)
  polygon(c(0,0,28,28), c(0.7, 0.5, 0.5 ,0.7), col=makeTransparent("gray80", 100), border="white")
  polygon(c(0,0,28,28), c(2.2, 2.4, 2.4 ,2.2), col=makeTransparent("gray80", 100), border="white")
  lapply(1:length(sumRes), function(x) polygon(c(x, rep(x,37)-rev(sumRes[[x]])/width,x, rep(x,37)+(sumRes[[x]])/width), c(c(0,allcstime), c(2.7, rev(allcstime))), border ="gray90",col=cols[x]))
  #lapply(1:length(sumRes), function(x) points(rep(x,37), seq(0.1,3.7,0.1), pch="|", col=dlnLColors(rev(sumRes[[x]]$slnL), "red")))
  lapply(1:length(sumRes), function(x) points(x, maxTimes[x], pch=21, col=dlnLColors(dlnLs[x],"black", max=7), bg=cols[x]))
  allcs <- apply(do.call(cbind, sumRes), 1, sum)
  allcs4 <- allcs - (max(allcs)-6)
  allcs4 <- ifelse(allcs4 < 0, 0, allcs4)
  #oldcumsum <- apply(sapply(sumRes, function(x) x$slnL),1,sum)
  #cs <-rev(oldcumsum)-20
  #cs[cs<0] <- 0
  polygon(c(28, rep(28,37)-rev(allcs4)/8, 28, rep(28,37)+allcs4/8), c(c(0,allcstime), c(2.7, rev(allcstime))), border = col,col=col)
}
shiftFigure(sumRes, width=10)
dev.off()


## Figures showing individual shifts
pdf("../output/ShiftMagnitudes.pdf", height=16, width=24)
o <- order(sapply(allRes, function(x) which(apply(x, 1, median)==max(apply(x,1,median)))), decreasing=TRUE)
o <- o[which(o!=2)]
par(ask=FALSE, mfrow=c(5,5), mar=c(1,1,2,0))
Rs <- do.call(rbind, lapply(sumprof, function(x) do.call(rbind, x)))
qRs <- apply(Rs[,2:3], 2,function(x) log(quantile(x, c(0.025, 0.975), na.rm=TRUE)))
pal <- colorRampPalette(c("#c0392b", "#3498db"))
for(i in o){
  cs <- apply(allRes[[i]], 1, mean)
  maxIndex <- which(cs==max(cs))
  maxTime <- allcstime[maxIndex]
  maxCS <- max(cs)
  rs <- apply(t(sapply(sumprof, function(x) unlist(x[[i]][maxIndex,2:3]))), 2, median)
  lwds <- abs(log(rs))*1.2+0.5
  cols <- pal(100)[floor(97*(log(rs)+7)/17)+1]
  plot(allcstime, cs, ylim=c(0, max(cs)+0.2*max(cs)), xlim=c(2.6,0),xaxt="n",xlab="",ylab="", yaxt="n", type="n", main=gsub("_", " ", newnames[i]), cex.main=2.5)
  abline(v=c(0.6, 1.7), lty=2, col="red", lwd=3)
  abline(v=seq(3,0,-1), lty=3, col="gray50")
  abline(h=seq(0,20, 1), lty=3, col="gray50")
  abline(h=2, col="red", lty=2)
  lapply(1:ncol(allRes[[i]]), function(x) lines(seqs[[x]], allRes[[i]][,x], col="gray80"))
  lines(allcstime, cs, lty=1, lwd=3, col="black")
  arrows(maxTime-0.1, maxCS-0.1*maxCS, maxTime-0.1, maxCS+0.1*maxCS, code=ifelse(rs[1] < 1, 2, 1), lwd=lwds[1], col=cols[1])
  arrows(maxTime+0.1, maxCS-0.1*maxCS, maxTime+0.1, maxCS+0.1*maxCS, code=ifelse(rs[2] < 1, 2, 1), lwd=lwds[2], col=cols[2])
}

dev.off()


## Ancestral state reconstruction figures
pdf("../output/ASRfigures.pdf", height=11, width=8.5)
par(mfrow=c(2,2))
for(t in 1:20){
  cols <- c("#3498db", "#e74c3c")
  i=1
  pars.sum <- lapply(sumprof, function(tr) lapply(tr, function(trait) cbind(trait[which(trait$value==max(trait$value)),2:5],time=seq(0.1,2.7, length.out=37)[which(trait$value==max(trait$value))])))
  tr <- Fns[[i]]$tdList[[t]]$phy
  dd <- Fns[[i]]$tdList[[t]]$dat
  dd <- data.frame(dd[[1]])
  colnames(dd)[1] <- names(dat)[t]
  rownames(dd) <- tr$tip.label
  plot.asrR2(unlist(pars.sum[[i]][[t]][1,]), tr, dd, Fns[[i]]$ARD.R2[[t]], show.tip.label=FALSE, cols=cols)
  time <- round(seq(0.1,2.7, length.out=37)[which(sumprof[[i]][[t]]$value==max(sumprof[[i]][[t]]$value))],2)
  tiplabels(pch=21, bg=cols[dd[,1]+1])
  text(2.7-time,2, time, pos=1, font=2)
}
dev.off()

pdf("../output/Morphology&Habitat.pdf", height=10, width=6)
## Morphology
cols <- c("#3498db", "#e74c3c")
i=1
t = 7
pars.sum <- lapply(sumprof, function(tr) lapply(tr, function(trait) cbind(trait[which(trait$value==max(trait$value)),2:5],time=seq(0.1,2.7, length.out=37)[which(trait$value==max(trait$value))])))
tr <- Fns[[i]]$tdList[[t]]$phy
dd <- Fns[[i]]$tdList[[t]]$dat
dd <- data.frame(dd[[1]])
colnames(dd)[1] <- names(dat)[t]
rownames(dd) <- tr$tip.label
plot.asrR2(unlist(pars.sum[[i]][[t]][1,]), tr, dd, Fns[[i]]$ARD.R2[[t]], show.tip.label=FALSE, cols=cols)
time <- round(seq(0.1,2.7, length.out=37)[which(sumprof[[i]][[t]]$value==max(sumprof[[i]][[t]]$value))],1)
tiplabels(pch=21, bg=cols[dd[,1]+1])
text(2.7-time, 1, time, pos=4, font=2)

## Freshwater
cols <- c("#3498db", "#e74c3c")
i = 1
t = 3
pars.sum <- lapply(sumprof, function(tr) lapply(tr, function(trait) cbind(trait[which(trait$value==max(trait$value)),2:5],time=seq(0.1,2.7, length.out=37)[which(trait$value==max(trait$value))])))
tr <- Fns[[i]]$tdList[[t]]$phy
dd <- Fns[[i]]$tdList[[t]]$dat
dd <- data.frame(dd[[1]])
colnames(dd)[1] <- names(dat)[t]
rownames(dd) <- tr$tip.label
plot.asrR2(unlist(pars.sum[[i]][[t]][1,]), tr, dd, Fns[[i]]$ARD.R2[[t]], show.tip.label=FALSE, cols=cols)
time <- round(seq(0.1,2.7, length.out=37)[which(sumprof[[i]][[t]]$value==max(sumprof[[i]][[t]]$value))],2)
tiplabels(pch=21, bg=cols[dd[,1]+1])
text(2.7-time,2, time, pos=1, font=2)

dev.off()