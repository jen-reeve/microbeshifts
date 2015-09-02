\documentclass{article}

\begin{document}
\SweaveOpts{concordance=TRUE}

\section{Evolution of oxygen metabolism in freshwater and marine habitats in the Archaea}
<<echo=FALSE, fig=FALSE>>=
##Load packages and prepare environment
require(phytools)
require(diversitree)
require(geiger)
setwd("~/Dropbox/Microbes")
##Load data
treNSOA <- read.tree("NSOA_treeClean.tre")
datNSOA <- read.table("NSOA_char.txt",row.names=1)
datNSOA[datNSOA=="?"] <- NA
colnames(datNSOA) <- c("Temp", "pH", "Habitat", "Meth", "SulfRed", "SulfOx", "SulfideOx", "SulfThioRed", "NifH", "NitRed", "O2Met")
tdNSOA <- treedata(treNSOA, datNSOA)
tdNSOA$dat <- tdNSOA$dat[tdNSOA$phy$tip.label,]
##Some functions to use later
na.rm.phy <- function(td, traits){
  dat <- td$dat[,traits]
  nas <- apply(dat, 1, function(x) sum(is.na(x)))
  nas.which <- which(nas > 0)
  if(length(nas.which) > 0){
    dat <- td$dat[-nas.which,]
    phy <- drop.tip(td$phy, td$phy$tip.label[nas.which])
  }
  return(list(phy=phy, dat=dat[phy$tip.label,]))
}
as.numeric.n <- function(x){
  nn <- names(x)
  y <- as.numeric(x)
  names(y) <- nn
  return(y)
}
@

\begin{figure}
<<echo=FALSE, fig=TRUE>>=
plot(c(0,12), c(0, 60), type="n", bty="n", xaxt="n", yaxt="n", xlab="",ylab="")
plotTree(tdNSOA$phy, fsize=0.5,add=TRUE)
for(i in 1:ncol(tdNSOA$dat)){
  tiplabels(pch=22, bg=1+as.numeric(as.factor(tdNSOA$dat[,i])),adj=c(5+i/3, 0.5))
  text(5+i/3-2/5+max(nodeHeights(tdNSOA$phy)),length(tdNSOA$phy$tip.l)+4, labels=colnames(datNSOA)[i], srt=75, cex=0.5)
}

@
\caption{Distribution of traits across the Archaeal phylogeny}
\label{NSOAfig}
\end{figure}

<<echo=FALSE, fig=FALSE>>=
##Prepare habitat and O2Met data
##Coding plastic phenotypes as freshwater, coding hypersaline as marine....seems like a bad idea, but first pass
tdHabO2 <- na.rm.phy(tdNSOA, c("Habitat", "O2Met"))
tdHabO2$dat[grep("&", tdHabO2$dat[,'Habitat']),'Habitat'] <- "0"
tdHabO2$dat[grep("2", tdHabO2$dat[,'Habitat']),'Habitat'] <- "1"
tdHabO2$dat <- tdHabO2$dat[,c("Habitat", "O2Met")]
tdHabO2$dat <- as.data.frame(tdHabO2$dat)
tdHabO2$dat[,1:2] <- apply(tdHabO2$dat, 2, as.numeric)
colnames(tdHabO2$dat) <- c("H", "O")
tdHabO2$phy <- multi2di(tdHabO2$phy)
mkn.lik <- make.mkn.multitrait(tdHabO2$phy, tdHabO2$dat)
##Start with ML
mkn.lik.c <- constrain(mkn.lik, qO01.0 ~ qO01.H, qO10.0 ~ qO10.H)
HabO2.mle <- find.mle(mkn.lik, rep(0.1, 8))
HabO2.mle.c <- find.mle(mkn.lik.c, rep(0.1, 6))
#round(HabO2.mle$par,5)
#round(HabO2.mle.c$par,5)
#anova(HabO2.mle, HabO2.mle.c)

mkn.lik.c.H <- constrain(mkn.lik, qH01.0 ~ qH01.O, qH10.0 ~ qH10.O)
HabO2.mle.c.H <- find.mle(mkn.lik.c.H, rep(0.1, 6))
#anova(HabO2.mle, HabO2.mle.c.H)

mkn.lik.c.H.ER <- constrain(mkn.lik, qH01.0 ~ qH01.O, qH10.0 ~ qH10.O, qO01.0 ~ qO10.0)
HabO2.mle.c.H.ER <- find.mle(mkn.lik.c.H.ER, rep(0.1, 5))
#anova(HabO2.mle.c.OH, HabO2.mle.c.H.ER)

mkn.lik.c.OH <- constrain(mkn.lik, qH01.0 ~ qH01.O, qH10.0 ~ qH10.O, qO01.0 ~ qO01.H, qO10.0 ~ qO10.H)
HabO2.mle.c.OH <- find.mle(mkn.lik.c.OH, rep(0.1, 4))
#anova(HabO2.mle, HabO2.mle.c.OH)

#nsteps=1000
#chain.HabO2 <- mcmc(mkn.lik, x.init=rep(0.1,8), nsteps=nsteps, start=floor(0.3*nsteps), end=nsteps, thin=10, w=1, print.every=50, save.file="chain.HabO2.csv")
#read.csv("chain.HabO2.csv")
@


<<echo=FALSE, fig=FALSE>>=
#round(HabO2.mle$par, 5)
#round(HabO2.mle.c$par, 5)
#round(HabO2.mle.c.H$par,5 )
#round(HabO2.mle.c.OH$par,5)


#anova(HabO2.mle, HabO2.mle.c)
#anova(HabO2.mle, HabO2.mle.c.H)
#anova(HabO2.mle, HabO2.mle.c.OH)

#nsteps=1000
#chain.HabO2 <- mcmc(mkn.lik, x.init=rep(0.1,8), nsteps=nsteps, start=floor(0.3*nsteps), end=nsteps, thin=10, w=1, print.every=50, save.file="chain.HabO2.csv")
#read.csv("chain.HabO2.csv")
@


No models in which the transition rates for Oxygen metabolism depended upon habitat were significantly better than models in which the two are independent (using likelihood ratio tests). Large effect sizes were often estimated, but the method has relatively little power with this many taxa (\Sexpr{length(tdHabO2$phy$tip.label)}). Will do more power analyses...

%\section{Power analysis for finding a relationship between O2met and Habitat}


%<<echo=FALSE, fig=FALSE>>=
%#habitat <- tdHabO2$dat[,1]
%#names(habitat) <- rownames(tdHabO2$dat)
%#mkn.hab <- make.mk2(tdHabO2$phy, habitat)
%#find.mle(mkn.hab, c(1,1))
%PM=NULL
%for(r in 1/(c(1, 10, 100))){
%PP=NULL%

  %for(i in 1:50){
%nn=100
%tree <- sim.bdtree(1,d=0, n=nn)
%rescale.factor <- max(nodeHeights(tree))/max(nodeHeights(tdNSOA$phy))
%tree$edge.length <- tree$edge.length/rescale.factor
%Qt1 <- list(rbind(c(-1.5, 1.5), c(1.5, -1.5)))
%Qt2 <- list(rbind(c(-0.458, .458), c(.458, -0.458)))
%%rownames(Qt1[[1]]) <- c(1,2)
%colnames(Qt1[[1]]) <- c(1,2)
%Dt1 <- sim.char(tree, Qt1, model="discrete", n=1)[,,1]
%true.history <- make.simmap(tree, Dt1, Q=Qt1[[1]])
%#plotSimmap(true.history)
%Dt2 <- sim.char(map.expand(r, true.history, state=2), Qt2, model="discrete", n=1)[,,1]
%
%#naive.FD <- fitDiscrete(tree, Dt2, model="ARD")
%#true.FD <- fitDiscrete(map.expand(r, true.history, state=2), Dt2, model="ARD")
%#plotSimmap(true.history,ftype="off")
%#tiplabels(pch=21, bg=Dt2+2, col=Dt2+2, cex=0.5)
%#tiplabels(pch=21, bg=Dt1+2, col=Dt1+2, cex=0.5, adj=c(1,0.5))

%sim.dat <- data.frame("A"=Dt1-1, "B"=Dt2-1)
%mkn.tmp <- make.mkn.multitrait(tree, sim.dat)
%full.model <- constrain(mkn.tmp, qA01.0 ~ qA01.B, qA10.0 ~ qA10.B)
%ind.model <- constrain(mkn.tmp, qA01.0 ~ qA01.B, qA10.0 ~ qA10.B, qB01.0 ~ qB01.A, qB10.0 ~ qB10.A)
%mle1<- (find.mle(full.model, rep(1, 6)))
%mle2<- (find.mle(ind.model, rep(1, 4)))
%PP <- c(PP, anova(mle1, mle2)$Pr[2])
%}
%PM <- c(PM, sum(PP<0.05))
%}
%HabO2.mle.c$par

%@

\section{Evolution of oxygen metabolism in different temperature regimes in the Archaea}
<<echo=FALSE, fig=FALSE>>=
##Hyperthermophiles vs. O2 metabolism or Sulfur reduction
##Coding plastic phenotypes as hyperthermophile or not
tdTempO2 <- na.rm.phy(tdNSOA, c("Temp", "O2Met"))
tdTempO2$dat[grep("1", tdTempO2$dat[,'Temp']),'Temp'] <- "0"
tdTempO2$dat[grep("2", tdTempO2$dat[,'Temp']),'Temp'] <- "1"
tdTempO2$dat <- tdTempO2$dat[,c("Temp", "O2Met")]
tdTempO2$dat <- as.data.frame(tdTempO2$dat)
tdTempO2$dat[,1:2] <- apply(tdTempO2$dat, 2, as.numeric)
colnames(tdTempO2$dat) <- c("T", "O")
tdTempO2$phy <- multi2di(tdTempO2$phy)
mkn.lik <- make.mkn.multitrait(tdTempO2$phy, tdTempO2$dat)
##Start with ML
mkn.lik.c <- constrain(mkn.lik, qO01.0 ~ qO01.T, qO10.0 ~ qO10.T)
TempO2.mle <- find.mle(mkn.lik, rep(0.1, 8))
TempO2.mle.c <- find.mle(mkn.lik.c, rep(0.1, 6))
mkn.lik.c.T <- constrain(mkn.lik, qT01.0 ~ qT01.O, qT10.0 ~ qT10.O)
TempO2.mle.c.T <- find.mle(mkn.lik.c.T, rep(0.1, 6))
mkn.lik.c.OT <- constrain(mkn.lik, qT01.0 ~ qT01.O, qT10.0 ~ qT10.O, qO01.0 ~ qO01.T, qO10.0 ~ qO10.T)
TempO2.mle.c.OT <- find.mle(mkn.lik.c.OT, rep(0.1, 4))
@

<<echo=FALSE, fig=FALSE>>=
#round(TempO2.mle$par, 5)
#round(TempO2.mle.c$par, 5 )
#round(TempO2.mle.c.T$par, 5)
#round(TempO2.mle.c.OT$par,5)


#anova(TempO2.mle, TempO2.mle.c)
#anova(TempO2.mle, TempO2.mle.c.T)
#anova(TempO2.mle, TempO2.mle.c.OT)

#nsteps=1000
#chain.HabO2 <- mcmc(mkn.lik, x.init=rep(0.1,8), nsteps=nsteps, start=floor(0.3*nsteps), end=nsteps, thin=10, w=1, print.every=50, save.file="chain.HabO2.csv")
#read.csv("chain.HabO2.csv")
@

Same results as for habitat...

<<echo=FALSE, fig=FALSE>>=
##Hyperthermophiles vs. Sulfur reduction
##Coding plastic phenotypes as hyperthermophile or not
tdTempSR <- na.rm.phy(tdNSOA, c("Temp", "SulfRed"))
tdTempSR$dat[grep("1", tdTempSR$dat[,'Temp']),'Temp'] <- "0"
tdTempSR$dat[grep("2", tdTempSR$dat[,'Temp']),'Temp'] <- "1"
tdTempSR$dat <- tdTempSR$dat[,c("Temp", "SulfRed")]
tdTempSR$dat <- as.data.frame(tdTempSR$dat)
tdTempSR$dat[,1:2] <- apply(tdTempSR$dat, 2, as.numeric)
colnames(tdTempSR$dat) <- c("T", "S")
tdTempSR$phy <- multi2di(tdTempSR$phy)
mkn.lik <- make.mkn.multitrait(tdTempSR$phy, tdTempSR$dat)
##Start with ML
mkn.lik.c <- constrain(mkn.lik, qS01.0 ~ qS01.T, qS10.0 ~ qS10.T)
TempSR.mle <- find.mle(mkn.lik, rep(0.1, 8))
TempSR.mle.c <- find.mle(mkn.lik.c, rep(0.1, 6))
mkn.lik.c.T <- constrain(mkn.lik, qT01.0 ~ qT01.S, qT10.0 ~ qT10.S)
TempSR.mle.c.T <- find.mle(mkn.lik.c.T, rep(0.1, 6))
mkn.lik.c.ST <- constrain(mkn.lik, qT01.0 ~ qT01.S, qT10.0 ~ qT10.S, qS01.0 ~ qS01.T, qS10.0 ~ qS10.T)
TempSR.mle.c.ST <- find.mle(mkn.lik.c.ST, rep(0.1, 4))
@

<<echo=FALSE, fig=FALSE>>=
#round(TempSR.mle$par, 5)
#round(TempSR.mle.c$par, 5 )
#round(TempSR.mle.c.T$par, 5)
#round(TempSR.mle.c.ST$par,5)


#anova(TempSR.mle, TempSR.mle.c)
#anova(TempSR.mle, TempSR.mle.c.T)
#anova(TempSR.mle, TempSR.mle.c.ST)

#nsteps=1000
#chain.HabO2 <- mcmc(mkn.lik, x.init=rep(0.1,8), nsteps=nsteps, start=floor(0.3*nsteps), end=nsteps, thin=10, w=1, print.every=50, save.file="chain.HabO2.csv")
#read.csv("chain.HabO2.csv")
@

%\section{Sulfur Reduction vs. NifH gene presence/absence}
<<echo=FALSE, fig=FALSE>>=
##NifH vs. Sulfur reduction
##Coding plastic phenotypes as hyperthermophile or not
tdNifHSR <- na.rm.phy(tdNSOA, c("NifH", "SulfRed"))
#tdNifHSR$dat[grep("1", tdNifHSR$dat[,'NifH']),'NifH'] <- "0"
#tdNifHSR$dat[grep("2", tdNifHSR$dat[,'NifH']),'NifH'] <- "1"
tdNifHSR$dat <- tdNifHSR$dat[,c("NifH", "SulfRed")]
tdNifHSR$dat <- as.data.frame(tdNifHSR$dat)
tdNifHSR$dat[,1:2] <- apply(tdNifHSR$dat, 2, as.numeric)
colnames(tdNifHSR$dat) <- c("N", "S")
tdNifHSR$phy <- multi2di(tdNifHSR$phy)
mkn.lik <- make.mkn.multitrait(tdNifHSR$phy, tdNifHSR$dat)
##Start with ML
mkn.lik.c <- constrain(mkn.lik, qS01.0 ~ qS01.N, qS10.0 ~ qS10.N)
NifHSR.mle <- find.mle(mkn.lik, rep(0.1, 8))
NifHSR.mle.c <- find.mle(mkn.lik.c, rep(0.1, 6))
mkn.lik.c.N <- constrain(mkn.lik, qN01.0 ~ qN01.S, qN10.0 ~ qN10.S)
NifHSR.mle.c.N <- find.mle(mkn.lik.c.N, rep(0.1, 6))
mkn.lik.c.SN <- constrain(mkn.lik, qN01.0 ~ qN01.S, qN10.0 ~ qN10.S, qS01.0 ~ qS01.N, qS10.0 ~ qS10.N)
NifHSR.mle.c.SN <- find.mle(mkn.lik.c.SN, rep(0.1, 4))

#full.model <- constrain(mkn.lik, qN01.0 ~ qN01.S, qN10.0 ~ qN10.S)
#ind.model <- constrain(mkn.lik, qN01.0 ~ qN01.S, qN10.0 ~ qN10.S, qS01.0 ~ qS01.N, qS10.0 ~ qS10.N)
#mle1<- (find.mle(full.model, rep(1, 6)))
#mle2<- (find.mle(ind.model, rep(1, 4)))
#anova(mle1, mle2)
#round(mle1$par,5)
#round(mle2$par,5)
@

<<echo=FALSE, fig=FALSE>>=
#round(NifHSR.mle$par, 5)
##round(NifHSR.mle.c$par, 5 )
#round(NifHSR.mle.c.N$par, 5)
#round(NifHSR.mle.c.SN$par,5)


#anova(NifHSR.mle, NifHSR.mle.c)
#anova(NifHSR.mle, NifHSR.mle.c.N)
#anova(NifHSR.mle, NifHSR.mle.c.SN)

#nsteps=1000
#chain.HabO2 <- mcmc(mkn.lik, x.init=rep(0.1,8), nsteps=nsteps, start=floor(0.3*nsteps), end=nsteps, thin=10, w=1, print.every=50, save.file="chain.HabO2.csv")
#read.csv("chain.HabO2.csv")
@

\section{Time Slices- Before or after Oxygenation event?}

<<echo=FALSE, fig=FALSE>>=
  tr <- multi2di(tdNSOA$phy)
seq.x <- seq(0, max(branching.times(tr)), length.out=100)
bd.lik <- make.bd(tr)
  b.lik <- constrain(bd.lik, mu~0)
  fit.b <- find.mle(b.lik, 1)
  lambda.est <- fit.b$par
  mu.est <- 0
  dat <- as.numeric.n(tdNSOA$dat[,'O2Met'])
  e1.lik <- make.bisse(tr,dat)
  e2.lik <- make.bisse.t(tr, dat, functions=c(rep("constant.t",4), rep("stepf.t", 2)))
  e3.lik <- make.bisse.t(tr, dat, functions=c(rep("constant.t",4), rep("stepf.t", 2)))
  #argnames(e2.lik)
  grO2event <- max(branching.times(tr))-2.32
  secO2event <- max(branching.times(tr))-(0.663+0.551)/2
  O2.lik <- list()
  O2.lik[[1]] <- constrain(e1.lik, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, q01~q10)
  O2.lik[[2]] <- constrain(e2.lik, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, q01.y0~q10.y0, q01.y1~q10.y1, 
                           q01.tc~grO2event, q10.tc~grO2event)
  O2.lik[[3]] <- constrain(e3.lik, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, q01.y0~q10.y0, q01.y1~q10.y1, 
                       q01.tc~secO2event, q10.tc~secO2event)
  O2.lik[[4]] <- constrain(e1.lik, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0)
  O2.lik[[5]] <- constrain(e2.lik, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                           q01.tc~grO2event, q10.tc~grO2event)
  O2.lik[[6]] <- constrain(e3.lik, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                           q01.tc~secO2event, q10.tc~secO2event)
  
  O2.fits <- lapply(1:length(O2.lik),function(x) find.mle(O2.lik[[x]], rep(1, length(argnames(O2.lik[[x]])))))
  O2.lnLiks <- sapply(O2.fits, function(x) x$lnLik)
  names(O2.fits) <- c("ER", "ER-GO2E", "ER-SO2E", "ARD","ARD-GO2E", "ARD-SO2E")
@
<<echo=FALSE, fig=FALSE>>=
  ##Akaike scores
  #sort(sapply(O2.fits, AIC))
  ##All Primary o2 events LRTs
  #anova( O2.fits[[1]],  O2.fits[[2]] , O2.fits[[4]], O2.fits[[5]])
  ##All Secondary o2 events LRTs
  #anova( O2.fits[[1]], O2.fits[[2]], O2.fits[[3]], O2.fits[[6]])
  #sapply(O2.fits, function(x) round(x$par,5))
@

There appears no evidence for a significant shift in rates of Oxygen metabolism at any specific in time in the last 3.8 by. Examining all traits, only one shows a significant shift (SulfThioRed). However, looking at the plots for each of them individually suggests that there may be some signal for all traits combined for either the great or secondary oxygenation events. The following figures show estimates for combined traits and each trait individually. In current implementations, the models are too complex to favor a shift (adding a time shift doubles the number of parameters and is very costly), but there is no reason why this must be so. It will just require some tinkering to fit simpler models that aren't so parameter rich. 

<<echo=FALSE, fig=TRUE>>=
traits <-  c("NitRed", "O2Met", "SulfOx", "SulfideOx")
tdOx <- na.rm.phy(tdNSOA,traits)
tr <- multi2di(tdOx$phy)
dat <- as.data.frame(tdOx$dat[,traits])
dat[,1:4] <- apply(dat, 2, as.numeric)
bd.lik <- make.bd(tr)
b.lik <- constrain(bd.lik, mu~0)
fit.b <- find.mle(b.lik, 1)
lambda.est <- fit.b$par
mu.est <- 0
list.dat <- as.list(dat)
list.dat <- lapply(list.dat, function(x){ names(x) <- rownames(dat); x})
bisse.fns <- lapply(list.dat, function(x) make.bisse.t(tr, x, functions=c(rep("constant.t",4), rep("stepf.t", 2))))
cons.fns <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                    q01.y0~q10.y0, q01.y1~q10.y1, q01.tc ~ q10.tc))
cons.fns.ARD <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                     q01.tc ~ q10.tc))
#cons.fns <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
#                                                   q01.y0~q01.y1 * 5, q10.y0~q10.y1 * 5, q01.tc ~ q10.tc))
.all.lik <- function(fxns, time, pars=list()){
    liks <- sapply(1:length(fxns), function(x) fxns[[x]](c(pars[[x]], time)))
    sum(liks)
}
all.lik <- function(pars, fxns=cons.fns){
  time <- pars[1]
  pars[2:length(pars)] <- exp(pars[2:length(pars)])
  new.pars <- list()
  lens <- sapply(fxns, function(x) length(argnames(x)))-1
  end <- 1
  for(i in 1:length(fxns)){
    start <- end+1
    end <- start+lens[i]-1
    new.pars[[i]] <- pars[start:end]
  }
  .all.lik(fxns, time, new.pars)
}
#.all.lik(cons.fns, 1, pars=list(c(1,1), c(1,1), c(1,1), c(1,1)))
#fit.Ox <- find.mle(all.lik, rep(1,9), method="nlminb", lower=rep(0, 9), upper=c(max(branching.times(tr)), rep(10,8)))
all.lik.ARD <- function(pars, fxns=cons.fns.ARD){
  time <- pars[1]
  pars[2:length(pars)] <- exp(pars[2:length(pars)])
  new.pars <- list()
  lens <- sapply(fxns, function(x) length(argnames(x)))-1
  end <- 1
  for(i in 1:length(fxns)){
    start <- end+1
    end <- start+lens[i]-1
    new.pars[[i]] <- pars[start:end]
  }
  .all.lik(fxns, time, new.pars)
}
#fit.Ox <- find.mle(all.lik, rep(1,9), method="optim", lower=c(0,rep(-500,8)), upper=c(max(branching.times(tr)), rep(500,8)))
#fit.Ox.ARD <- find.mle(all.lik.ARD, rep(1,17), method="optim",lower=c(0,rep(-20,8)), upper=c(max(branching.times(tr)), rep(10,8)))
#save(fit.Ox, file="fit.Ox.RData")
#save(fit.Ox.ARD, file="fit.Ox.ARD.RData")
load("fit.Ox.RData")
load("fit.Ox.ARD.RData")

lik.prof <- sapply(seq(0, max(branching.times(tr)), length.out=100), function(x) all.lik(c(x, fit.Ox$par[2:9])))
lik.prof.ARD <- sapply(seq(0, max(branching.times(tr)), length.out=100), function(x) all.lik.ARD(c(x, fit.Ox.ARD$par[2:17])))


ylim.min=-270
plot(seq(0, max(branching.times(tr)),length.out=100), lik.prof,type="n", ylim=c(ylim.min, -250),xaxt="n",ylab="Log Likelihood", xlab="Time before present (by)", main="NitRed, O2Met, SulfOx, SulfideOx")
axis(1, at=c(0, 1 ,2 ,3, 3.8), labels=c(3.8, 3, 2, 1, 0))
#lines(max(branching.times(tr))-seq(0, max(branching.times(tr)),length.out=100), lik.prof,col="red")
lines(max(branching.times(tr))-seq(0, max(branching.times(tr)),length.out=100), lik.prof.ARD, col="blue")
#lines(max(branching.times(tr))-seq(0, max(branching.times(tr)),length.out=100), lik.prof, col="red")
#abline(v=grO2event,lty=2,lwd=2)
text(grO2event, ylim.min, "Great Oxygenation Event",cex=0.7)
#abline(v=secO2event,lty=2,lwd=2)
text(secO2event, ylim.min, "Secondary Oxygenation Event",cex=0.7)
abline(v=max(branching.times(tr))-fit.Ox.ARD$par[1], col="red", lwd=2)
points(max(branching.times(tr))-fit.Ox.ARD$par[1], fit.Ox.ARD$lnLik,pch=21, bg="red")
rect(-2, max(lik.prof.ARD-2), 4, max(lik.prof.ARD), col=rgb(0,0,255, maxColorValue=255, alpha=50), border=NA)
rect(max(branching.times(tr))- 2.45, -500, max(branching.times(tr))-1.85, 0, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(max(branching.times(tr))- 0.663, -500, max(branching.times(tr))-0.551, 0, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
@

%\includegraphics{"Oxygenation.png"}

<<echo=FALSE, fig=TRUE>>=
##All data
traits <-  c("Habitat","Temp","NitRed", "O2Met", "SulfOx", "SulfideOx")
tdOx <- na.rm.phy(tdNSOA,traits)
tr <- multi2di(tdOx$phy)
dat <- as.data.frame(tdOx$dat[,traits])
dat[dat=="0&1"] <- "0"
dat[dat[,1]=="2",1] <- "1"
dat[dat[,2]=="2",2] <- "1"
dat[,1:ncol(dat)] <- apply(dat, 2, as.numeric)
bd.lik <- make.bd(tr)
b.lik <- constrain(bd.lik, mu~0)
fit.b <- find.mle(b.lik, 1)
lambda.est <- fit.b$par
mu.est <- 0
list.dat <- as.list(dat)
list.dat <- lapply(list.dat, function(x){ names(x) <- rownames(dat); x})
bisse.fns <- lapply(list.dat, function(x) make.bisse.t(tr, x, functions=c(rep("constant.t",4), rep("stepf.t", 2)))) 
cons.fns <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                    q01.y0~q10.y0, q01.y1~q10.y1, q01.tc ~ q10.tc))
cons.fns.ARD <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                     q01.tc ~ q10.tc))
#cons.fns <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
#                                                   q01.y0~q01.y1 * 5, q10.y0~q10.y1 * 5, q01.tc ~ q10.tc))
.all.lik <- function(fxns, time, pars=list()){
    liks <- sapply(1:length(fxns), function(x) fxns[[x]](c(pars[[x]], time)))
    sum(liks)
}
all.lik <- function(pars, fxns=cons.fns){
  time <- pars[1]
  pars[2:length(pars)] <- exp(pars[2:length(pars)])
  new.pars <- list()
  lens <- sapply(fxns, function(x) length(argnames(x)))-1
  end <- 1
  for(i in 1:length(fxns)){
    start <- end+1
    end <- start+lens[i]-1
    new.pars[[i]] <- pars[start:end]
  }
  .all.lik(fxns, time, new.pars)
}
all.lik.ARD <- function(pars, fxns=cons.fns.ARD){
  time <- pars[1]
  pars[2:length(pars)] <- exp(pars[2:length(pars)])
  new.pars <- list()
  lens <- sapply(fxns, function(x) length(argnames(x)))-1
  end <- 1
  for(i in 1:length(fxns)){
    start <- end+1
    end <- start+lens[i]-1
    new.pars[[i]] <- pars[start:end]
  }
  .all.lik(fxns, time, new.pars)
}
#fit.Ox2 <- find.mle(all.lik, rep(1,13), method="optim", lower=c(0,rep(-20,12)), upper=c(max(branching.times(tr)), rep(10,12)))

#ind.fits <- lapply(1:ncol(dat), function(x) find.mle(cons.fns[[x]], c(1,1,1)))
#ind.fits.ARD <- lapply(1:ncol(dat), function(x) find.mle(cons.fns.ARD[[x]], c(1,1,1,1,1)))

#fit.pars <- unlist(lapply(ind.fits.ARD,function(x) x$par[1:4]))
#fit.Ox.ARD <- find.mle(all.lik.ARD, c(2, fit.pars), method="optim", lower=rep(0, 4*6+1), upper=c(max(branching.times(tr)), rep(20,4*6)))
#fit.Ox2.ARD <- find.mle(all.lik.ARD, rep(1, 25), method="optim", lower=c(0,rep(-15,24)), upper=c(max(branching.times(tr)), rep(5,24)))
#save(fit.Ox2, file="fit.Ox2.RData")
#save(fit.Ox2.ARD, file="fit.Ox2.ARD.RData")
load("fit.Ox2.RData")
load("fit.Ox2.ARD.RData")
lik.prof <- sapply(seq(0, max(branching.times(tr)), length.out=100), function(x) all.lik(c(x, fit.Ox2$par[2:13])))
lik.prof.ARD <- sapply(seq(0, max(branching.times(tr)), length.out=100), function(x) all.lik.ARD(c(x, fit.Ox2.ARD$par[2:25])))

ylim.min=-410
plot(seq(0, max(branching.times(tr)),length.out=100), lik.prof,type="n",ylim=c(ylim.min, -395),xaxt="n",ylab="Log Likelihood", xlab="Time before present (by)", main="Hab, Temp, O2met, SulfOx, SulfideOx & NitRed")
axis(1, at=c(0, 1 ,2 ,3, 3.8), labels=c(3.8, 3, 2, 1, 0))
#lines(max(branching.times(tr))-seq(0, max(branching.times(tr)),length.out=100), lik.prof,col="red")
lines(max(branching.times(tr))-seq(0, max(branching.times(tr)),length.out=100), lik.prof.ARD, col="blue")
text(grO2event, ylim.min, "Great Oxygenation Event",cex=0.7)
text(secO2event, ylim.min, "Secondary Oxygenation Event",cex=0.7)
abline(v=max(branching.times(tr))-fit.Ox2.ARD$par[1], col="red", lwd=2)
points(max(branching.times(tr))-fit.Ox2.ARD$par[1], fit.Ox.ARD$lnLik,pch=21, bg="red")
rect(-2, max(lik.prof.ARD-2), 4, max(lik.prof.ARD), col=rgb(0,0,255, maxColorValue=255, alpha=50), border=NA)
rect(max(branching.times(tr))- 2.45, -500, max(branching.times(tr))-1.85, 0, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(max(branching.times(tr))- 0.663, -500, max(branching.times(tr))-0.551, 0, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
#rect(-2, max(lik.prof-2), 4, max(lik.prof), col=rgb(255,0,0, maxColorValue=255, alpha=50), border=NA)
#fit.pars <- sapply(1:ncol(dat), function(x) round(exp(fit.Ox.ARD$par[2:length(fit.Ox.ARD$par)]),2)[(x*1):(x*1+3)])
#fit.pars
@

%<<echo=FALSE, fig=TRUE>>=
%#ind.fits.ARD <- lapply(1:ncol(dat), function(x) find.mle(cons.fns.ARD[[x]], c(1,1,1,1,1)))
%ind.fits.pars <- lapply(ind.fits.ARD, function(x) x$par[1:4])
%ind.fits.times <- round(sapply(ind.fits.ARD, function(x) x$par[5]),2)
%abline(v=max(branching.times(tr))-ind.fits.times, col=2:8, lty=2)
%legend(0, ylim.min+5, legend=colnames(dat), lty=rep(2, ncol(dat)), col=2:(2+ncol(dat)),cex=0.75)
%round(sapply(ind.fits.ARD, function(x) x$par[1:4]),2)
%fit.pars

%rbind(sapply(ind.fits.ARD,function(x) x$lnLik), sapply(1:ncol(dat), function(x) cons.fns.ARD[[x]](c(fit.pars[[x]],fit.Ox.ARD$par[1]))))

%par(mfrow=c(2,3))
%for(y in 1:ncol(dat)){
%  #l1 <- sapply(seq.x, function(x) cons.fns.ARD[[y]](c(fit.pars[,y], x)))
%  l2 <- sapply(seq.x, function(x) cons.fns.ARD[[y]](c(ind.fits.pars[[y]], x)))
%  l1 <- l2
%  plot(max(branching.times(tr))-seq.x, l2, type="n", ylim=c(max(c(l1,l2))-5, max(c(l1,l2))),xaxt="n",xlab="Time before %present (by)", main=colnames(dat)[y])
%  axis(1, at=c(0, 1 ,2 ,3, 3.8), labels=c(3.8, 3, 2, 1, 0))
%  #lines(max(branching.times(tr))-seq.x, l1, col="blue")
%  lines(max(branching.times(tr))-seq.x, l2, col="red",lwd=2)
%  rect(max(branching.times(tr))- 2.45, -500, max(branching.times(tr))-1.85, 0, col=rgb(100,100,100, maxColorValue=255%,alpha=25), border=NA)
%rect(max(branching.times(tr))- 0.663, -500, max(branching.times(tr))-0.551, 0, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
%}

%@

<<echo=FALSE, fig=TRUE>>=
##All data
#traits <-  colnames(tdNSOA$dat)
#tdOx <- na.rm.phy(tdNSOA,traits)
#tr <- tdOx$phy
tr <- multi2di(tdNSOA$phy)
dat <- data.frame(tdNSOA$dat)#as.data.frame(tdOx$dat[,traits])
dat[dat=="0&1"] <- "0"
dat[which(dat[,1]=="2"),1] <- "1"
dat[which(dat[,2]=="2"),2] <- "1"
dat[which(dat[,3]=="2"),3] <- "1"
dat[,1:ncol(dat)] <- apply(dat, 2, as.numeric)
bd.lik <- make.bd(tr)
b.lik <- constrain(bd.lik, mu~0)
fit.b <- find.mle(b.lik, 1)
lambda.est <- fit.b$par
mu.est <- 0
list.dat <- as.list(dat)
list.dat <- lapply(list.dat, function(x){ names(x) <- rownames(dat); x})
bisse.fns <- lapply(list.dat, function(x) make.bisse.t(tr, x, functions=c(rep("constant.t",4), rep("stepf.t", 2)))) 
notime.fns <- lapply(list.dat, function(x) make.bisse(tr, x))
cons.fns <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                    q01.y0~q10.y0, q01.y1~q10.y1, q01.tc ~ q10.tc))
cons.notime.fns.ER <- lapply(notime.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, q01~q10))
cons.fns.ARD <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                     q01.tc ~ q10.tc))
cons.notime.fns.ARD <- lapply(notime.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0))

#ind.fits.ER <- lapply(1:ncol(dat), function(x) find.mle(cons.fns[[x]], c(1,1,1)))
#ind.fits.ARD <- lapply(1:ncol(dat), function(x) find.mle(cons.fns.ARD[[x]], c(1,1,1,1,1)))
#ind.fits.notime <- lapply(1:ncol(dat), function(x) find.mle(cons.notime.fns.ARD[[x]], c(1,1)))
#ind.fits.notime.ER <- lapply(1:ncol(dat), function(x) find.mle(cons.notime.fns.ER[[x]], c(1)))
#save(ind.fits.ER, file="ind.fits.ER.RData")
#save(ind.fits.ARD, file="ind.fits.ARD.RData")
#save(ind.fits.notime, file="ind.fits.notime.RData")
#save(ind.fits.notime.ER, file="ind.fits.notime.ER.RData")
load(file="ind.fits.ER.RData")
load(file="ind.fits.ARD.RData")
load(file="ind.fits.notime.RData")
load(file="ind.fits.notime.ER.RData")
AIC.table <- rbind(sapply(ind.fits.ER, AIC), sapply(ind.fits.ARD, AIC), sapply(ind.fits.notime, AIC), sapply(ind.fits.notime.ER, AIC))
colnames(AIC.table) <- (colnames(dat))
rownames(AIC.table) = c('ER + time', 'ARD + time', 'ARD', 'ER')
#AIC.table
#lapply(1:ncol(dat), function(x) anova(ind.fits.ARD[[x]], ind.fits.notime[[x]]))
#lapply(1:ncol(dat), function(x) anova(ind.fits.ER[[x]], ind.fits.notime.ER[[x]]))
ind.fits.pars <- lapply(ind.fits.ARD, function(x) x$par[1:4])
ind.fits.times <- round(sapply(ind.fits.ARD, function(x) x$par[5]),2)
#pdf("IndividualTimeSplitsFits.pdf")
par(mfrow=c(2,3))
for(y in 1:6){
  #l1 <- sapply(seq.x, function(x) cons.fns.ARD[[y]](c(fit.pars[,y], x)))
  l2 <- sapply(seq.x, function(x) cons.fns.ARD[[y]](c(ind.fits.pars[[y]], x)))
  l1 <- l2
  plot(max(branching.times(tr))-seq.x, l2, type="n", ylim=c(max(c(l1,l2))-5, max(c(l1,l2))),xaxt="n",xlab="Time before present (by)", main=colnames(dat)[y])
  axis(1, at=c(0, 0.8 ,1.8 ,2.8, 3.8), labels=c(3.8, 3, 2, 1, 0))
  #lines(max(branching.times(tr))-seq.x, l1, col="blue")
  lines(max(branching.times(tr))-seq.x, l2, col="red",lwd=2)
  rect(max(branching.times(tr))- 2.45, -500, max(branching.times(tr))-1.85, 0, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(max(branching.times(tr))- 0.663, -500, max(branching.times(tr))-0.551, 0, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
  abline(v=max(branching.times(tr))-seq.x[which(l2==max(l2))],col="red", lty=2)
}
#dev.off()
@

<<echo=FALSE, fig=TRUE>>=
par(mfrow=c(2,3))
for(y in 7:11){
  #l1 <- sapply(seq.x, function(x) cons.fns.ARD[[y]](c(fit.pars[,y], x)))
  l2 <- sapply(seq.x, function(x) cons.fns.ARD[[y]](c(ind.fits.pars[[y]], x)))
  l1 <- l2
  plot(max(branching.times(tr))-seq.x, l2, type="n", ylim=c(max(c(l1,l2))-5, max(c(l1,l2))),xaxt="n",xlab="Time before present (by)", main=colnames(dat)[y])
  axis(1, at=c(0, 0.8 ,1.8 ,2.8, 3.8), labels=c(3.8, 3, 2, 1, 0))
  #lines(max(branching.times(tr))-seq.x, l1, col="blue")
  lines(max(branching.times(tr))-seq.x, l2, col="red",lwd=2)
  rect(max(branching.times(tr))- 2.45, -500, max(branching.times(tr))-1.85, 0, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(max(branching.times(tr))- 0.663, -500, max(branching.times(tr))-0.551, 0, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
  abline(v=max(branching.times(tr))-seq.x[which(l2==max(l2))],col="red", lty=2)
}
@

Direction of changes:

<<echo=FALSE, fig=TRUE>>=
names(ind.fits.pars) <- colnames(dat)
rate.changes <- lapply(ind.fits.pars, function(x) matrix(x, ncol=2, nrow=2))
par(mfrow=c(3,4))
tmp <- lapply(1:11,function(x){ image(rate.changes[[x]], main=names(rate.changes)[x], col=cm.colors(12), xaxt="n", yaxt="n");
                         axis(1, at=c(0, 1), labels=c("T0", "T1"));
                         axis(2, at=c(0, 1), labels=c("0->1","1->0"), las=2)
                         })
@


<<echo=FALSE, fig=TRUE>>=
#cons.fns <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
#                                                   q01.y0~q01.y1 * 5, q10.y0~q10.y1 * 5, q01.tc ~ q10.tc))
.all.lik <- function(fxns, time, pars=list()){
    liks <- sapply(1:length(fxns), function(x) fxns[[x]](c(pars[[x]], time)))
    sum(liks)
}
all.lik <- function(pars, fxns=cons.fns){
  time <- pars[1]
  pars[2:length(pars)] <- exp(pars[2:length(pars)])
  new.pars <- list()
  lens <- sapply(fxns, function(x) length(argnames(x)))-1
  end <- 1
  for(i in 1:length(fxns)){
    start <- end+1
    end <- start+lens[i]-1
    new.pars[[i]] <- pars[start:end]
  }
  .all.lik(fxns, time, new.pars)
}
all.lik.ARD <- function(pars, fxns=cons.fns.ARD){
  time <- pars[1]
  pars[2:length(pars)] <- exp(pars[2:length(pars)])
  new.pars <- list()
  lens <- sapply(fxns, function(x) length(argnames(x)))-1
  end <- 1
  for(i in 1:length(fxns)){
    start <- end+1
    end <- start+lens[i]-1
    new.pars[[i]] <- pars[start:end]
  }
  .all.lik(fxns, time, new.pars)
}
#fit.Ox3 <- find.mle(all.lik, rep(1,23), method="optim", lower=c(0,rep(-15,22)), upper=c(max(branching.times(tr)), rep(5,22)))

#ind.fits <- lapply(1:ncol(dat), function(x) find.mle(cons.fns[[x]], c(1,1,1)))
#ind.fits.ARD <- lapply(1:ncol(dat), function(x) find.mle(cons.fns.ARD[[x]], c(1,1,1,1,1)))

#fit.pars <- unlist(lapply(ind.fits.ARD,function(x) x$par[1:4]))
#fit.Ox.ARD <- find.mle(all.lik.ARD, c(2, fit.pars), method="optim", lower=rep(0, 4*6+1), upper=c(max(branching.times(tr)), rep(20,4*6)))
#fit.Ox3.ARD <- find.mle(all.lik.ARD, rep(1,45), method="optim", lower=c(0,rep(-10,44)), upper=c(max(branching.times(tr)), rep(3,44)))
#save(fit.Ox3, file="fit.Ox3.RData")
#save(fit.Ox3.ARD, file="fit.Ox3.ARD.RData")
load("fit.Ox3.RData")
load("fit.Ox3.ARD.RData")
lik.prof <- sapply(seq(0, max(branching.times(tr)), length.out=100), function(x) all.lik(c(x, fit.Ox3$par[2:23])))
lik.prof.ARD <- sapply(seq(0, max(branching.times(tr)), length.out=100), function(x) all.lik.ARD(c(x, fit.Ox3.ARD$par[2:45])))
ylim.min=-760
plot(seq(0, max(branching.times(tr)),length.out=100), lik.prof,type="n",ylim=c(ylim.min, -750),xaxt="n",ylab="Log Likelihood", xlab="Time before present (by)", main="All traits, one time shift")
axis(1, at=c(0, 1 ,2 ,3, 3.8), labels=c(3.8, 3, 2, 1, 0))
#lines(max(branching.times(tr))-seq(0, max(branching.times(tr)),length.out=100), lik.prof,col="red")
lines(max(branching.times(tr))-seq(0, max(branching.times(tr)),length.out=100), lik.prof.ARD, col="blue")
text(grO2event, ylim.min, "Great Oxygenation Event",cex=0.7)
text(secO2event, ylim.min, "Secondary Oxygenation Event",cex=0.7)
abline(v=max(branching.times(tr))-fit.Ox3.ARD$par[1], col="red", lwd=2)
points(max(branching.times(tr))-fit.Ox3.ARD$par[1], fit.Ox3.ARD$lnLik,pch=21, bg="red")
rect(-2, max(lik.prof.ARD-2), 4, max(lik.prof.ARD), col=rgb(0,0,255, maxColorValue=255, alpha=50), border=NA)
rect(max(branching.times(tr))- 2.45, -1000, max(branching.times(tr))-1.85, 0, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(max(branching.times(tr))- 0.663, -1000, max(branching.times(tr))-0.551, 0, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
@



\end{document}