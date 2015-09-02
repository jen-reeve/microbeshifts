## # Evolution of oxygen metabolism in freshwater and marine habitats in the Archaea
## Load packages and prepare environment
require(phytools)
require(diversitree)
require(geiger)
require(devtools)
require(testthat)
require(roxygen2)
require(foreach)
require(doParallel)
require(aRbor)
setwd("~/repos/microbeshifts/analysis")
source("./shiftfunctions.R")


## # Analysis of the Archaea

##Load data
treNSOA <- multi2di(read.tree("../data/NSOA_treeClean.tre"))
datNSOA <- read.table("../data/NSOA_char.txt",row.names=1)
datNSOA[datNSOA=="?"] <- NA
colnames(datNSOA) <- c("Temp", "pH", "Habitat", "Meth", "SulfRed", "SulfOx", "SulfideOx", "SulfThioRed", "NifH", "NitRed", "O2Met")
tdNSOA <- make.treedata(treNSOA, datNSOA)


## Distribution of traits across the Archaeal phylogeny
plot(c(0,12), c(0, 60), type="n", bty="n", xaxt="n", yaxt="n", xlab="",ylab="")
plotTree(tdNSOA$phy, fsize=0.5,add=TRUE)
for(i in 1:ncol(tdNSOA$dat)){
  tiplabels(pch=22, bg=1+as.numeric(as.factor(tdNSOA$dat[,i])),adj=c(5+i/3, 0.5))
  text(5+i/3-2/5+max(nodeHeights(tdNSOA$phy)),length(tdNSOA$phy$tip.l)+4, labels=colnames(datNSOA)[i], srt=75, cex=0.5)
}

#traits <-  c("NitRed", "O2Met", "SulfOx", "SulfideOx")
tdOx <- select(tdNSOA, Habitat, Temp,  NitRed, O2Met, SulfOx, SulfideOx)
tdOx <- filter(tdOx, !(apply(apply(tdNSOA$dat, 1, is.na), 2, any)))

bd.lik <- make.bd(tdOx$phy)
b.lik <- constrain(bd.lik, mu~0)
fit.b <- find.mle(b.lik, 1)
lambda.est <- fit.b$par
mu.est <- 0
bisse.fns <- lapply(1:ncol(tdOx$dat), function(x) make.bisse.t(tdOx$phy, setNames(tdOx$dat[,x], rownames(tdOx$dat)), functions=c(rep("constant.t",4), rep("stepf.t", 2))))

cons.fns <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                    q01.y0~q10.y0, q01.y1~q10.y1, q01.tc ~ q10.tc))
cons.fns.ARD <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                        q01.tc ~ q10.tc))
multi.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                      q01.y0~q01.y1 * r, q10.y0~q10.y1 * r, q01.tc ~ q10.tc, extra="r"))
## Test to make sure multi.fns work:
testpars <- c(3, 2, 3, 2)
bisse.fns[[1]](c(lambda.est, lambda.est, 0, 0, testpars[1]*testpars[2], testpars[2], testpars[4], testpars[1]*testpars[3], testpars[3], testpars[4]))==multi.fns[[1]](testpars)

## Will optimize the following vector pars = c(r, q01.y1, q10.y1, time)
.mlik <- function(rs, times, qs, fxns){
  pars.list <- lapply(fxns, function(x) rep(1, length(argnames(x))))
  pars.list <- lapply(1:length(fxns), function(x) setNames(pars.list[[x]], argnames(fxns[[x]])))
  ## Set R values
  if(length(rs)==1){
    pars.list <- lapply(pars.list, function(x){x['r'] <- rs; x})
  } else {
    pars.list <- lapply(1:length(rs), function(x){pars.list[[x]]['r'] <- rs[x]; pars.list[[x]]})
  }
  ## Set shift times
  if(length(times)==1){
    pars.list <- lapply(pars.list, function(x){x['q10.tc'] <- times; x})
  } else {
    pars.list <- lapply(1:length(times), function(x){pars.list[[x]]['q10.tc'] <- times[x]; pars.list[[x]]})
  }
  start=1
  end=start+1
  for(i in 1:length(pars.list)){
    pars.list[[i]][2:3] <- qs[start:end]
    start=end+1
    end=start+1
  }
  sum(sapply(1:length(fxns), function(x) fxns[[x]](pars.list[[x]])))
}
mlik <- function(pars, fxns=multi.fns){
  rs <- pars[1]
  times <- pars[2]
  qs <- exp(pars[3:length(pars)])
  .mlik(rs, times, qs, fxns)
}
mlik.Rs<- function(pars, fxns=multi.fns){
  rs <- pars[1:4]
  times <- pars[5]
  qs <- exp(pars[6:length(pars)])
  .mlik(rs, times, qs, fxns)
}

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

fit.Ox <- find.mle(all.lik, rep(1,9), method="nlminb", lower=rep(0, 9), upper=c(max(branching.times(tr)), rep(10,8)))
fit.Ox.R <- find.mle(mlik, rep(1,10), method="optim", lower=c(0, 0, rep(-20, 8)), upper=c(10000, max(branching.times(tr)), rep(10,8)))
fit.Ox.Rs <- find.mle(mlik.Rs, rep(1,13), method="optim", lower=c(rep(0,4), 0, rep(-20, 8)), upper=c(rep(10000,4), max(branching.times(tr)), rep(10,8)))
expand.par <- c(fit.Ox.R$par[2], log(unlist(lapply(fit.Ox.R$par[3:10], function(x){y <- c(fit.Ox.R$par[1]*exp(x), exp(x)); y}))))
fit.Ox.ARD <- find.mle(all.lik.ARD, expand.par, method="optim",lower=c(0,rep(-20,8)), upper=c(max(branching.times(tr)), rep(10,8)))

all.lik.ARD(expand.par) - mlik(fit.Ox.R$par) < 1^-13

fit.Ox
plot.likProf(all.lik.ARD, tdOx, fit.Ox.ARD$par, 1, c(0,3.8), -300, n=100, lcol="darkgreen")
plot.likProf(mlik.Rs, tdOx, fit.Ox.Rs$par, 5, c(0,3.8), -300, n=100, add=TRUE, lcol="purple")
plot.likProf(mlik, tdOx, fit.Ox.R$par, 2, c(0,3.8), -300, n=100,  add=TRUE, lcol="blue")
plot.likProf(all.lik, tdOx, fit.Ox$par, 1, c(0,3.8), -300, n=100,  add=TRUE, lcol="red")

## Great Oxygenation Event traits - > Traits = c("Temp","Meth", "SulfOx", "SulfideOx", "NifH", "pH")

tdGO1 <- select(tdNSOA, Temp, Meth, SulfOx, SulfideOx, NifH, pH)
tdGO1 <- filter(tdGO1, !(apply(apply(tdGO1$dat, 1, is.na), 2, any)))
tdGO1$dat[tdGO1$dat==2] <- 1

bd.lik <- make.bd(tdGO1$phy)
b.lik <- constrain(bd.lik, mu~0)
fit.b <- find.mle(b.lik, 1)
lambda.est <- fit.b$par
mu.est <- 0

list.dat <- as.list(tdGO1$dat)
list.dat <- lapply(list.dat, function(x){ names(x) <- rownames(tdGO1$dat); x})
bisse.fns <- lapply(list.dat, function(x) make.bisse.t(tdGO1$phy, x, functions=c(rep("constant.t",4), rep("stepf.t", 2)))) 

cons.fns <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                    q01.y0~q10.y0, q01.y1~q10.y1, q01.tc ~ q10.tc))
cons.fns.ARD <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                        q01.tc ~ q10.tc))
multi.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                      q01.y0~q01.y1 * r, q10.y0~q10.y1 * r, q01.tc ~ q10.tc, extra="r"))

.mlik <- function(rs, times, qs, fxns){
  pars.list <- lapply(fxns, function(x) rep(1, length(argnames(x))))
  pars.list <- lapply(1:length(fxns), function(x) setNames(pars.list[[x]], argnames(fxns[[x]])))
  ## Set R values
  if(length(rs)==1){
    pars.list <- lapply(pars.list, function(x){x['r'] <- rs; x})
  } else {
    pars.list <- lapply(1:length(rs), function(x){pars.list[[x]]['r'] <- rs[x]; pars.list[[x]]})
  }
  ## Set shift times
  if(length(times)==1){
    pars.list <- lapply(pars.list, function(x){x['q10.tc'] <- times; x})
  } else {
    pars.list <- lapply(1:length(times), function(x){pars.list[[x]]['q10.tc'] <- times[x]; pars.list[[x]]})
  }
  start=1
  end=start+1
  for(i in 1:length(pars.list)){
    pars.list[[i]][2:3] <- qs[start:end]
    start=end+1
    end=start+1
  }
  sum(sapply(1:length(fxns), function(x) fxns[[x]](pars.list[[x]])))
}
mlik <- function(pars, fxns=multi.fns){
  rs <- pars[1]
  times <- pars[2]
  qs <- exp(pars[3:length(pars)])
  .mlik(rs, times, qs, fxns)
}
mlik.Rs<- function(pars, fxns=multi.fns){
  rs <- pars[1:6]
  times <- pars[7]
  qs <- exp(pars[8:length(pars)])
  .mlik(rs, times, qs, fxns)
}
.mlikbl <- function(rs, times, qs, fxns){
  pars.list <- lapply(fxns, function(x) rep(1, length(argnames(x))))
  pars.list <- lapply(1:length(fxns), function(x) setNames(pars.list[[x]], argnames(fxns[[x]])))
  ## Set R values
  if(length(rs)==1){
    pars.list <- lapply(pars.list, function(x){x['r'] <- rs; x})
  } else {
    pars.list <- lapply(1:length(rs), function(x){pars.list[[x]]['r'] <- rs[x]; pars.list[[x]]})
  }
  ## Set shift times
  if(length(times)==1){
    pars.list <- lapply(pars.list, function(x){x['q10.tc'] <- times; x})
  } else {
    pars.list <- lapply(1:length(times), function(x){pars.list[[x]]['q10.tc'] <- times[x]; pars.list[[x]]})
  }
  start=1
  end=start+1
  for(i in 1:length(pars.list)){
    pars.list[[i]][2:3] <- qs[start:end]
    start=end+1
    end=start+1
  }
  sum(sapply(1:length(fxns), function(x) fxns[[x]](pars.list[[x]])))
}
mlik.Rblock <- function(pars, fxns=multi.fns){
  rs <- pars[1:2]
  times <- pars[3]
  qs <- exp(pars[4:length(pars)])
  .mlikbl(rs, times, qs, fxns)
}



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

fit.mlik <- find.mle(mlik, rep(1,14), method="optim", lower=c(0,0,rep(-20,12)), upper=c(10000, max(branching.times(tdGO1$phy)), rep(10,12)))
fit.mlik.Rs <- find.mle(mlik.Rs, rep(1,19), method="optim", lower=c(rep(0,6),0,rep(-20,12)), upper=c(rep(10000,6), max(branching.times(tdGO1$phy)), rep(10,12)))

plot.likProf(mlik, tdGO1, fit.mlik$par, 2, c(3.8,0), c(-450, -370), n=100, lcol="darkgreen")
plot.likProf(mlik.Rs, tdGO1, fit.mlik.Rs$par, 7, c(3.8,0), c(-450, -370), n=100, lcol="blue", add=TRUE)
plot.likProf(mlik, tdGO1, c(1, fit.mlik$par[-1]), 2, c(3.8,0), c(-450, -370), n=100, lcol="red", add=TRUE)



## ALL Traits = c("Temp","Meth", "SulfOx", "SulfideOx", "NifH", "pH")

tdALL <- tdNSOA
tdALL$dat[tdALL$dat==2] <- 1
tdALL$dat[tdALL$dat=="0&1"] <- 1
tdALL$dat[tdALL$dat==3] <- 1
tdALL$dat$Habitat <- as.numeric(as.character(tdALL$dat$Habitat))
tdALL <- filter(tdALL, !(apply(apply(tdALL$dat, 1, is.na), 2, any)))
unique(unlist(tdALL$dat))

bd.lik <- make.bd(tdALL$phy)
b.lik <- constrain(bd.lik, mu~0)
fit.b <- find.mle(b.lik, 1)
lambda.est <- fit.b$par
mu.est <- 0

list.dat <- as.list(tdALL$dat)
list.dat <- lapply(list.dat, function(x){ names(x) <- rownames(tdALL$dat); x})
bisse.fns <- lapply(list.dat, function(x) make.bisse.t(tdALL$phy, x, functions=c(rep("constant.t",4), rep("stepf.t", 2)))) 
notime.fns <- lapply(list.dat, function(x) make.bisse(tdALL$phy, x))
cons.fns <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                    q01.y0~q10.y0, q01.y1~q10.y1, q01.tc ~ q10.tc))
cons.fns.ARD <- lapply(bisse.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                        q01.tc ~ q10.tc))
multi.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, 
                                                      q01.y0~q01.y1 * r, q10.y0~q10.y1 * r, q01.tc ~ q10.tc, extra="r"))
notime.fns <- lapply(notime.fns, function(x) constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0))

.mlik <- function(rs, times, qs, fxns){
  pars.list <- lapply(fxns, function(x) rep(1, length(argnames(x))))
  pars.list <- lapply(1:length(fxns), function(x) setNames(pars.list[[x]], argnames(fxns[[x]])))
  ## Set R values
  if(length(rs)==1){
    pars.list <- lapply(pars.list, function(x){x['r'] <- rs; x})
  } else {
    pars.list <- lapply(1:length(rs), function(x){pars.list[[x]]['r'] <- rs[x]; pars.list[[x]]})
  }
  ## Set shift times
  if(length(times)==1){
    pars.list <- lapply(pars.list, function(x){x['q10.tc'] <- times; x})
  } else {
    pars.list <- lapply(1:length(times), function(x){pars.list[[x]]['q10.tc'] <- times[x]; pars.list[[x]]})
  }
  start=1
  end=start+1
  for(i in 1:length(pars.list)){
    pars.list[[i]][2:3] <- qs[start:end]
    start=end+1
    end=start+1
  }
  sum(sapply(1:length(fxns), function(x) fxns[[x]](pars.list[[x]])))
}
mlik <- function(pars, fxns=multi.fns){
  rs <- pars[1]
  times <- pars[2]
  qs <- exp(pars[3:length(pars)])
  .mlik(rs, times, qs, fxns)
}
mlik.Rs<- function(pars, fxns=multi.fns){
  rs <- pars[1:11]
  times <- pars[12]
  qs <- exp(pars[13:length(pars)])
  .mlik(rs, times, qs, fxns)
}
.mlikbl <- function(rs, times, qs, fxns){
  pars.list <- lapply(fxns, function(x) rep(1, length(argnames(x))))
  pars.list <- lapply(1:length(fxns), function(x) setNames(pars.list[[x]], argnames(fxns[[x]])))
  ## Set R values
  if(length(rs)==1){
    pars.list <- lapply(pars.list, function(x){x['r'] <- rs; x})
  } else {
    pars.list <- lapply(1:length(rs), function(x){pars.list[[x]]['r'] <- rs[x]; pars.list[[x]]})
  }
  ## Set shift times
  if(length(times)==1){
    pars.list <- lapply(pars.list, function(x){x['q10.tc'] <- times; x})
  } else {
    pars.list <- lapply(1:length(times), function(x){pars.list[[x]]['q10.tc'] <- times[x]; pars.list[[x]]})
  }
  start=1
  end=start+1
  for(i in 1:length(pars.list)){
    pars.list[[i]][2:3] <- qs[start:end]
    start=end+1
    end=start+1
  }
  sum(sapply(1:length(fxns), function(x) fxns[[x]](pars.list[[x]])))
}
mlik.Rblock <- function(pars, fxns=multi.fns){
  rs <- pars[1:2]
  times <- pars[3]
  qs <- exp(pars[4:length(pars)])
  .mlikbl(rs, times, qs, fxns)
}

nlik <- function(pars, fxns=notime.fns){
  ind <- seq(1, length(pars), 2)
  pars <- exp(pars)
  liks <- sapply(1:length(ind), function(x) notime.fns[[x]](c(pars[ind[x]:(ind[x]+1)])))
  sum(liks)
}

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

registerDoMC(cores=8)
fit.nlik <- foreach(i=1:5) %dopar% find.mle(nlik, runif(22, -10, 10), method="optim", lower=rep(-20, 22), upper=rep(10, 22))
nliks <- sapply(fit.nlik, function(x) x$lnLik)
fit.mlik <- foreach(i=1:5) %dopar% find.mle(mlik, c(exp(runif(1,-1,1)),runif(1, 0.5, 3.7), fit.nlik[[i]]$par), method="optim", lower=c(0,0,rep(-20,22)), upper=c(10000, max(branching.times(tdALL$phy)), rep(10,22)))
mliks <- sapply(fit.mlik, function(x) x$lnLik)
fit.mlik.Rs <- foreach(i=1:5) %dopar% find.mle(mlik.Rs, c(exp(runif(11,-5,5)),runif(1, 0.5, 3.7), fit.nlik[[i]]$par), method="optim", lower=c(rep(0,11),0,rep(-20,22)), upper=c(rep(10000,11), max(branching.times(tdALL$phy)), rep(10,22)))
mlikRs <- sapply(fit.mlik.Rs, function(x) x$lnLik)

cat("AIC for no shift model\n",
    2*length(fit.nlik[[1]]$par) - 2*max(nliks),
    "\nAIC for shift model \n",
    2*length(fit.mlik[[1]]$par) - 2*max(mliks), 
    "\nAIC for the shift separate R model\n",
    2*length(fit.mlik.Rs[[1]]$par) - 2*max(mlikRs))

parests <- cbind(exp(fit.nlik[[5]]$par)*TH, exp(fit.mlik[[5]]$par[3:24])*TH, TH*fit.mlik[[5]]$par[1]*exp(fit.mlik[[5]]$par[3:24]),fit.mlik[[5]]$par[1],
                 TH* exp(fit.mlik.Rs[[4]]$par[13:34]),TH*as.vector(sapply(1:11, function(x) rep(fit.mlik.Rs[[4]]$par[x], 2)))* exp(fit.mlik.Rs[[4]]$par[13:34]),as.vector(sapply(1:11, function(x) rep(fit.mlik.Rs[[4]]$par[x], 2))))
rownames(parests) <- unlist(sapply(1:11, function(x) c(paste(colnames(tdALL$dat[x]), "01"), paste(colnames(tdALL$dat[x]), "10"))))
parests
ymin = -716
ymax= -713
plot.likProf(mlik, tdALL, fit.mlik[[which(mliks==max(mliks))]]$par, 2, c(3.8,0), c(ymin, ymax), n=100, lcol="darkgreen")
plot.likProf(mlik.Rs, tdALL, fit.mlik.Rs[[which(mlikRs==max(mlikRs))]]$par, 12, c(3.8,0), c(ymin, ymax), n=100, lcol="blue", add=TRUE)
abline(h=fit.nlik[[which(nliks==max(nliks))]]$lnLik, col="red")


## # Individual fits 0&1 -> 0; >=2 -> 1
tdList <- lapply(1:11, function(x) select(tdNSOA, x))
tdList <- lapply(1:11, function(x) filter(tdList[[x]], !is.na(tdList[[x]]$dat[,1])))
tdList <- lapply(1:11, function(x){tdList[[x]]$dat[,1] <- binarify(tdList[[x]]$dat[,1]); tdList[[x]]})

## Set global birth-death parameters
bd.lik <- make.bd(tdNSOA$phy)
bd.est <- find.mle(bd.lik, x.init=c(1,1))
lambda <- bd.est$par[1]
mu <- bd.est$par[2]

## Make individual bisse functions for each using multipliers
bisse.fns <- lapply(1:11, function(x) make.bisse.t(tdList[[x]]$phy, setNames(tdList[[x]]$dat[,1], rownames(tdList[[x]]$dat)), functions=c(rep("constant.t",4), rep("stepf.t", 2)))) 
notime.fns <- lapply(1:11, function(x) make.bisse(tdList[[x]]$phy, setNames(tdList[[x]]$dat[,1], rownames(tdList[[x]]$dat))))
ARD.notime.fns <- lapply(notime.fns, function(x) constrain(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu))
mk.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu, q01.y1~r1*q01.y0, q10.y1~r2*q10.y0, q01.tc ~ q10.tc, extra=c('r1', 'r2')))
ER.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu, q10.y1~r*q10.y0, q01.y0~q10.y0, q01.y1~r*q10.y0, q01.tc ~ q10.tc, extra=c('r')))
ARD.R1.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu, q10.y1~r*q10.y0, q01.y1~r*q10.y0, q01.tc ~ q10.tc, extra=c('r')))
ARD.R2.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu, q10.y1~r1*q10.y0, q01.y1~r2*q01.y0, q01.tc ~ q10.tc, extra=c('r1', 'r2')))


## Begin by fitting the full model to each one. This gives the maximal likelihood I should be expecting.
registerDoMC(cores=11)
for(i in 1:20){
  Ind.Fits <- foreach(i=1:11) %dopar% {find.mle(ARD.R2.fns[[i]], x.init=start.gen(ARD.R2.fns[[i]], tdList[[i]]), method="optim", lower=lower.gen(ARD.R2.fns[[i]], tdList[[i]]), upper=upper.gen(ARD.R2.fns[[i]], tdList[[i]]) )}
  ## Create a summary table
  parests <- do.call(rbind, lapply(Ind.Fits, function(x) as.data.frame(matrix(c(x$par, x$lnLik, x$message), nrow=1))))
  colnames(parests) <- c(argnames(ARD.R2.fns[[1]]), "lnL", "message")
  ## Save only the best-fitting indpendent runs
  if(i==1){
    IndFit.best <- parests
  } else {
    replace <- which(defactor(parests$lnL) > defactor(IndFit.best$lnL))
    if(length(replace) > 0){
      IndFit.best[-ncol(IndFit.best)] <- apply(IndFit.best[-ncol(IndFit.best)], 2, defactor)
      parests[-ncol(parests)] <- apply(parests[-ncol(parests)], 2, defactor)
      IndFit.best[replace, ] <- parests[replace, ]
    }
  }
}


## Compare to the no time-shift model. See what the maximum likelihood difference we could obtain is.
for(i in 1:20){
  Ind.Fits.notime <- foreach(i=1:11) %dopar% {find.mle(ARD.notime.fns[[i]], x.init=start.gen(ARD.notime.fns[[i]], tdList[[i]]), method="optim", lower=lower.gen(ARD.notime.fns[[i]], tdList[[i]]), upper=upper.gen(ARD.notime.fns[[i]], tdList[[i]]) )}
  ## Create a summary table
  parests.notime <- do.call(rbind, lapply(Ind.Fits.notime, function(x) as.data.frame(matrix(c(x$par, x$lnLik, x$message), nrow=1))))
  colnames(parests.notime) <- c(argnames(ARD.notime.fns[[1]]), "lnL", "message")
  ## Save only the best-fitting indpendent runs
  if(i==1){
    IndFit.best.notime <- parests.notime
  } else {
    replace <- which(defactor(parests.notime$lnL) > defactor(IndFit.best.notime$lnL))
    if(length(replace) > 0){
      IndFit.best.notime[-ncol(IndFit.best.notime)] <- apply(IndFit.best.notime[-ncol(IndFit.best.notime)], 2, defactor)
      parests.notime[-ncol(parests.notime)] <- apply(parests.notime[-ncol(parests.notime)], 2, defactor)
      IndFit.best.notime[replace, ] <- parests.notime[replace, ]
    }
  }
}

## Maximum amount the full model could be better than the no-time shift model
sum(defactor(IndFit.best$lnL)) - sum(defactor(IndFit.best.notime$lnL))

## Compare to the GO1 constrained model. See what the maximum likelihood difference we could obtain is.
for(i in 1:50){
  Ind.Fits.GO1 <- foreach(i=1:11) %dopar% {find.mle(ARD.R2.fns[[i]], x.init=start.genGO1(ARD.R2.fns[[i]], tdList[[i]]), method="optim", lower=c(0, 0, 0, 0, 1.85), upper=c(10000, 10000, 10000, 10000, 2.45))}
  ## Create a summary table
  parests.GO1 <- do.call(rbind, lapply(Ind.Fits.GO1, function(x) as.data.frame(matrix(c(x$par, x$lnLik, x$message), nrow=1))))
  colnames(parests.GO1) <- c(argnames(ARD.R2.fns[[1]]), "lnL", "message")
  ## Save only the best-fitting indpendent runs
  if(i==1){
    #IndFit.best.GO1 <- parests.GO1
  } else {
    replace <- which(defactor(parests.GO1$lnL) > defactor(IndFit.best.GO1$lnL))
    if(length(replace) > 0){
      IndFit.best.GO1[-ncol(IndFit.best.GO1)] <- apply(IndFit.best.GO1[-ncol(IndFit.best.GO1)], 2, defactor)
      parests.GO1[-ncol(parests.GO1)] <- apply(parests.GO1[-ncol(parests.GO1)], 2, defactor)
      IndFit.best.GO1[replace, ] <- parests.GO1[replace, ]
    }
    print(replace)
  }
}

## Compare to the GO2 constrained model. See what the maximum likelihood difference we could obtain is.
for(i in 1:50){
  Ind.Fits.GO2 <- foreach(i=1:11) %dopar% {find.mle(ARD.R2.fns[[i]], x.init=start.genGO2(ARD.R2.fns[[i]], tdList[[i]]), method="optim", lower=c(0, 0, 0, 0, 0.551), upper=c(10000, 10000, 10000, 10000, 0.663))}
  ## Create a summary table
  parests.GO2 <- do.call(rbind, lapply(Ind.Fits.GO2, function(x) as.data.frame(matrix(c(x$par, x$lnLik, x$message), nrow=1))))
  colnames(parests.GO2) <- c(argnames(ARD.R2.fns[[1]]), "lnL", "message")
  ## Save only the best-fitting indpendent runs
  if(i==1){
    #IndFit.best.GO2 <- parests.GO2
  } else {
    replace <- which(defactor(parests.GO2$lnL) > defactor(IndFit.best.GO2$lnL))
    if(length(replace) > 0){
      IndFit.best.GO2[-ncol(IndFit.best.GO2)] <- apply(IndFit.best.GO2[-ncol(IndFit.best.GO2)], 2, defactor)
      parests.GO2[-ncol(parests.GO2)] <- apply(parests.GO2[-ncol(parests.GO2)], 2, defactor)
      IndFit.best.GO2[replace, ] <- parests.GO2[replace, ]
    }
    print(replace)
  }
}

## Examine the independent time-shift models
summary.GOevents <- cbind(IndFit.best.notime$lnL, IndFit.best$lnL, IndFit.best.GO1$lnL, IndFit.best.GO2$lnL,round(IndFit.best.GO1$lnL-IndFit.best.notime$lnL,2),round(IndFit.best.GO2$lnL-IndFit.best.notime$lnL,2), round(IndFit.best$lnL-IndFit.best.notime$lnL,2), IndFit.best$q10.tc)
rownames(IndFit.best.notime) <- rownames(IndFit.best.GO1) <- rownames(IndFit.best.GO2) <- rownames(summary.GOevents) <- colnames(tdNSOA$dat)
colnames(summary.GOevents) <- c("notime", "full", "GO1", "GO2", "GO1-notime", "GO2-notime", "full-notime", "time.est")
summary.GOevents
apply(summary.GOevents[,5:7], 2, sum)
IndFit.best.GO2[!(rownames(IndFit.best.GO2) %in% c("Habitat", "Meth", "SulfThioRed", "NifH")),]
sum(apply(cbind(IndFit.best.GO1$lnL, IndFit.best.GO2$lnL), 1, max))

## Individual profile plots:
par(mfrow=c(3,3), ask=TRUE)
lapply(1:11, function(x) plot.likProf(ARD.R2.fns[[x]], tdList[[1]], unlist(IndFit.best[x, 1:5]), 5, c(TH, 0), ylim=NULL))
