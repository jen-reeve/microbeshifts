## # Evolution of oxygen metabolism in freshwater and marine habitats in the Archaea
## Load packages and prepare environment
require(phytools)
require(diversitree)
require(geiger)
require(devtools)
require(testthat)
require(roxygen2)
require(foreach)
require(doMC)
load_all("~/repos/aRbor/aRbor")
setwd("~/repos/microbeshifts/analysis")
source("./shiftfunctions.R")


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



## # Same Analysis for Cyanobacteria
## Data preperation stuff
cyanodat <- read.table("../data/charactermatrix.txt")
celldiam <- read.table("../data/celldiameter.msq")
rownames(celldiam) <- celldiam[,1]
celldiam <- celldiam[,-1]
rownames(cyanodat) <- cyanodat[,1]
cyanodat <- cyanodat[,-1]
states <- readLines("../data/charactersandstates.txt")
states[grep(".", states)]
heads <- strsplit(states[grep(".", states)], "\\. ", perl=TRUE)
colnames(cyanodat) <- sapply(gsub("-", "", gsub(" ", "_", as.data.frame(do.call(rbind, heads[sapply(heads, length)==2]))[,2])), function(x) substr(x, 1, nchar(x)-2))
cyanodat$celldiam_min <- celldiam[,1]
cyanodat$celldiam_mean <- celldiam[,2]
cyanodat$celldiam_max <- celldiam[,3]
cyanodat[cyanodat=="?"] <- NA

## # Read in the phylogenies and match to names
#lapply(td1$dat,function(x) table(factor(x)))
whichcol <- c('Thermophilic', 'Nonfreshwater_habitat', 'Akinetes', 'Heterocysts', 'Nitrogen_fixation', 'Morphology',
              'Habit', 'Freeliving', 'Mats', 'Epi/Endolithic',  'Epiphytic', 'Periphytic', 'Motility', 'Hormogonia', 
              'Gas_vesicles', 'False_Branching', 'True_Branching', 'Fission_in_multiple_planes', 'Uniseriate_trichome', 
              'Multiseriate_trichomes', 'Baeocytes', 'Extracellular_sheath', 'Mucilage', 'celldiam_mean')
dat <- cyanodat
dat$Morphology[cyanodat$Morphology=="0&2"] <- 0
dat$Motility[cyanodat$Motility=="2"] <- 1
dat$Multiseriate_trichomes[cyanodat$Multiseriate_trichomes!=0] <- 1
dat$Mucilage[cyanodat$Mucilage!=0] <- 1
dat$Habit[cyanodat$Habit=="0&1"] <- 0
dat <- dat[,whichcol]
dat$Pelagic <- as.numeric(dat$Nonfreshwater_habitat==1 & dat$Habit==0)
dat$celldiam_mean <- as.numeric(as.numeric(as.character(cyanodat$celldiam_mean))>=3.5)

tree1 <- read.tree("../data/Tree1MrB.tre")
tree1$edge.length <- tree1$edge.length/1000

tdcy <- make.treedata(tree1, dat, name_column=0)
nc <- ncol(tdcy$dat)
tdcyList <- lapply(1:nc, function(x) select(tdcy, x))
tdcyList <- lapply(1:nc, function(x) filter(tdcyList[[x]], !is.na(tdcyList[[x]]$dat[,1])))
#tdcyList[[23]] <- mutate(tdcyList[[23]], logCD = log(celldiam_mean))
#tdcyList[[23]] <- select(tdcyList[[23]], 2)
#tdcyList <- lapply(1:nc, function(x){tdcyList[[x]]$dat[,1] <- binarify(tdcyList[[x]]$dat[,1]); tdcyList[[x]]})

## Check to see that all characters are binary (except for cell diameter)
all(sapply(tdcy$dat[,-ncol(tdcy$dat)], function(x) levels(factor(x))==c("0", "1")))

## Visualize the distribution of traits:
asrDiscrete <- aceArbor(tdcy, charType="discrete", aceType="marginal", discreteModelType="ARD", na.rm="bytrait")
par(mfrow=c(2,2))
plot(asrDiscrete, type="fan", show.tip.label=FALSE)

asrCont <- aceArbor(tdcyList[[23]], charType="continuous")
plot(asrCont, cex=1, show.tip.label=FALSE, adj=c(0.5, 8))


## Set global birth-death parameters
bd.lik <- make.bd(tdcy$phy)
bd.est <- find.mle(bd.lik, x.init=c(5,0))
lambda <- bd.est$par[1]
mu <- bd.est$par[2]

## Make the functions
bisse.fns <- lapply(1:nc, function(x) make.bisse.t(tdcyList[[x]]$phy, setNames(tdcyList[[x]]$dat[,1], rownames(tdcyList[[x]]$dat)), functions=c(rep("constant.t",4), rep("stepf.t", 2)))) 
notime.fns <- lapply(1:nc, function(x) make.bisse(tdcyList[[x]]$phy, setNames(tdcyList[[x]]$dat[,1], rownames(tdcyList[[x]]$dat))))
ARD.notime.fns <- lapply(notime.fns, function(x) constrain(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu))
#mk.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu, q01.y1~r1*q01.y0, q10.y1~r2*q10.y0, q01.tc ~ q10.tc, extra=c('r1', 'r2')))
#ER.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu, q10.y1~r*q10.y0, q01.y0~q10.y0, q01.y1~r*q10.y0, q01.tc ~ q10.tc, extra=c('r')))
ARD.R1.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu, q10.y1~r*q10.y0, q01.y1~r*q10.y0, q01.tc ~ q10.tc, extra=c('r')))
ARD.R2.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu, q10.y1~r1*q10.y0, q01.y1~r2*q01.y0, q01.tc ~ q10.tc, extra=c('r1', 'r2')))

## Begin by fitting the full model to each one. This gives the maximal likelihood I should be expecting.
registerDoMC(cores=11)
for(i in 1:20){
  cyFits <- foreach(i=1:22) %dopar% {find.mle(ARD.R2.fns[[i]], x.init=start.gen(ARD.R2.fns[[i]], tdcyList[[i]]), method="optim", lower=lower.gen(ARD.R2.fns[[i]], tdcyList[[i]]), upper=upper.gen(ARD.R2.fns[[i]], tdcyList[[i]]) )}
  ## Create a summary table
  parests <- do.call(rbind, lapply(cyFits, function(x) as.data.frame(matrix(c(x$par, x$lnLik, x$message), nrow=1))))
  colnames(parests) <- c(argnames(ARD.R2.fns[[1]]), "lnL", "message")
  ## Save only the best-fitting indpendent runs
  if(i==Inf){
    #cyFit.best <- parests
  } else {
    replace <- which(defactor(parests$lnL) > defactor(cyFit.best$lnL))
    if(length(replace) > 0){
      cyFit.best[-ncol(cyFit.best)] <- apply(cyFit.best[-ncol(cyFit.best)], 2, defactor)
      parests[-ncol(parests)] <- apply(parests[-ncol(parests)], 2, defactor)
      cyFit.best[replace, ] <- parests[replace, ]
    }
    print(replace)
  }
}

## Compare to the no time-shift model. See what the maximum likelihood difference we could obtain is.
for(i in 1:20){
  cyFits.notime <- foreach(i=1:22) %dopar% {find.mle(ARD.notime.fns[[i]], x.init=start.gen(ARD.notime.fns[[i]], tdcyList[[i]]), method="optim", lower=lower.gen(ARD.notime.fns[[i]], tdcyList[[i]]), upper=upper.gen(ARD.notime.fns[[i]], tdcyList[[i]]) )}
  ## Create a summary table
  parests.notime <- do.call(rbind, lapply(cyFits.notime, function(x) as.data.frame(matrix(c(x$par, x$lnLik, x$message), nrow=1))))
  colnames(parests.notime) <- c(argnames(ARD.notime.fns[[1]]), "lnL", "message")
  ## Save only the best-fitting indpendent runs
  if(i==Inf){
    cyFit.best.notime <- parests.notime
  } else {
    replace <- which(defactor(parests.notime$lnL) > defactor(cyFit.best.notime$lnL))
    if(length(replace) > 0){
      cyFit.best.notime[-ncol(cyFit.best.notime)] <- apply(cyFit.best.notime[-ncol(cyFit.best.notime)], 2, defactor)
      parests.notime[-ncol(parests.notime)] <- apply(parests.notime[-ncol(parests.notime)], 2, defactor)
      cyFit.best.notime[replace, ] <- parests.notime[replace, ]
    }
    print(replace)
  }
}

## Compare to the GO1 constrained model. See what the maximum likelihood difference we could obtain is.
for(i in 1:10){
  cy.Fits.GO1 <- foreach(i=1:22) %dopar% {find.mle(ARD.R2.fns[[i]], x.init=start.genGO1(ARD.R2.fns[[i]], tdcyList[[i]]), method="optim", lower=c(0, 0, 0, 0, 1.85), upper=c(10000, 10000, 10000, 10000, 2.45))}
  ## Create a summary table
  parests.GO1 <- do.call(rbind, lapply(cy.Fits.GO1, function(x) as.data.frame(matrix(c(x$par, x$lnLik, x$message), nrow=1))))
  colnames(parests.GO1) <- c(argnames(ARD.R2.fns[[1]]), "lnL", "message")
  ## Save only the best-fitting indpendent runs
  if(i==1){
    cyFit.best.GO1 <- parests.GO1
  } else {
    replace <- which(defactor(parests.GO1$lnL) > defactor(cyFit.best.GO1$lnL))
    if(length(replace) > 0){
      cyFit.best.GO1[-ncol(cyFit.best.GO1)] <- apply(cyFit.best.GO1[-ncol(cyFit.best.GO1)], 2, defactor)
      parests.GO1[-ncol(parests.GO1)] <- apply(parests.GO1[-ncol(parests.GO1)], 2, defactor)
      cyFit.best.GO1[replace, ] <- parests.GO1[replace, ]
    }
    print(replace)
  }
}

## Compare to the GO2 constrained model. See what the maximum likelihood difference we could obtain is.
for(i in 1:10){
  cy.Fits.GO2 <- foreach(i=1:22) %dopar% {find.mle(ARD.R2.fns[[i]], x.init=start.genGO2(ARD.R2.fns[[i]], tdcyList[[i]]), method="optim", lower=c(0, 0, 0, 0, 0.551), upper=c(10000, 10000, 10000, 10000, 0.663))}
  ## Create a summary table
  parests.GO2 <- do.call(rbind, lapply(cy.Fits.GO2, function(x) as.data.frame(matrix(c(x$par, x$lnLik, x$message), nrow=1))))
  colnames(parests.GO2) <- c(argnames(ARD.R2.fns[[1]]), "lnL", "message")
  ## Save only the best-fitting indpendent runs
  if(i==1){
    cyFit.best.GO2 <- parests.GO2
  } else {
    replace <- which(defactor(parests.GO2$lnL) > defactor(cyFit.best.GO2$lnL))
    if(length(replace) > 0){
      cyFit.best.GO2[-ncol(cyFit.best.GO2)] <- apply(cyFit.best.GO2[-ncol(cyFit.best.GO2)], 2, defactor)
      parests.GO2[-ncol(parests.GO2)] <- apply(parests.GO2[-ncol(parests.GO2)], 2, defactor)
      cyFit.best.GO2[replace, ] <- parests.GO2[replace, ]
    }
    print(replace)
  }
}

## With only 1 R
for(i in 1:20){
  cy.Fits.GO2.R1 <- foreach(i=1:22) %dopar% {find.mle(ARD.R1.fns[[i]], x.init=start.genGO2(ARD.R1.fns[[i]], tdcyList[[i]]), method="optim", lower=c(0, 0, 0, 0.551), upper=c(10000, 10000, 10000, 0.663))}
  ## Create a summary table
  parests.GO2.R1 <- do.call(rbind, lapply(cy.Fits.GO2.R1, function(x) as.data.frame(matrix(c(x$par, x$lnLik, x$message), nrow=1))))
  colnames(parests.GO2.R1) <- c(argnames(ARD.R1.fns[[1]]), "lnL", "message")
  ## Save only the best-fitting indpendent runs
  if(i==1){
    cyFit.best.GO2.R1 <- parests.GO2.R1
  } else {
    replace <- which(defactor(parests.GO2.R1$lnL) > defactor(cyFit.best.GO2.R1$lnL))
    if(length(replace) > 0){
      cyFit.best.GO2.R1[-ncol(cyFit.best.GO2.R1)] <- apply(cyFit.best.GO2.R1[-ncol(cyFit.best.GO2.R1)], 2, defactor)
      parests.GO2.R1[-ncol(parests.GO2.R1)] <- apply(parests.GO2.R1[-ncol(parests.GO2.R1)], 2, defactor)
      cyFit.best.GO2.R1[replace, ] <- parests.GO2.R1[replace, ]
    }
    print(replace)
  }
}

cySummary <- cbind(defactor(cyFit.best$lnL), defactor(cyFit.best.notime$lnL), defactor(cyFit.best.GO1$lnL),defactor(cyFit.best.GO2$lnL), round(defactor(cyFit.best$lnL)- defactor(cyFit.best.notime$lnL),2),round(defactor(cyFit.best.GO1$lnL)- defactor(cyFit.best.notime$lnL),2),round(defactor(cyFit.best.GO2$lnL)- defactor(cyFit.best.notime$lnL),2), cyFit.best$q10.tc)
rownames(cySummary) <- colnames(tdcy$dat)[1:22]
cySummary
sum(cySummary[,3])
par(mfrow=c(1,2))

plot(cySummary[,4], cySummary[,3], xlim=c(0,3.8))
abline(v= c(0.663, 2.25))

cutoff=0
plot(density(cySummary[cySummary[,3]>cutoff,4],adj=0.1), xlim=c(3.8,0))
rect(2.45, -1000, 1.85, 1000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(0.663, -1000, 0.551, 1000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)

cySummary[cySummary[,3] > cutoff,]
## Individual profile plots:
par(mfrow=c(3,3), ask=TRUE)
lapply(1:22, function(x) plot.likProf(ARD.R2.fns[[x]], tdcyList[[1]], unlist(cyFit.best[x, 1:5]), 5, c(TH, 0), ylim=NULL,main=colnames(tdcy$dat)[x]) )

## Fit model with a single, estimated time shift and a single rate multiplier
## Define joint likelihood function
#Parameter order c(t, r, q01.1, q1.10.1, ..., q01.22, q10.22); total of 46 parameters
TH <- max(branching.times(tdcy$phy))
pars <- c(2.5, rep(0, 51))
lik.ARD.R1.T1 <- function(pars, fns= ARD.R1.fns){
  time <- pars[1]
  r <- exp(pars[2])
  qs <- exp(matrix(pars[3:length(pars)], ncol=2, byrow=TRUE))
  pars <- cbind(r, qs, time)
  lnLs <- sapply(1:length(fns), function(x) fns[[x]](pars[x,]))
  res <- sum(lnLs)
  return(res)
}
## Set function attributes
attributes(lik.ARD.R1.T1)$varnames <- c("t", "r", paste(c("q01", "q10"), sort(rep(1:25,2)), sep="."))
attributes(lik.ARD.R1.T1)$converter <- function(pars){
  pars[2:length(pars)] <- exp(pars[2:length(pars)])
  names(pars) <- attributes(lik.ARD.R1.T1)$varnames
  pars
}
attributes(lik.ARD.R1.T1)$starter <- function(){
  c(runif(1,0, TH), rnorm(53, 0, 2))
}
attributes(lik.ARD.R1.T1)$smartstart <- function(heat=1){
  mean <- log(c(1, as.vector(t(cyFit.best.notime[,1:2]))))
  time <- 0.5*TH
  pars <- c(time, rnorm(45, mean, heat))
}
attributes(lik.ARD.R1.T1)$lower <- c(0, rep(-20, 45))
attributes(lik.ARD.R1.T1)$upper <- c(TH, rep(10, 45))

## Optimize joint likelihood function
fn <- lik.ARD.R1.T1
registerDoMC(cores=3)
res2 <- foreach(i=1:5) %dopar% find.mle(fn, x.init=attributes(fn)$smartstart(0.5), method="subplex", lower=attributes(fn)$lower, upper=attributes(fn)$upper)
res.liks <- sapply(res, function(x) x$lnLik)
best.id <- which(res.liks == max(res.liks))
best.res <- res[[best.id]]

attributes(fn)$converter(best.res$par)
AIC.notime <- (2*44 - 2*sum(cySummary[,2]))
AIC.ARD.R1.T1 <- (2*46 - 2*best.res$lnLik)
AIC.notime
AIC.ARD.R1.T1
AIC.notime-AIC.ARD.R1.T1

plot.likProf(fn, tdcy, best.res$par, 1, c(TH, 0))
abline(v=sapply(res, function(x) x$par[1]), col="red", lty=2)
points(best.res$par[1], best.res$lnLik, pch=21, bg="red", cex=1.25)
abline(h=sum(cySummary[,2]), col="blue")

colnames(cySummary) <- c("full", "notime", "GO1", "GO2","full-notime", "GO1-notime", "GO2-notime", "time.est")
cbind(sapply(1:22, function(x) ARD.R1.fns[[x]](c(exp(best.res$par[2]), exp(best.res$par[(2*x+1):(2*x+2)]), best.res$par[1]))), cySummary)
image(t(as.matrix(round(log(cyFit.best[,1:2]+0.01),3))))
plot(c(0, 4), c(-10, 10), type="n", xlab="time")
points(cyFit.best$q10.tc, log(cyFit.best$r1), col='red')
points(cyFit.best$q10.tc, log(cyFit.best$r2), pch=21, bg="red")

plot(log(cyFit.best$r1), log(cyFit.best$r2), pch=21, bg="green", xlim=c(-10, 10), ylim=c(-10, 10), cex=cySummary[,5]+.5, col=factor(rownames(cySummary)),lwd=3)
curve(1*x, add=TRUE)
abline(v=0)
abline(h=0)
points(log(cyFit.best.GO1$r1), log(cyFit.best.GO1$r2), pch=21, bg="blue", col=factor(rownames(cySummary)), cex=cySummary[,6]+.5,lwd=3)
points(log(cyFit.best.GO2$r1), log(cyFit.best.GO2$r2), pch=21, bg="red", col=factor(rownames(cySummary)), cex=cySummary[,7]+.5,lwd=3)
points(best.res$par[2], best.res$par[2], pch=23, bg="pink", cex=2)


## Get likelihood profiles for across entire region for each trait.
seq1 <- seq(0.2, 3.6, 0.2)
seq2 <- seq(0.1, 3.7, 0.2)

## Start only over seq1, add seq2 if it looks productive
registerDoMC(cores=8)
fns <- ARD.R2.fns
profres <- foreach(i=1:22) %dopar% {
  tmp <- list()
  for(j in 1:length(seq1)) {
    fn <- fns[[i]]
    ft <<- seq1[j]
    fn <- constrain2(fn, q10.tc~ft, extra=c("r1", "r2"))
    if(j==1){
      startx <- runif(4, 0.5, 2)*c(1, 1, cyFit.best.notime[i, 1], cyFit.best.notime[i, 2])
    } else {
      startx <- runif(4, 0.5, 2)*c(1, 1, cyFit.best.notime[i, 1], cyFit.best.notime[i, 2]) #c(tmp[[j-1]]$par)
    }
    tmp[[j]] <- find.mle(fn, x.init=startx, method="optim")
  }
  tmp
}

save(profres, file="profres.rds")
load(file="profres.rds")

profliks <- sapply(1:22, function(x) sapply(profres[[x]], function(y) y$lnLik))
profliksN <- apply(profliks, 2, function(x) x-min(x))
plot(c(3.8, 0), c(0, 6), type="n")
for(i in 1:22){
  lines(seq1, profliksN[,i], col=i)
}

plot(seq1, apply(profliks, 1, sum))
lines(seq1,  apply(profliks, 1, sum))
rect(2.45, -100000, 1.85, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(0.663, -100000, 0.551, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)


plot(seq1, sapply(tmp, function(x) x$lnLik))
save.image("cyanoarchea.Rda")

seq2 <- seq(3.7, 0.1, -0.1)
profresFull <- foreach(i=1:22) %dopar% {
  tmp <- list()
  for(j in 1:length(seq2)) {
    fn <- fns[[i]]
    ft <<- seq2[j]
    fn <- constrain2(fn, q10.tc~ft, extra=c("r1", "r2"))
    if(j==1){
      startx <- runif(4, 0.5, 2)*c(1, 1, cyFit.best.notime[i, 1], cyFit.best.notime[i, 2])
    } else {
      startx <- runif(4, 0.5, 2)*c(1, 1, cyFit.best.notime[i, 1], cyFit.best.notime[i, 2]) #c(tmp[[j-1]]$par)
    }
    tmp[[j]] <- find.mle(fn, x.init=startx, method="optim")
  }
  tmp
}
#tmpfn <- function(i, j, ...) {
  tmp <- list()
  fn <- fns[[i]]
  ft <<- seq2[j]
  fn <- constrain2(fn, q10.tc~ft, extra=c("r1", "r2"))
  if(j==1){
    startx <- runif(4, 0.5, 2)*c(1, 1, cyFit.best.notime[i, 1], cyFit.best.notime[i, 2])
  } else {
    startx <- runif(4, 0.5, 2)*c(1, 1, cyFit.best.notime[i, 1], cyFit.best.notime[i, 2]) #c(tmp[[j-1]]$par)
  }
  tmp[[j]] <- find.mle(fn, x.init=startx, ...)
  tmp
}

profliksFull <- sapply(1:22, function(x) sapply(profresFull[[x]], function(y) y$lnLik))
plot(seq1, apply(profliks, 1, sum), ylim=c(-2960, -2910))
lines(seq1,  apply(profliks, 1, sum))
#plot(seq2, apply(profliksFull, 1, sum))
lines(seq2,  apply(profliksFull, 1, sum))
rect(2.45, -100000, 1.85, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(0.663, -100000, 0.551, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)

## Research and replace better log-likelihood fits using the subplex routine.
seq2 <- seq(3.7, 0.1, -0.1)
researchFull <- foreach(i=1:22) %dopar% {
  tmp <- list()
  for(j in 1:length(seq2)) {
    fn <- fns[[i]]
    ft <<- seq2[j]
    fn <- constrain2(fn, q10.tc~ft, extra=c("r1", "r2"))
    if(j==1){
      startx <- runif(4, 0.5, 2)*c(1, 1, cyFit.best.notime[i, 1], cyFit.best.notime[i, 2])
    } else {
      startx <- runif(4, 0.5, 2)*c(1, 1, cyFit.best.notime[i, 1], cyFit.best.notime[i, 2]) #c(tmp[[j-1]]$par)
    }
    ftmp <- find.mle(fn, x.init=startx, method="subplex")
    if(ftmp$lnLik > profresFull[[i]][[j]]$lnLik){
      tmp[[j]] <- ftmp
    } else {
      tmp[[j]] <- profresFull[[i]][[j]]
    }
  }
  tmp
}
## And again:
fns <- ARD.R2.fns
registerDoMC(cores=8)
researchFull <- foreach(i=1:22) %dopar% {
  tmp <- list()
  for(j in 1:length(seq2)) {
    fn <- fns[[i]]
    ft <<- seq2[j]
    fn <- constrain2(fn, q10.tc~ft, extra=c("r1", "r2"))
    if(j==1){
      startx <- unlist(cyFit.best[i,1:4])
      } else {
      startx <- unlist(cyFit.best[i,1:4])
    }
    ftmp <- find.mle(fn, x.init=startx, method="subplex")
    if(ftmp$lnLik > researchFull[[i]][[j]]$lnLik){
      tmp[[j]] <- ftmp
    } else {
      tmp[[j]] <- researchFull[[i]][[j]]
    }
  }
  tmp
}

## And again:
fns <- ARD.R2.fns
registerDoMC(cores=8)
researchFull <- foreach(i=1:22) %dopar% {
  tmp <- list()
  for(j in 1:length(seq2)) {
    fn <- fns[[i]]
    ft <<- seq2[j]
    fn <- constrain2(fn, q10.tc~ft, extra=c("r1", "r2"))
    if(j==1){
      startx <- unlist(cyFit.best.GO1[i,1:4])
    } else {
      startx <- unlist(cyFit.best.GO1[i,1:4])
    }
    ftmp <- find.mle(fn, x.init=startx, method="subplex")
    if(ftmp$lnLik > researchFull[[i]][[j]]$lnLik){
      tmp[[j]] <- ftmp
    } else {
      tmp[[j]] <- researchFull[[i]][[j]]
    }
  }
  tmp
}

## And again:
fns <- ARD.R2.fns
registerDoMC(cores=8)
researchFull <- foreach(i=1:22) %dopar% {
  tmp <- list()
  for(j in 1:length(seq2)) {
    fn <- fns[[i]]
    ft <<- seq2[j]
    fn <- constrain2(fn, q10.tc~ft, extra=c("r1", "r2"))
    if(j==1){
      startx <- unlist(cyFit.best.GO2[i,1:4])
    } else {
      startx <- unlist(cyFit.best.GO2[i,1:4])
    }
    ftmp <- find.mle(fn, x.init=startx, method="subplex")
    if(ftmp$lnLik > researchFull[[i]][[j]]$lnLik){
      tmp[[j]] <- ftmp
    } else {
      tmp[[j]] <- researchFull[[i]][[j]]
    }
  }
  tmp
}


resliksFull <- sapply(1:22, function(x) sapply(researchFull[[x]], function(y) y$lnLik))


plot(seq2, apply(resliksFull, 1, sum), xlim=c(3.8, 0), ylim=c(-2920, -2900))
#lines(seq1,  apply(profliks, 1, sum))
##plot(seq2, apply(profliksFull, 1, sum))
#lines(seq2,  apply(profliksFull, 1, sum))
lines(seq2, apply(resliksFull, 1, sum), col="red", lwd=2)
rect(2.45, -100000, 1.85, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(0.663, -100000, 0.551, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)

png("cumlnLcyano4.png", width=1000)
## Plot for all traits
plot(seq2, apply(resliksFull, 1, sum),type="n", xlim=c(3.8, -.8), ylim=c(0,35), bty="l", xlab="Billion years before present", ylab="Support")
#lines(seq1,  apply(profliks, 1, sum))
##plot(seq2, apply(profliksFull, 1, sum))
#lines(seq2,  apply(profliksFull, 1, sum))
resfullN <- apply(resliksFull, 2, function(x) x-min(x))
colnames(resfullN) <- colnames(tdcy$dat)[1:22]
resfullN <- resfullN[,order(apply(resfullN, 2, max))]
cumres <- sapply(1:nrow(resfullN), function(x) cumsum(resfullN[x, ]))
lapply(1:22, function(y) lines(seq2,cumres[y,], col=y, lwd=2))
rect(2.45, -100000, 1.85, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(0.663, -100000, 0.551, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
text(rep(0.1, 22), cumres[,ncol(cumres)], labels=colnames(resfullN), cex=0.5, pos=4)
max.lik <- seq2[which(cumres[22,]==max(cumres[22,]))]
abline(v=max.lik)
abline(v=0.6)
dev.off()

## Try ASR's with fitted ARD.R2 model
require(plotrix)
require(grid)
require(bayou)
pdf("asrRateShifts2.pdf", width=12, height=8)
asrFigure <- function(researchFull) {## Plot for all traits

  resliksFull <- sapply(1:nc, function(x) sapply(researchFull[[x]], function(y) y$lnLik))
  plot(seq2, apply(resliksFull, 1, sum),type="n", xlim=c(3.8, -.8), ylim=c(0,40), bty="l", xlab="Billion years before present", ylab="Support")
  #lines(seq1,  apply(profliks, 1, sum))
  ##plot(seq2, apply(profliksFull, 1, sum))
  #lines(seq2,  apply(profliksFull, 1, sum))
  resfullN <- apply(resliksFull, 2, function(x) x-min(x))
  colnames(resfullN) <- colnames(tdcy$dat)
  o <- order(apply(resfullN, 2, max))
  resfullN <- resfullN[,o]
  cumres <- sapply(1:nrow(resfullN), function(x) cumsum(resfullN[x, ]))
  lapply(1:nc, function(y) lines(seq2,cumres[y,], col=y, lwd=2))
  rect(2.45, -100000, 1.85, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
  rect(0.663, -100000, 0.551, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
  text(rep(0.1, 22), cumres[,ncol(cumres)], labels=colnames(resfullN), cex=0.5, pos=4)
  max.lik <- seq2[which(cumres[22,]==max(cumres[22,]))]
  abline(v=max.lik)
  abline(v=0.6)

  par(mfrow=c(1,2))
  for(i in rev(o)){
    basicARDfn <- make.bisse.t(tdcyList[[i]]$phy, setNames(tdcyList[[i]]$dat[,1], tdcyList[[i]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
    pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
    asr <- list(t(asr.marginal(basicARDfn, pars)))
    class(asr) <- c("asrArbor", "list")
    attributes(asr)$charType <- "discrete"
    attributes(asr)$td <- tdcyList[[i]]
    attributes(asr)$na.drop <- NULL
    attributes(asr)$charStates <- list(c("0","1"))
    plot(asr, type="fan", cex=0.5, label.offset=0.5, pal=colorRampPalette(c("green", "blue")))
    draw.circle(0, 0, TH-2.15, lwd=20, border=makeTransparent('gray80', alpha=100))
    draw.circle(0, 0, TH-0.6, lwd=4, border=makeTransparent('gray80', alpha=100))
    draw.circle(0, 0, TH-seq2[which(resliksFull[,i]==max(resliksFull[,i]))], lwd=3, lty=2, border='red')
    plot(seq2, resliksFull[,i],type="n", xlim=c(3.8, -.8), bty="l", xlab="Billion years before present", ylab="Support")
    lines(seq2, resliksFull[,i], col="red", lwd=2)
    rect(2.45, -1000, 1.85, 1000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
    rect(0.663, -1000, 0.551, 1000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
    text(3, max(resliksFull[,i]), label=paste("max(dlnL) =", round(max(resliksFull[,i]-min(resliksFull[,i])),2)), cex=1.25)
    researchFull[[i]]
  }
}
asrFigure(researchFull)
dev.off()

par(mfrow=c(2,2))
plot(asrDiscrete, show.tip.label=FALSE)
rect(3.8-2.45, -100000, 3.8-1.85, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(3.8-0.663, -100000,3.8- 0.551, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)

Rs <- lapply(researchFull, function(x) t(sapply(x, function(y) c(y$par.full[c(7, 5:6, 8:9)], round(y$par.full[5]/y$par.full[6],2), round(y$par.full[8]/y$par.full[9],2), y$lnL))))

best.Rs <- t(sapply(Rs, function(x) x[which(x[,8]==max(x[,8])),]))
tmp <- cbind(best.Rs[,1], log(best.Rs[,6]+0.000001), log(best.Rs[,7]+0.0000001))
rownames(tmp) <- colnames(tdcy$dat)[1:26]
subset(tmp, tmp[,1] > 1.5)
subset(tmp, tmp[,1] < 0.8)


## Now add in the 3 extra traits
names(researchFull) <- c('Thermophilic', 'Nonfreshwater_habitat', 'Akinetes', 'Heterocysts', 'Nitrogen_fixation', 'Morphology',
                         'Habit', 'Freeliving', 'Mats', 'Epi/Endolithic',  'Epiphytic', 'Periphytic', 'Motility', 'Hormogonia', 
                         'Gas_vesicles', 'False_Branching', 'True_Branching', 'Fission_in_multiple_planes', 'Seriate_trichomes',
                         'Baeocytes', 'Extracellular_sheath', 'Mucilage')
new.researchFull <- lapply(1:25, function(x) lapply(1:37, function(y) list(lnL=NA)))
names(new.researchFull) <- colnames(tdcy$dat)

new.researchFull[names(researchFull)] <- researchFull
starts <- lapply(1:25, function(x) lapply(1:37, function(y) c(exp(runif(2, -3, 3)), exp(runif(2, -3, 3)))))

newtraits <- which(!(names(new.researchFull) %in% names(researchFull)))

registerDoMC(cores=3)
fns <- ARD.R2.fns
sequence <- seq2
trait.seq <- sample(1:25, 3, replace=FALSE)
previous <- new.researchFull
previousliks <- sapply(1:nc, function(x) sapply(previous[[x]], function(y) y$lnL))
extra <- c("r1", "r2")
starts <- lapply(1:nc, function(x) lapply(1:37, function(y) previous[[x]][[1]]$par))

na.comp <- function(a, b){
  if(is.na(b)){
    TRUE
  } else {
    a > b
  }
}
res <- foreach(i=trait.seq) %dopar% {
  tmp <- list()
  for(j in 1:length(sequence)) {
    fn <- fns[[i]]
    ft <<- sequence[j]
    fn <- constrain2(fn, q10.tc~ft, extra=extra)
    startx <- starts[[i]][[j]]
    if(is.null(previous)){
      tmp[[j]] <- find.mle(fn, x.init=startx, method="subplex")
    } else{
      ftmp <- find.mle(fn, x.init=startx, method="subplex")
      if(na.comp(ftmp$lnL, previous[[i]][[j]]$lnL)){
        tmp[[j]] <- ftmp
        print(paste("Found a better lnL for trait", i, "at", sequence[j], "bybp", sep=" "))
      } else {
        tmp[[j]] <- previous[[i]][[j]]
      }
    }
  }
  tmp
}

new.researchFull[trait.seq] <- res
pdf("ASRFigure25traits.pdf", width=16, height=8)
asrFigure(new.researchFull)
dev.off()

sapply(tdcyList, function(x) colnames(x$dat)[1])

save(new.researchFull, file="new.researchFull.rds")


stackedLiks <- function(researchFull, pal=rainbow, ...){
  resliksFull <- sapply(1:nc, function(x) sapply(researchFull[[x]], function(y) y$lnLik))
  plot(seq2, apply(resliksFull, 1, sum),type="n", xlim=c(3.8, 0), ylim=c(0,40), xlab="Billion years before present", ylab=expression(paste(Delta, " Ln likelihood", sep="")),...)
  rect(2.45, 0, 1.85, 41, col="gray40", border=NA)
  rect(0.663, 0, 0.551, 41, col="gray40", border=NA)
  #lines(seq1,  apply(profliks, 1, sum))
##plot(seq2, apply(profliksFull, 1, sum))
#lines(seq2,  apply(profliksFull, 1, sum))
  resfullN <- apply(resliksFull, 2, function(x) x-min(x))
  colnames(resfullN) <- colnames(tdcy$dat)
  o <- order(apply(resfullN, 2, max))
  resfullN <- resfullN[,o]
  cumres <- sapply(1:nrow(resfullN), function(x) cumsum(resfullN[x, ]))
  lapply(1:nc, function(y) lines(seq2,cumres[y,], col=pal(nc)[y], lwd=2))
  rect(2.45, -100000, 1.85, 10000, col=rgb(255,255,255, maxColorValue=255,alpha=75), border=NA)
  rect(0.663, -100000, 0.551, 10000, col=rgb(255,255,255, maxColorValue=255,alpha=75), border=NA)
  #text(rep(0.1, 25), cumres[,ncol(cumres)]+seq(0, 1, length.out=25), labels=colnames(resfullN), cex=0.5, pos=4, col=pal(nc))
  max.lik <- seq2[which(cumres[25,]==max(cumres[25,]))]
  abline(v=1.6, col="orange", lwd=2, lty=2)
  abline(v=0.6, col="orange", lwd=2,lty=2)
}

pdf("../docs/images/Stackedliks.pdf", width=10, height=8)
par(bg="black", mar=c(5,5, 3,3))
stackedLiks(new.researchFull, colorRampPalette(c("green", "lightgreen", "cyan", "violet")), cex.axis=1.5, cex.lab=2, col.axis="lightgreen", col.lab="lightgreen", col.sub="lightgreen", fg="lightgreen")
dev.off()

findpeaks <- function(researchFull, window){
  resliksFull <- sapply(1:nc, function(x) sapply(researchFull[[x]], function(y) y$lnLik))
  resfullN <- apply(resliksFull, 2, function(x) x-min(x))
  colnames(resfullN) <- colnames(tdcy$dat)
  maxes <- apply(resfullN, 2, max)
  whichmax <- apply(resfullN, 2, function(x) which(x==max(x)))
  maxwindow <-sapply(whichmax, function(x) as.numeric(seq2[x] >= window[1] & seq2[x] <= window[2]))
  return(maxwindow)
}

peak0.6 <- findpeaks(new.researchFull, c(0.4, 0.8))
peak0.0 <- findpeaks(new.researchFull, c(0, 0.3))
peak1.8 <- findpeaks(new.researchFull, c(1.5, 2))
peak1.0 <- findpeaks(new.researchFull, c(1.1, 1.4))
bimodal <- peak0.6
bimodal[1:length(bimodal)] <- 0
bimodal[c('Motility', 'celldiam_mean', 'Gas_vesicles', 'Nitrogen_fixation', 'Nonfreshwater_habitat', 'True_Branching', 'False_Branching')] <- 1

resliksFull <- sapply(1:nc, function(x) sapply(researchFull[[x]], function(y) y$lnLik))
resfullN <- apply(resliksFull, 2, function(x) x-min(x))
colnames(resfullN) <- colnames(tdcy$dat)

require(ggplot2)
require(reshape2)
peaks <- bimodal
peakplot<- function(peaks, pal=colorRampPalette(c("lightgreen", "green", "cyan", "blue","yellow", "orange", "pink","darkred"))){
  resliksN <- resfullN
  colnames(resliksN) <- gsub("_", " ", colnames(tdcy$dat))
  o <- rev(order(apply(resliksN, 2, max)[as.logical(peaks)]))
  o.all <- order(apply(resliksN, 2, max))
  resliksN <- resliksN[,o.all]
  peaks <- peaks[o.all]
  df <- melt(resliksN[,as.logical(peaks)])
  df[,1] <- seq2
  cols <- pal(nc)
  df$grid <- factor(df$Var2, levels=rev(factor((colnames(resliksN)[as.logical(peaks)]))))
  p <- ggplot(df, aes(Var1, value, colour=factor(Var2, levels=(colnames(resliksN)[as.logical(peaks)]))))
  p <- p + scale_colour_manual(values=cols[as.logical(peaks)])
  p <- p + geom_vline(xintercept=0.6, linetype=2, colour="lightgreen", size=1)
  p <- p + geom_vline(xintercept=1.8, linetype=2, colour="lightgreen", size=1)
  p <- p + geom_line(size=1.5)
  p <- p + xlim(c(3.8,0))
  p <- p + theme_bw()
  p <- p + theme(strip.background=element_rect(fill="black", colour="lightgreen", size=1),
               line= element_line(colour="lightgreen", size=1),
               legend.position = "none",
               plot.background=element_rect("black"),
               panel.background=element_rect("black"),
               panel.grid.major=element_blank(),
               panel.border=element_rect(colour="lightgreen", size=1),
               axis.line=element_line(colour="lightgreen", size=1),
               text=element_text(colour="lightgreen", size=14, face="bold"),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank()
              )
  p <- p + xlab("Gyr")
  p <- p + ylab("Support")
  p <- p + facet_grid(grid~.,  scales="free_y")
  p
}
pal=colorRampPalette(c("green", "lightgreen", "cyan", "violet"))
pdf("../docs/images/GO1shifts.pdf", height=sum(peak1.8)*2, width=8)
peaks <- peak1.8
{resliksN <- resfullN
 colnames(resliksN) <- gsub("_", " ", colnames(tdcy$dat))
 o <- rev(order(apply(resliksN, 2, max)[as.logical(peaks)]))
 o.all <- order(apply(resliksN, 2, max))
 resliksN <- resliksN[,o.all]
 peaks <- peaks[o.all]
 df <- melt(resliksN[,as.logical(peaks)])
 df[,1] <- seq2
 cols <- pal(nc)
 df$grid <- factor(df$Var2, levels=rev(factor((colnames(resliksN)[as.logical(peaks)]))))
 p <- ggplot(df, aes(Var1, value, colour=factor(Var2, levels=(colnames(resliksN)[as.logical(peaks)]))))
 p <- p + scale_colour_manual(values=cols[as.logical(peaks)])
 p <- p + geom_vline(xintercept=0.6, linetype=2, colour="lightgreen", size=1)
 p <- p + geom_vline(xintercept=1.8, linetype=2, colour="lightgreen", size=1)
 p <- p + geom_line(size=1.5)
 p <- p + xlim(c(3.8,0))
 p <- p + theme_bw()
 p <- p + theme(strip.background=element_rect(fill="black", colour="lightgreen", size=1),
                text=element_text(colour="lightgreen", size=14, face="bold"),
                line= element_line(colour="lightgreen", size=1),
                legend.position = "none",
                plot.background=element_rect(fill="black", colour="black"),
                panel.background=element_rect("black"),
                panel.grid.major=element_blank(),
                axis.title=element_text(size=24),
                panel.border=element_rect(colour="lightgreen", size=1),
                axis.line=element_line(colour="lightgreen", size=1),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()
 )
 p <- p + xlab("Gyr")
 p <- p + ylab("Support")
 p <- p + facet_grid(grid~.,  scales="free_y")
 p
}
dev.off()

pdf("../docs/images/GO2shifts.pdf", height=sum(peak0.6)*2, width=8)
peaks <- peak0.6
{resliksN <- resfullN
 colnames(resliksN) <- gsub("_", " ", colnames(tdcy$dat))
 o <- rev(order(apply(resliksN, 2, max)[as.logical(peaks)]))
 o.all <- order(apply(resliksN, 2, max))
 resliksN <- resliksN[,o.all]
 peaks <- peaks[o.all]
 df <- melt(resliksN[,as.logical(peaks)])
 df[,1] <- seq2
 cols <- pal(nc)
 df$grid <- factor(df$Var2, levels=rev(factor((colnames(resliksN)[as.logical(peaks)]))))
 p <- ggplot(df, aes(Var1, value, colour=factor(Var2, levels=(colnames(resliksN)[as.logical(peaks)]))))
 p <- p + scale_colour_manual(values=cols[as.logical(peaks)])
 p <- p + geom_vline(xintercept=0.6, linetype=2, colour="lightgreen", size=1)
 p <- p + geom_vline(xintercept=1.8, linetype=2, colour="lightgreen", size=1)
 p <- p + geom_line(size=1.5)
 p <- p + xlim(c(3.8,0))
 p <- p + theme_bw()
 p <- p + theme(strip.background=element_rect(fill="black", colour="lightgreen", size=1),
                text=element_text(colour="lightgreen", size=14, face="bold"),
                line= element_line(colour="lightgreen", size=1),
                legend.position = "none",
                plot.background=element_rect(fill="black", colour="black"),
                panel.background=element_rect("black"),
                panel.grid.major=element_blank(),
                axis.title=element_text(size=24),
                panel.border=element_rect(colour="lightgreen", size=1),
                axis.line=element_line(colour="lightgreen", size=1),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()
 )
 p <- p + xlab("Gyr")
 p <- p + ylab("Support")
 p <- p + facet_grid(grid~.,  scales="free_y")
 p
}
dev.off()

pdf("../docs/images/bimodalshifts.pdf", height=sum(bimodal)*2, width=8)
peaks <- bimodal
{resliksN <- resfullN
 colnames(resliksN) <- gsub("_", " ", colnames(tdcy$dat))
 o <- rev(order(apply(resliksN, 2, max)[as.logical(peaks)]))
 o.all <- order(apply(resliksN, 2, max))
 resliksN <- resliksN[,o.all]
 peaks <- peaks[o.all]
 df <- melt(resliksN[,as.logical(peaks)])
 df[,1] <- seq2
 cols <- pal(nc)
 df$grid <- factor(df$Var2, levels=rev(factor((colnames(resliksN)[as.logical(peaks)]))))
 p <- ggplot(df, aes(Var1, value, colour=factor(Var2, levels=(colnames(resliksN)[as.logical(peaks)]))))
 p <- p + scale_colour_manual(values=cols[as.logical(peaks)])
 p <- p + geom_vline(xintercept=0.6, linetype=2, colour="lightgreen", size=1)
 p <- p + geom_vline(xintercept=1.8, linetype=2, colour="lightgreen", size=1)
 p <- p + geom_line(size=1.5)
 p <- p + xlim(c(3.8,0))
 p <- p + theme_bw()
 p <- p + theme(strip.background=element_rect(fill="black", colour="lightgreen", size=1),
                text=element_text(colour="lightgreen", size=14, face="bold"),
                line= element_line(colour="lightgreen", size=1),
                legend.position = "none",
                plot.background=element_rect(fill="black", colour="black"),
                panel.background=element_rect("black"),
                panel.grid.major=element_blank(),
                axis.title=element_text(size=24),
                panel.border=element_rect(colour="lightgreen", size=1),
                axis.line=element_line(colour="lightgreen", size=1),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()
 )
 p <- p + xlab("Gyr")
 p <- p + ylab("Support")
 p <- p + facet_grid(grid~.,  scales="free_y")
 p
}
dev.off()


ratebox <- function(rates, center, halfwidth, pal, n=100, log=TRUE, range=NULL){
  cols <- pal(n+2)
  rr <- log(rates)
  if(is.null(range)){
    range <- c(min(rr), max(rr))
  }
  rr <- rr-range[1]+1
  rr <- round(rr/(range[2]-range[1]+1)*100,0)
  rect(center[1], center[2], center[1]+halfwidth, center[2]+halfwidth, col=cols[rr[1]], border=cols[rr[1]])
  rect(center[1]-halfwidth, center[2], center[1], center[2]+halfwidth, col=cols[rr[2]], border=cols[rr[2]])
  rect(center[1], center[2]-halfwidth, center[1]+halfwidth, center[2], col=cols[rr[3]], border=cols[rr[3]])
  rect(center[1]-halfwidth, center[2]-halfwidth, center[1], center[2], col=cols[rr[4]], border=cols[rr[4]])
  rect(center[1]-halfwidth, center[2]-halfwidth, center[1]+halfwidth, center[2]+halfwidth, border="lightgreen", lwd=2)
}
circlerates <- function(rates, cex, center, dists, pal,adj,yadj, n=100,widthscale=10, log=TRUE, range=NULL,...){
  cols <- pal(n+2)
  rr <- log(rates)
  if(is.null(range)){range <- c(min(rr), max(rr))}
  rr <- rr-range[1]+1
  rr <- ifelse(rr < 1, 1, rr)
  rr <- ifelse(rr > (range[2]-range[1]), range[2], rr)
  rr <- round(rr/(range[2]-range[1]+1)*100,0)
  arrows(center[1]-1*dists[1] + adj,center[2]+yadj,  center[1]-1*dists[1]+dists[2]-adj, center[2]+yadj, col=cols[rr[2]], lwd=rr[2]/(n/widthscale+0.1), ...)
  arrows(center[1]-1*dists[1]+dists[2]-adj,center[2]-yadj, center[1]-1*dists[1] + adj, center[2]-yadj, col=cols[rr[4]], lwd=rr[4]/(n/widthscale+0.1), ...)
  arrows(center[1]+1*dists[1]-adj,center[2]-yadj,  center[1]+1*dists[1]-dists[2]+adj, center[2]-yadj, col=cols[rr[3]], lwd=rr[3]/(n/widthscale+0.1), ...)
  arrows(center[1]+1*dists[1]-dists[2]+adj,center[2]+yadj, center[1]+1*dists[1]-adj, center[2]+yadj, col=cols[rr[1]], lwd=rr[1]/(n/widthscale+0.1), ...)
  #points(c(center[1]-1*dists[1], center[1]-1*dists[1]+dists[2], center[1]+dists[1], center[1]+dists[1]-dists[2]), rep(center[2],4), pch=21, cex=cex, bg=c("lightgreen", "red", "red", "lightgreen"), col=c("lightgreen", "red", "red", "lightgreen"))
  text(c(center[1]-1*dists[1], center[1]-1*dists[1]+dists[2], center[1]+dists[1], center[1]+dists[1]-dists[2]), rep(center[2],4), labels=c("0", "1", "1", "0"), col="lightgreen", cex=cex)
}

rates <- sapply(new.researchFull[as.logical(peak1.8)], function(x) x[[20]]$par.full)[c(5,6,8,9),]
nn=ncol(rates)
      pdf("../docs/images/rateboxesGO1.pdf", height=sum(peak1.8)*2, width=8)
      par(bg="black")
      plot(c(-0.8,0.8), c(0,nn), type="n",xlab="",ylab="", xaxt="n", yaxt="n", bty="n")
      #lapply(1:nn, function(x) ratebox(GO1rates[,x], c(0, x), 0.1, colorRampPalette(c("blue","blue", "white", "violet")), n=100, log=TRUE, range=c(min(log(GO1rates)), max(log(GO1rates)))))
      ymarks <- 1:nn
      lapply(1:nn, function(x) circlerates(rates[,x],cex=4, c(0, ymarks[x]), dists=c(0.4, 0.3), colorRampPalette(c("blue","blue", "violet", "hotpink")), adj=0.08, yadj=0.1, widthscale=6, range=c(min(log(rates)), max(log(rates)))))
      
dev.off()
  
resliksFull <- sapply(1:25, function(x) sapply(new.researchFull[[x]], function(y) y$lnLik))
resfullN <- apply(resliksFull, 2, function(x) x-min(x))
colnames(resfullN) <- colnames(tdcy$dat)
maxes <- apply(resfullN, 2, max)
o <- order(maxes, decreasing=FALSE)
rates <- sapply(new.researchFull[o][as.logical(peak0.6[o])], function(x) x[[32]]$par.full)[c(5,6,8,9),]
nn=ncol(rates)
pdf("../docs/images/rateboxesGO2.pdf", height=sum(peak0.6)*2, width=8)
par(bg="black")
plot(c(-0.8,0.8), c(0,nn), type="n",xlab="",ylab="", xaxt="n", yaxt="n", bty="n")
#lapply(1:nn, function(x) ratebox(GO1rates[,x], c(0, x), 0.1, colorRampPalette(c("blue","blue", "white", "violet")), n=100, log=TRUE, range=c(min(log(GO1rates)), max(log(GO1rates)))))
ymarks <- 1:nn
lapply(1:nn, function(x) circlerates(rates[,x],cex=4, c(0, ymarks[x]), dists=c(0.4, 0.3), colorRampPalette(c("blue", "violet", "hotpink", "white")), adj=0.08, yadj=0.1, widthscale=10, range=c(-3, 2)))
dev.off()


## ASR's for extracellular sheath and Planktonic vs. Benthic/Sessile
pdf("../docs/images/ecsheath.pdf", height=10, width=8)
par(bg="black")
trait <- "Extracellular_sheath"
{
i = which(names(tdcy$dat)==trait)
resliksFull <- sapply(1:25, function(x) sapply(new.researchFull[[x]], function(y) y$lnLik))
colnames(resliksFull) <- colnames(tdcy$dat)
basicARDfn <- make.bisse.t(tdcyList[[i]]$phy, setNames(tdcyList[[i]]$dat[,1], tdcyList[[i]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
asr <- list(t(asr.marginal(basicARDfn, pars)))
tmptd <- tdcyList[[i]]
colnames(tmptd$dat)[1] <- ""
charType <- "discrete"
na.drop <- NULL
charStates <- (c("0","1"))
plot(tmptd$phy, cex=0.5, label.offset=0.1, edge.color="lightgreen", edge.width=3, tip.color="lightgreen", x.lim=c(0, 6))
pal=colorRampPalette(c("green", "red"))
nodelabels(pie=asr[[1]], piecol=pal(length(charStates)), cex=0.6, frame="circle")
tiplabels(pch=21, bg=pal(length(charStates))[as.numeric(factor(tmptd$dat[,1], levels=charStates))]) 
#legend("bottomleft", fill=pal(length(charStates)), legend=charStates)
abline(v=TH-0.5, col="orange", lwd=3, lty=2)
}
dev.off()
pdf("../docs/images/habit.pdf", height=10, width=8)
par(bg="black")
trait <- "Habit"
{
  i = which(names(tdcy$dat)==trait)
  resliksFull <- sapply(1:25, function(x) sapply(new.researchFull[[x]], function(y) y$lnLik))
  colnames(resliksFull) <- colnames(tdcy$dat)
  basicARDfn <- make.bisse.t(tdcyList[[i]]$phy, setNames(tdcyList[[i]]$dat[,1], tdcyList[[i]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
  pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
  asr <- list(t(asr.marginal(basicARDfn, pars)))
  tmptd <- tdcyList[[i]]
  colnames(tmptd$dat)[1] <- ""
  charType <- "discrete"
  na.drop <- NULL
  charStates <- (c("0","1"))
  plot(tmptd$phy, cex=0.5, label.offset=0.1, edge.color="lightgreen", edge.width=3, tip.color="lightgreen", x.lim=c(0, 6))
  pal=colorRampPalette(c("green", "red"))
  nodelabels(pie=asr[[1]], piecol=pal(length(charStates)), cex=0.6, frame="circle")
  tiplabels(pch=21, bg=pal(length(charStates))[as.numeric(factor(tmptd$dat[,1], levels=charStates))]) 
  #legend("bottomleft", fill=pal(length(charStates)), legend=charStates)
  abline(v=TH-0.5, col="orange", lwd=3, lty=2)
}
dev.off()
pdf("../docs/images/pelagic.pdf", height=10, width=8)
par(bg="black")
trait <- "Pelagic"
{
  i = which(names(tdcy$dat)==trait)
  resliksFull <- sapply(1:25, function(x) sapply(new.researchFull[[x]], function(y) y$lnLik))
  colnames(resliksFull) <- colnames(tdcy$dat)
  basicARDfn <- make.bisse.t(tdcyList[[i]]$phy, setNames(tdcyList[[i]]$dat[,1], tdcyList[[i]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
  pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
  asr <- list(t(asr.marginal(basicARDfn, pars)))
  tmptd <- tdcyList[[i]]
  colnames(tmptd$dat)[1] <- ""
  charType <- "discrete"
  na.drop <- NULL
  charStates <- (c("0","1"))
  plot(tmptd$phy, cex=0.5, label.offset=0.1, edge.color="lightgreen", edge.width=3, tip.color="lightgreen", x.lim=c(0, 6))
  pal=colorRampPalette(c("green", "red"))
  nodelabels(pie=asr[[1]], piecol=pal(length(charStates)), cex=0.6, frame="circle")
  tiplabels(pch=21, bg=pal(length(charStates))[as.numeric(factor(tmptd$dat[,1], levels=charStates))]) 
  #legend("bottomleft", fill=pal(length(charStates)), legend=charStates)
  abline(v=TH-0.5, col="orange", lwd=3, lty=2)
}
dev.off()

pdf("../docs/images/fwhabitat.pdf", height=10, width=8)
par(bg="black")
trait <- "Nonfreshwater_habitat"
{
  i = which(names(tdcy$dat)==trait)
  resliksFull <- sapply(1:25, function(x) sapply(new.researchFull[[x]], function(y) y$lnLik))
  colnames(resliksFull) <- colnames(tdcy$dat)
  basicARDfn <- make.bisse.t(tdcyList[[i]]$phy, setNames(tdcyList[[i]]$dat[,1], tdcyList[[i]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
  pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
  asr <- list(t(asr.marginal(basicARDfn, pars)))
  tmptd <- tdcyList[[i]]
  colnames(tmptd$dat)[1] <- ""
  charType <- "discrete"
  na.drop <- NULL
  charStates <- (c("0","1"))
  plot(tmptd$phy, cex=0.5, label.offset=0.1, edge.color="lightgreen", edge.width=3, tip.color="lightgreen", x.lim=c(0, 6))
  pal=colorRampPalette(c("green", "red"))
  nodelabels(pie=asr[[1]], piecol=pal(length(charStates)), cex=0.6, frame="circle")
  tiplabels(pch=21, bg=pal(length(charStates))[as.numeric(factor(tmptd$dat[,1], levels=charStates))]) 
  #legend("bottomleft", fill=pal(length(charStates)), legend=charStates)
  abline(v=TH-0.5, col="orange", lwd=3, lty=2)
}
dev.off()

  pdf("../docs/images/heterocysts.pdf", height=10, width=8)
  par(bg="black")
  trait <- "Heterocysts"
  {
    i = which(names(tdcy$dat)==trait)
    resliksFull <- sapply(1:25, function(x) sapply(new.researchFull[[x]], function(y) y$lnLik))
    colnames(resliksFull) <- colnames(tdcy$dat)
    basicARDfn <- make.bisse.t(tdcyList[[i]]$phy, setNames(tdcyList[[i]]$dat[,1], tdcyList[[i]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
    pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
    asr <- list(t(asr.marginal(basicARDfn, pars, root=ROOT.FLAT)))
    tmptd <- tdcyList[[i]]
    colnames(tmptd$dat)[1] <- ""
    charType <- "discrete"
    na.drop <- NULL
    charStates <- (c("0","1"))
    plot(tmptd$phy, cex=0.5, label.offset=0.1, edge.color="lightgreen", edge.width=3, tip.color="lightgreen", x.lim=c(0, 6))
    pal=colorRampPalette(c("green", "red"))
    nodelabels(pie=asr[[1]], piecol=pal(length(charStates)), cex=0.6, frame="circle")
    tiplabels(pch=21, bg=pal(length(charStates))[as.numeric(factor(tmptd$dat[,1], levels=charStates))]) 
    #legend("bottomleft", fill=pal(length(charStates)), legend=charStates)
    abline(v=TH-1.8, col="orange", lwd=3, lty=2)
  }
  dev.off()

pdf("../docs/images/Morphology.pdf", height=10, width=8)
par(bg="black")
trait <- "Morphology"
{
  i = which(names(tdcy$dat)==trait)
  resliksFull <- sapply(1:25, function(x) sapply(new.researchFull[[x]], function(y) y$lnLik))
  colnames(resliksFull) <- colnames(tdcy$dat)
  basicARDfn <- make.bisse.t(tdcyList[[i]]$phy, setNames(tdcyList[[i]]$dat[,1], tdcyList[[i]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
  pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
  asr <- list(t(asr.marginal(basicARDfn, pars, root=ROOT.FLAT)))
  tmptd <- tdcyList[[i]]
  colnames(tmptd$dat)[1] <- ""
  charType <- "discrete"
  na.drop <- NULL
  charStates <- (c("0","1"))
  plot(tmptd$phy, cex=0.5, label.offset=0.1, edge.color="lightgreen", edge.width=3, tip.color="lightgreen", x.lim=c(0, 6))
  pal=colorRampPalette(c("green", "red"))
  nodelabels(pie=asr[[1]], piecol=pal(length(charStates)), cex=0.6, frame="circle")
  tiplabels(pch=21, bg=pal(length(charStates))[as.numeric(factor(tmptd$dat[,1], levels=charStates))]) 
  #legend("bottomleft", fill=pal(length(charStates)), legend=charStates)
  abline(v=TH-1.8, col="orange", lwd=3, lty=2)
}
dev.off()

pdf("../docs/images/Motility.pdf", height=10, width=8)
par(bg="black")
trait <- "Motility"
{
  i = which(names(tdcy$dat)==trait)
  resliksFull <- sapply(1:25, function(x) sapply(new.researchFull[[x]], function(y) y$lnLik))
  colnames(resliksFull) <- colnames(tdcy$dat)
  basicARDfn <- make.bisse.t(tdcyList[[i]]$phy, setNames(tdcyList[[i]]$dat[,1], tdcyList[[i]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
  pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
  asr <- list(t(asr.marginal(basicARDfn, pars, root=ROOT.FLAT)))
  tmptd <- tdcyList[[i]]
  colnames(tmptd$dat)[1] <- ""
  charType <- "discrete"
  na.drop <- NULL
  charStates <- (c("0","1"))
  plot(tmptd$phy, cex=0.5, label.offset=0.1, edge.color="lightgreen", edge.width=3, tip.color="lightgreen", x.lim=c(0, 6))
  pal=colorRampPalette(c("green", "red"))
  nodelabels(pie=asr[[1]], piecol=pal(length(charStates)), cex=0.6, frame="circle")
  tiplabels(pch=21, bg=pal(length(charStates))[as.numeric(factor(tmptd$dat[,1], levels=charStates))]) 
  #legend("bottomleft", fill=pal(length(charStates)), legend=charStates)
  abline(v=TH-1.8, col="orange", lwd=3, lty=2)
}
dev.off()

pdf("../docs/images/Hormogonia.pdf", height=10, width=8)
par(bg="black")
trait <- "Hormogonia"
{
  i = which(names(tdcy$dat)==trait)
  resliksFull <- sapply(1:25, function(x) sapply(new.researchFull[[x]], function(y) y$lnLik))
  colnames(resliksFull) <- colnames(tdcy$dat)
  basicARDfn <- make.bisse.t(tdcyList[[i]]$phy, setNames(tdcyList[[i]]$dat[,1], tdcyList[[i]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
  pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
  asr <- list(t(asr.marginal(basicARDfn, pars, root=ROOT.FLAT)))
  tmptd <- tdcyList[[i]]
  colnames(tmptd$dat)[1] <- ""
  charType <- "discrete"
  na.drop <- NULL
  charStates <- (c("0","1"))
  plot(tmptd$phy, cex=0.5, label.offset=0.1, edge.color="lightgreen", edge.width=3, tip.color="lightgreen", x.lim=c(0, 6))
  pal=colorRampPalette(c("green", "red"))
  nodelabels(pie=asr[[1]], piecol=pal(length(charStates)), cex=0.6, frame="circle")
  tiplabels(pch=21, bg=pal(length(charStates))[as.numeric(factor(tmptd$dat[,1], levels=charStates))]) 
  #legend("bottomleft", fill=pal(length(charStates)), legend=charStates)
  abline(v=TH-1.8, col="orange", lwd=3, lty=2)
}
dev.off()


pdf("../docs/images/Nitrogen_fixation.pdf", height=10, width=8)
par(bg="black")
trait <- "Nitrogen_fixation"
{
  i = which(names(tdcy$dat)==trait)
  resliksFull <- sapply(1:25, function(x) sapply(new.researchFull[[x]], function(y) y$lnLik))
  colnames(resliksFull) <- colnames(tdcy$dat)
  basicARDfn <- make.bisse.t(tdcyList[[i]]$phy, setNames(tdcyList[[i]]$dat[,1], tdcyList[[i]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
  pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
  asr <- list(t(asr.marginal(basicARDfn, pars, root=ROOT.FLAT)))
  tmptd <- tdcyList[[i]]
  colnames(tmptd$dat)[1] <- ""
  charType <- "discrete"
  na.drop <- NULL
  charStates <- (c("0","1"))
  plot(tmptd$phy, cex=0.5, label.offset=0.1, edge.color="lightgreen", edge.width=3, tip.color="lightgreen", x.lim=c(0, 6))
  pal=colorRampPalette(c("green", "red"))
  nodelabels(pie=asr[[1]], piecol=pal(length(charStates)), cex=0.6, frame="circle")
  tiplabels(pch=21, bg=pal(length(charStates))[as.numeric(factor(tmptd$dat[,1], levels=charStates))]) 
  #legend("bottomleft", fill=pal(length(charStates)), legend=charStates)
  abline(v=TH-1.8, col="orange", lwd=3, lty=2)
}
dev.off()

pdf("../docs/images/Mats.pdf", height=10, width=8)
par(bg="black")
trait <- "Mats"
{
  i = which(names(tdcy$dat)==trait)
  resliksFull <- sapply(1:25, function(x) sapply(new.researchFull[[x]], function(y) y$lnLik))
  colnames(resliksFull) <- colnames(tdcy$dat)
  basicARDfn <- make.bisse.t(tdcyList[[i]]$phy, setNames(tdcyList[[i]]$dat[,1], tdcyList[[i]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
  pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
  asr <- list(t(asr.marginal(basicARDfn, pars, root=ROOT.FLAT)))
  tmptd <- tdcyList[[i]]
  colnames(tmptd$dat)[1] <- ""
  charType <- "discrete"
  na.drop <- NULL
  charStates <- (c("0","1"))
  plot(tmptd$phy, cex=0.5, label.offset=0.1, edge.color="lightgreen", edge.width=3, tip.color="lightgreen", x.lim=c(0, 6))
  pal=colorRampPalette(c("green", "red"))
  nodelabels(pie=asr[[1]], piecol=pal(length(charStates)), cex=0.6, frame="circle")
  tiplabels(pch=21, bg=pal(length(charStates))[as.numeric(factor(tmptd$dat[,1], levels=charStates))]) 
  #legend("bottomleft", fill=pal(length(charStates)), legend=charStates)
  abline(v=TH-0.1, col="orange", lwd=3, lty=2)
}
dev.off()

library(animation)
require(ggplot2)
require(reshape2)
require(grid)
#vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
col1 <- "darkolivegreen1"
names(col1) <- 1
col2 <- "lightcyan"
sequence = seq2
animfn <- function(X) {
  phyplot <- function(i) {
    par(bg="black", mar=c(5,8,0,1))
    basicARDfn <- make.bisse.t(tdcyList[[X]]$phy, setNames(tdcyList[[X]]$dat[,1], tdcyList[[X]]$phy$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
    pars <- researchFull[[X]][[i]]$par.full
    asr <- list(t(asr.marginal(basicARDfn, pars)))
    class(asr) <- c("asrArbor", "list")
    attributes(asr)$charType <- "discrete"
    attributes(asr)$td <- tdcyList[[X]]
    attributes(asr)$na.drop <- NULL
    attributes(asr)$charStates <- list(c("0","1"))
    plot(attributes(asr)$td$phy, edge.col=col1, edge.width=8, show.tip.label=FALSE, label.offset=0.5, x.lim=c(0,4.4))
    nodelabels(pie=asr[[1]], piecol=c("green", "red"), cex=1, frame="circle")
    tiplabels(pch=21, bg=pal(length(attributes(asr)$charStates))[as.numeric(factor(attributes(asr)$td$dat, levels=attributes(asr)$charStates))]) 
    #legend("bottomleft", fill=pal(length(charStates)), legend=charStates)
    #plot(asr, cex=0.5, label.offset=0.5, pal=colorRampPalette(c("green", "red")), edge.color=col1, edge.width=8, show.tip.label=FALSE, x.lim=c(0,4))
    #plot(ladderize(tdcy$phy), show.tip.label=FALSE, edge.color=col1, edge.width=8, tip.color=col1, cex=1, x.lim=c(0,4))
    tiplabels(pch=21, cex=3, col =1+as.numeric(as.factor(tdcyList[[X]]$dat[,1])), bg=1+as.numeric(as.factor(tdcy$dat$Motility)),adj=c(0.5, 0.5))
    abline(v=3.8-sequence[i], col="orange", lwd=15, lty=2)
  }
  profplot <- function(i){
    par(bg="black", mar=c(5,8,0,1))
    plot(resliksFull[,X], type="n", xlab="Gyr", ylab=expression(paste(Delta, "Ln lik", sep=" ")), xlim=c(3.8,-0.65), col.axis=col1, col.lab=col1, cex.axis=5, cex.lab=8, bty="n", xaxt="n", yaxt="n")
    lines(sequence[1:i], resliksFull[1:i,X], lwd=8, col=col1)
    abline(v=sequence[i], lty=2, lwd=15, col="orange")
  }
  circlerates <- function(rates, cex, center, dists, pal,adj,yadj, n=100,widthscale=10, log=TRUE, range=NULL,...){
    cols <- pal(n+2)
    rr <- log(rates)
    if(is.null(range)){range <- c(min(rr), max(rr))}
    rr <- rr-range[1]+1
    rr <- ifelse(rr < 1, 1, rr)
    rr <- ifelse(rr > (range[2]-range[1]), range[2], rr)
    rr <- round(rr/(range[2]-range[1]+1)*100,0)
    arrows(center[1]-1*dists[1] + adj,center[2]+yadj,  center[1]-1*dists[1]+dists[2]-adj, center[2]+yadj, col=cols[rr[2]], lwd=rr[2]/(n/widthscale+0.3), ...)
    arrows(center[1]-1*dists[1]+dists[2]-adj,center[2]-yadj, center[1]-1*dists[1] + adj, center[2]-yadj, col=cols[rr[4]], lwd=rr[4]/(n/widthscale+0.3), ...)
    arrows(center[1]+1*dists[1]-adj,center[2]-yadj,  center[1]+1*dists[1]-dists[2]+adj, center[2]-yadj, col=cols[rr[3]], lwd=rr[3]/(n/widthscale+0.3), ...)
    arrows(center[1]+1*dists[1]-dists[2]+adj,center[2]+yadj, center[1]+1*dists[1]-adj, center[2]+yadj, col=cols[rr[1]], lwd=rr[1]/(n/widthscale+0.3), ...)
    points(c(center[1]-1*dists[1], center[1]-1*dists[1]+dists[2], center[1]+dists[1], center[1]+dists[1]-dists[2]), rep(center[2],4), pch=21, cex=cex, bg=c("red", "green","green", "red"), col=c("red", "lightgreen","lightgreen", "red"))
    #text(c(center[1]-1*dists[1], center[1]-1*dists[1]+dists[2], center[1]+dists[1], center[1]+dists[1]-dists[2]), rep(center[2],4), labels=c("0", "1", "1", "0"), col="lightgreen", cex=cex)
  }
  for(i in c(1:37, rep(37,5))){
    #plot.new()
    layout(mat=matrix(c(1,1, 2), ncol=1, nrow=3))
    phyplot(i) 
    profplot(i)
    circlerates(researchFull[[X]][[i]]$par.full[c(6, 5, 9, 8)], 16, c(seq2[i], min(resliksFull[,X])+0.5),dists=c(0.8, 0.6), pal=cm.colors, adj=0.2, yadj=0.12, widthscale=25, range=c(-5, 4))
  }
}
animfn(6)

saveGIF(animfn(6), ani.height=1600, ani.width=1200, interval=0.15)



p <- ggplot(subset(df, Var2 <= 26), aes(Var1, value, colour=factor(Var2)))
p <- p + geom_line()
p <- p + geom_vline(xintercept=0.6, linetype=2)
p <- p + geom_vline(xintercept=1.8, linetype=2)
p <- p + xlim(c(3.8,0))
p <- p + theme_bw()
p <- p + xlab("Gyr")
p <- p + ylab("Support")
p <- p + facet_grid(grid2~grid1)
p


## Threshold analysis
require(Rphylip)
dtf <- tdcy$dat
dtf[is.na(dtf)] <- "?"
fit <- Rthreshml(tdcy$phy, dtf,types=rep("d", ncol(tdcy$dat)), path="~/Bioinf_programs/threshml/.", cleanup=FALSE)
save(fit, file="NSOAthreshmlFit.Rda")




## Check to see if the same ML's were found:
plot(apply(resliksFull, 2, max),cyFit.best$lnL)
plot(apply(resliksFull, 2, max)-cyFit.best$lnL)

## Do the same analysis for Archaea
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
registerDoMC(cores=8)
seq2 <- seq(3.7, 0.1, -0.1)
fns <- ARD.R2.fns
archresFull <- foreach(i=1:11) %dopar% {
  tmp <- list()
  for(j in 1:length(seq2)) {
    fn <- fns[[i]]
    ft <<- seq2[j]
    fn <- constrain2(fn, q10.tc~ft, extra=c("r1", "r2"))
    if(j==1){
      startx <- runif(4, 0.5, 2)*c(1, 1, IndFit.best.notime[i, 1], IndFit.best.notime[i, 2])
    } else {
      startx <- runif(4, 0.5, 2)*c(1, 1, IndFit.best.notime[i, 1], IndFit.best.notime[i, 2]) #c(tmp[[j-1]]$par)
    }
    tmp[[j]] <- find.mle(fn, x.init=startx, method="subplex")
  }
  tmp
}
archliksFull <- sapply(1:11, function(x) sapply(archresFull[[x]], function(y) y$lnLik))

## Do a more substantial search, using better starting points
registerDoMC(cores=11)
seq2 <- seq(3.7, 0.1, -0.1)
fns <- ARD.R2.fns
#archresFull2 <- foreach(i=1:11) %dopar% {
  tmp <- list()
  for(j in 1:length(seq2)) {
    fn <- fns[[i]]
    ft <<- seq2[j]
    fn <- constrain2(fn, q10.tc~ft, extra=c("r1", "r2"))
    if(j==1){
      startx <- unlist(IndFit.best[i,][,1:4])
    } else {
      startx <- unlist(IndFit.best[i,][,1:4])
    }
    ftmp <- find.mle(fn, x.init=startx, method="subplex")
    if(ftmp$lnL > archresFull[[i]][[j]]$lnL){
      tmp[[j]] <- ftmp
    } else {
      tmp[[j]] <- archresFull[[i]][[j]]
    }
  }
  tmp
}

archliksFull2 <- sapply(1:11, function(x) sapply(archresFull2[[x]], function(y) y$lnLik))

plot(seq2, apply(archliksFull, 1, sum), xlim=c(3.8, 0))
#lines(seq1,  apply(profliks, 1, sum))
##plot(seq2, apply(profliksFull, 1, sum))
#lines(seq2,  apply(profliksFull, 1, sum))
lines(seq2, apply(archliksFull, 1, sum), col="red", lwd=2)
lines(seq2, apply(archliksFull2, 1, sum), col="red", lwd=2)

rect(2.45, -100000, 1.85, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(0.663, -100000, 0.551, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)

png("cumlnLarch.png", width=1000)
## Plot for all traits
plot(seq2, apply(archliksFull2, 1, sum),type="n", xlim=c(3.8, -.8), ylim=c(0,15), bty="l", xlab="Billion years before present", ylab="Support")
#lines(seq1,  apply(profliks, 1, sum))
##plot(seq2, apply(profliksFull, 1, sum))
#lines(seq2,  apply(profliksFull, 1, sum))
aresfullN <- apply(archliksFull2, 2, function(x) x-min(x))
colnames(aresfullN) <- colnames(tdALL$dat)[1:11]
aresfullN <- aresfullN[,order(apply(aresfullN, 2, max))]
acumres <- sapply(1:nrow(aresfullN), function(x) cumsum(aresfullN[x, ]))
lapply(1:11, function(y) lines(seq2,acumres[y,], col=y, lwd=2))
rect(2.45, -100000, 1.85, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(0.663, -100000, 0.551, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
text(rep(0.1, 11), acumres[,ncol(acumres)], labels=colnames(aresfullN), cex=0.5, pos=4)
amax.lik <- seq2[which(acumres[11,]==max(acumres[11,]))]
abline(v=amax.lik)
abline(v=0.5)
dev.off()

plot(apply(archliksFull, 2, max), IndFit.best$lnL)
plot(apply(archliksFull, 2, max)- IndFit.best$lnL)
sum(IndFit.best$lnL)
sum(apply(archliksFull, 2, max))


## Combined archaea and cyanos
nn=37
plot(seq2, apply(archliksFull2, 1, sum),type="n", xlim=c(3.8, -.8), ylim=c(0,35), bty="l", xlab="Billion years before present", ylab="Support")
combresfullN <- cbind(resfullN, aresfullN)
combresfullN <- combresfullN[, order(apply(combresfullN, 2, max))]
ccumres <- sapply(1:nrow(combresfullN), function(x) cumsum(combresfullN[x, ]))
lapply(1:33, function(y) lines(seq2[1:nn],ccumres[y,1:nn], col=y, lwd=2))
rect(2.45, -100000, 1.85, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
rect(0.663, -100000, 0.551, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
text(rep(0.1, 33), ccumres[,ncol(ccumres)], labels=colnames(combresfullN), cex=0.5, pos=4)
cmax.lik <- seq2[which(ccumres[33,1:35]==max(ccumres[33,1:nn]))]
abline(v=0.5)
abline(v=1.7)

abline(v=seq2[apply(combresfullN, 2, function(x) which(x == max(x)))])

##Should simulate data under the phylogeny and reestimate, see if we get the same bimodality

#Subplex is a lot faster...
#system.time(for(i in 1:2) ttt <- tmpfn(22, 10, method="subplex"))
#ttt[[10]]$lnLik
#profliksFull[10, 22]
#profres[[22]][[14]]$lnLik
#ttt[[10]]$par
#profresFull[[1]][[10]]$par


#plot.likProf(mlik, tdGO1, c(1, fit.mlik$par[-1]), 2, c(3.8,0), c(-800, -700), n=100, lcol="red", add=TRUE)
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