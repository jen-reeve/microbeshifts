## Cyanofunctions
## Cleaning and preparing the data
## This function should likely receive some attention from Carrine!
clean.data <- function(){
  ## Read in character matrix and cell diameter data
  cyanodat <- read.table("../data/charactermatrix_5_2015.txt")
  celldiam <- read.table("../data/celldiameter.msq")
  rownames(celldiam) <- celldiam[,1]
  celldiam <- celldiam[,-1]
  rownames(cyanodat) <- cyanodat[,1]
  cyanodat <- cyanodat[,-1]
  ## Read in the characters and states metadata
  states <- readLines("../data/charactersandstates_5_2015.txt")
  ## Clean the text so that they can be used as column names
  states[grep(".", states)]
  heads <- strsplit(states[grep(".", states)], "\\. ", perl=TRUE)
  colnames(cyanodat) <- sapply(gsub("-", "", gsub(" ", "_", as.data.frame(do.call(rbind, heads[sapply(heads, length)==2]))[,2])), function(x) substr(x, 1, nchar(x)-2))
  cyanodat$celldiam_min <- celldiam[,1]
  cyanodat$celldiam_mean <- celldiam[,2]
  cyanodat$celldiam_max <- celldiam[,3]
  cyanodat[cyanodat=="?"] <- NA
  ## Determine what set of traits will be kept
  whichcol <- c('Thermophilic', 'Freshwater_habitat','temp', 'Akinetes', 'Heterocysts', 'Nitrogen_fixation', 'Morphology',
              'Habit', 'Freeliving', 'Mats', 'Epi/Endolithic',  'Epiphytic', 'Periphytic', 'Motility', 'Hormogonia', 
              'Gas_vesicles', 'False_Branching', 'True_Branching', 'Fission_in_multiple_planes', 'Uniseriate_trichome', 
              'Multiseriate_trichomes', 'Baeocytes', 'Extracellular_sheath', 'Mucilage', 'celldiam_mean')
  dat <- cyanodat
  ## Decisions on how to group ambiguous or multistate characters
  dat$Morphology[cyanodat$Morphology=="0&2"] <- 0
  dat$Motility[cyanodat$Motility=="2"] <- 1
  dat$Multiseriate_trichomes[cyanodat$Multiseriate_trichomes!=0] <- 1
  dat$Mucilage[cyanodat$Mucilage!=0] <- 1
  dat$Habit[cyanodat$Habit=="0&1"] <- 1
  dat$temp[cyanodat$temp=="0&1"] <- NA
  dat <- dat[,whichcol]
  dat$Pelagic <- as.numeric(dat$temp==1 & dat$Habit==0)
  dat$celldiam_mean <- as.numeric(as.numeric(as.character(cyanodat$celldiam_mean))>=3.5)
  names(dat)[which(names(dat)=="temp")] <- "Nonfreshwater_habitat"
  # change Epi/Endolithic to remove slash
  colnames(dat)[colnames(dat)=="Epi/Endolithic"]<-"Epi_Endolithic"
  return(dat)
}



## 
make.bisse.fns <- function(tree, dat){
  ## set branch lengths to by
  #tree$edge.length <- tree$edge.length/1000
  
  ## Combine tree 
  tdcy <- make.treedata(tree, dat, name_column=0)
  rownames(tdcy$dat) <- tdcy$phy$tip.label
  colnames(tdcy$dat) <- gsub("/", "_", colnames(tdcy$dat), fixed=TRUE)
  nc <- ncol(tdcy$dat)
  tdcyList <- lapply(1:nc, function(x) select(tdcy, x))
  tdcyList <- lapply(1:nc, function(x) filter_(tdcyList[[x]], paste("!is.na(",names(tdcyList[[x]]$dat),")", sep="")))
  
  ## Set global birth-death parameters
  bd.lik <- make.bd(tdcy$phy)
  bd.est <- find.mle(bd.lik, x.init=c(5,0))
  lambda <<- bd.est$par[1]
  mu <<- bd.est$par[2]
  
  ## Make the functions
  bisse.fns <- lapply(1:nc, function(x) make.bisse.t(tdcyList[[x]]$phy, setNames(tdcyList[[x]]$dat[[1]], attributes(tdcyList[[x]])$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))) 
  notime.fns <- lapply(1:nc, function(x) make.bisse(tdcyList[[x]]$phy, setNames(tdcyList[[x]]$dat[[1]], attributes(tdcyList[[x]])$tip.label)))
  ARD.notime.fns <- lapply(notime.fns, function(x) constrain(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu))
  ARD.R1.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu, q10.y1~r*q10.y0, q01.y1~r*q01.y0, q01.tc ~ q10.tc, extra=c('r')))
  ARD.R2.fns <- lapply(bisse.fns, function(x) constrain2(x, lambda0~lambda, lambda1~lambda, mu0~mu, mu1~mu, q10.y1~r1*q10.y0, q01.y1~r2*q01.y0, q01.tc ~ q10.tc, extra=c('r1', 'r2')))
  
  fns <- list(bisse=bisse.fns, notime=notime.fns, ARD.notime=ARD.notime.fns, ARD.R1=ARD.R1.fns, ARD.R2=ARD.R2.fns, tdList=tdcyList)
  return(fns)
  
}

fullProfiles <- function(n, fns, res=NULL){
  
}

fitFns <- function(n, fns, tds, res=NULL){
  nc <- length(fns)
  for(j in 1:n) {
    fits <- foreach(i=1:nc) %dopar% {
      find.mle(fns[[i]], x.init=start.gen(fns[[i]], tds[[i]]), method="optim", lower=lower.gen(fns[[i]], tds[[i]]), upper=upper.gen(fns[[i]], tds[[i]]) )
    }
    ## Create a summary table
    parests <- do.call(rbind, lapply(fits, function(x) as.data.frame(matrix(c(x$par, x$lnLik, x$message), nrow=1))))
    colnames(parests) <- c(argnames(fns[[1]]), "lnL", "message")
    ## Save only the best-fitting indpendent runs
    if(is.null(res)){
      res <- parests
    } else {
      replace <- which(defactor(parests$lnL) > defactor(res$lnL))
      if(length(replace) > 0){
        res[-ncol(res)] <- apply(res[-ncol(res)], 2, defactor)
        parests[-ncol(parests)] <- apply(parests[-ncol(parests)], 2, defactor)
        res[replace, ] <- parests[replace, ]
        #cyFit.best[-ncol(cyFit.best)] <- apply(cyFit.best[-ncol(cyFit.best)], 2, defactor)
        #parests[-ncol(parests)] <- apply(parests[-ncol(parests)], 2, defactor)
        #cyFit.best[replace, ] <- parests[replace, ]
      }
      print(replace)
      
    }
  }
  return(res)
}

profiles <- function(n, fns, tds, starts, res=NULL, seq=seq(0.1, 3.7, 0.2), cores=1, start4=FALSE){
  ## Start only over seq1, add seq2 if it looks productive
  nc <- length(fns)
  if(is.null(res)){
      res <- list()
      seqFns <- lapply(fns, function(x) lapply(seq, function(y) {ft <<- y; constrain2(x, q10.tc~ft, extra=c("r1", "r2"))}))
      for(i in 1:nc){
        tmpfns <- seqFns[[i]]
        tmp <-  mclapply(1:length(seq), function(j) {
          fn <- tmpfns[[j]];
          ft <<- seq[j];
          #fn <- constrain2(fn, q10.tc~ft, extra=c("r1", "r2"))
          if(start4){
            startx <- starts[[i]][,1:4];
            startx[which(startx[,1:2]>=10000)] <- 10000
            startx[which(startx[,1:2]<=0.01)] <- 0.01+0.0000001
            startx[which(startx[,3:4]>=0.5*100)] <- 0.5*100-0.1
            startx[2*nrow(startx)+which(startx[,3:4]<=0.01/100)] <- 0.01/100+0.0000000001
          } else{
            startx <- cbind(1, 1, starts[i, 1], starts[i, 2]);
            startx[which(startx[,3:4]>=0.5*100)] <- 0.5*100-0.1
            startx[2*nrow(startx)+which(startx[,3:4]<=0.01/100)] <- 0.01/100+0.0000000001
          }
          if(start4) {
            optimx(startx[j,], fn, method=c("L-BFGS-B", "nlminb", "spg", "bobyqa", "hjkb"), lower=c(0.001,0.001,0.01/100, 0.01/100), upper=c(1000000, 1000000, 0.5*100, 0.5*100), control=list(maximize=TRUE))
          } else {
            optimx(as.vector(startx), fn, method=c("L-BFGS-B", "nlminb", "spg", "bobyqa","nmkb", "hjkb"), lower=c(0.001,0.001,0.01/100, 0.01/100), upper=c(1000000, 1000000, 0.5*100, 0.5*100), control=list(maximize=TRUE))
          }
          #find.mle(fn, x.init=startx, method="subplex")
        }, mc.preschedule=FALSE, mc.cores=cores)
        res[[i]] <- tmp
        rm(tmp)     
        gc()
      }
  }# else {
  #  for(i in 1:nc){
  #    tmp <- foreach(j=1:length(seq)) %dopar% {
  #      fn <- fns[[i]]
  #      ft <<- seq[j]
  #      fn <- constrain2(fn, q10.tc~ft, extra=c("r1", "r2"))
  #      if(j==1){
  #        startx <- runif(4, 0.5, 2)*c(1, 1, starts[i, 1], starts[i, 2])
  #      } else {
  #        startx <- runif(4, 0.5, 2)*c(1, 1, starts[i, 1], starts[i, 2]) #c(tmp[[j-1]]$par)
  #      }
  #      find.mle(fn, x.init=startx, method="subplex")
  #    }
  #    replace <- which(sapply(tmp, function(x) x$lnLik) > sapply(res[[i]], function(x) x$lnLik))
  #    if(length(replace)>0){
  #      res[[i]][replace] <- tmp[replace]
  #    }
  #  }
  #}
  res
}

smooth.profiles <- function(fns, tds, prof, seq, cores=1){
  nc <- length(fns)
  seqFns <- lapply(fns, function(x) lapply(seq, function(y) {ft <<- y; constrain2(x, q10.tc~ft, extra=c("r1", "r2"))}))
  nprof <- lapply(prof, function(x) lapply(x, function(y) y[which(y$value==max(y$value)), ]))
  for(i in 2:length(seq)){
    tmpfns <- lapply(seqFns, function(x) x[[i]])
    tmp <-  mclapply(1:nc, function(j) {
    for(j in 1:nc){
    fn <- function(x) tmpfns[[j]](exp(x));
    startx <- log(unlist(nprof[[j]][[i-1]][,1:4])+1e-12)
    #snlminb(startx, function(x) -1*fn(x), lower=c(0,0,0,0))
    x <- optimx(startx, fn, method="nlminb", control=list(maximize=TRUE))
    tmp[[j]] <- x
    }
    #find.mle(fn, x.init=startx, method="subplex")
    }, mc.preschedule=TRUE, mc.cores=cores)
    replace <- which(sapply(nprof, function(x) x[[i]]$value) < sapply(tmp, function(x) x$value))
    if(length(replace) > 0){
      res[[replace]][[i]] <- tmp[[replace]]
    }
    rm(tmp)     
    gc()
  }
  return(res)
}

read.trees.from.mesquite <- function(filename){
  text <- scan(filename, what="character")
  trids <- grep("coccus", text)
  trees <- text[trids]
  trees <- lapply(trees, function(x) gsub("\\","", x ,fixed=TRUE))
  trees <- lapply(trees, function(x) read.tree(text = x))
  return(trees)
}


summarizeProfile <- function(prof){
  dlnL <- max(prof[,1]-min(prof[37,1]))
  maxT <- seq(0.1,3.7,0.1)[which(prof[,1]==max(prof[,1]))]
  slnL <- prof[,1]-min(prof[37,1])[1]
  list(dlnL=dlnL, maxT=maxT, slnL=slnL)
}

makeTransparent <- function (someColor, alpha = 100) {
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata) {
    rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], 
        alpha = alpha, maxColorValue = 255)
  })
}

dlnLColors <- function(x, col, max=6, n=255){
  ind <- floor(x/max*n)
  sapply(ind, function(x) makeTransparent(col, alpha=x))
}

## Process a list of fit objects for each trait
processFits <- function(fits, names="all"){
  if(names=="all"){
    traits <- 1:length(fits)
  } else {
    traits <- names
  }
  results <- fits[traits]
  res <- lapply(results, function(trait) do.call(rbind, lapply(1:length(trait), function(time) c(trait[[time]]$lnLik, trait[[time]]$par))))
  cums <- apply(sapply(res, function(x) x[,1]), 1, sum)
  cums <- cums-min(cums)
  sumRes <- lapply(res, summarizeProfile)
  names(sumRes) <- gsub("/", "_", names(sumRes), fixed=TRUE)
  starts <- lapply(res, function(x) x[,-1])
  names(starts) <- gsub("/", "_", names(starts), fixed=TRUE)
  out <- list(profile=cums,  traitSums=sumRes, estimates=starts)
  return(out)
}

##Process a list of optimx results for each trait
processProfiles <- function(profiles, names="all"){
  if(names=="all"){
    traits <- 1:length(profiles)
  } else {
    traits <- names
  }
  results <- profiles[traits]
  lnLs <- lapply(results, function(x) sapply(x, function(y) max(y$value, na.rm=TRUE)))
  slnLs <- lapply(lnLs, function(x) x-min(x))
  maxT <- lapply(slnLs, function(x) which(x==max(x)))
  dlnL <- lapply(slnLs, function(x) max(x))
  traitSums <- lapply(1:length(traits), function(x) list(dlnL=dlnL[[x]], maxT=maxT[[x]], slnL=slnLs[[x]]))
  cums <- sapply(1:length(slnLs[[1]]), function(x) sum(sapply(slnLs, function(y) y[x])))
   
}

plot.asrR2 <- function(pars, tree, dat, fn, cols=c("gray10", "gray90"), ...){
  fns <- make.bisse.fns(tree, dat)
  fn <- fns$ARD.R2[[1]]
  td <- fns$tdList[[1]]
  #basicARDfn <- make.bisse.t(tree, setNames(dat, tree$tip.label), functions=c(rep("constant.t",4), rep("stepf.t", 2)))
  #pars <- researchFull[[i]][[which(resliksFull[,i]==max(resliksFull[,i]))]]$par.full
  asr <- list(t(asr.marginal(fn, pars)))
  TH <- max(branching.times(tree))
  class(asr) <- c("asrArbor", "list")
  attributes(asr)$charType <- "discrete"
  attributes(asr)$aceType <- "poopboobs"
  attributes(asr)$td <- td
  attributes(asr)$na.drop <- NULL
  attributes(asr)$charStates <- list(c("0","1"))
  plot(asr, cex=0.5, label.offset=0.1, pal=colorRampPalette(cols), ...)
  abline(v=TH-pars[5], lwd=2, lty=2)
  #draw.circle(0, 0, TH-pars[5], lwd=1, border=makeTransparent('gray10', alpha=100))
  
}
