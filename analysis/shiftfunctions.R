## Constants
#grO2event <- max(branching.times(tr))-2.32
#secO2event <- max(branching.times(tr))-(0.663+0.551)/2

##New constrain function that constrains to a multiplier

#constrain(x, lambda0~lambda.est, lambda1~lambda.est, mu0~0, mu1~0, q01.y0~ r, q10.y0~ r, q01.tc ~ q10.tc, extra="r")

#f <- fn1
#names <- argnames(f)
#constrain.parse <- diversitree:::constrain.parse
#extra="r"

#dotgetter <- function(...){
#  formulae <- list(...)
#}

constrain2 <- function (f, ..., formulae = NULL, names = argnames(f), extra = NULL) {
  if (is.constrained(f)) {
    formulae <- c(attr(f, "formulae"), formulae)
    f <- attr(f, "func")
  }
  formulae <- c(formulae, list(...))
  names.lhs <- names.rhs <- names
  rels <- list()
  for (formula in formulae) {
    res <- constrain.parse2(formula, names.lhs, names.rhs, extra)
    if (attr(res, "lhs.is.target")) {
      i <- try(which(sapply(rels, function(x) identical(x, 
                                                        res[[1]]))), silent = TRUE)
      if (inherits(i, "try-error")) 
        stop(sprintf("Error parsing constraint with %s on lhs", 
                     as.character(res[[1]])))
      rels[i] <- res[[2]]
      lhs.txt <- as.character(res[[1]])
      if (any(sapply(rels, function(x) lhs.txt %in% all.vars(x)))) 
        stop(sprintf("lhs (%s) is in an expression and can't be constrained", 
                     lhs.txt))
    }
    names.lhs <- setdiff(names.lhs, unlist(lapply(res, all.vars)))
    names.rhs <- setdiff(names.rhs, as.character(res[[1]]))
    rels <- c(rels, structure(res[2], names = as.character(res[[1]])))
  }
  final <- c(extra, names.rhs)
  npar <- length(final)
  free <- setdiff(names.rhs, names(rels))
  free.i <- match(free, names)
  free.j <- match(free, final)
  target.i <- match(names(rels), names)
  pars.out <- rep(NA, length(names))
  names(pars.out) <- names
  g <- function(pars, ..., pars.only = FALSE) {
    if (length(pars) != npar) 
      stop(sprintf("Incorrect parameter length: expected %d, got %d", 
                   npar, length(pars)))
    pars.out[free.i] <- pars[free.j]
    e <- structure(as.list(pars), names = final)
    pars.out[target.i] <- unlist(lapply(rels, eval, e))
    if (pars.only) 
      pars.out
    else f(pars.out, ...)
  }
  class(g) <- c("constrained", class(f))
  attr(g, "argnames") <- final
  attr(g, "formulae") <- formulae
  attr(g, "extra") <- extra
  attr(g, "func") <- f
  g
}

#pars <- c(10, 2, 3, 2.5)
#g(pars)         
#formula <- q01.y0 ~ q01.y1 * r
#names.lhs <- argnames(f)
#names.rhs <- argnames(f)
#extra="r"
          
#formula2 <- q01.y0 ~ r
#constrain.parse2(formula2, names.lhs, names.rhs, extra)

constrain.parse2 <- function (formula, names.lhs, names.rhs, extra = NULL) {
  formula <- as.formula(formula)
  if (length(formula) != 3L) 
    stop("Invalid formula")
  lhs <- formula[[2]]
  rhs <- formula[[3]]
  if (!is.name(lhs)) 
    stop("Invalid target on LHS of formula")
  lhs.is.target <- is.na(match(as.character(lhs), names.lhs))
  if (is.language(rhs)) {
    vars <- all.vars(rhs)
    ok <- (all(vars %in% names.rhs) || length(vars) == 1 && vars %in% extra)
    if (!ok && length(vars) == 1) {
      e <- parent.frame()
      if (exists(vars, e)) {
        rhs <- get(vars, e)
        ok <- TRUE
      }
    }
    if (!ok && length(vars) > 1){
      ok <- TRUE
    }
    if (!ok) 
      stop("Invalid RHS of formula:\n\t", as.character(rhs))
    if (as.character(lhs) %in% vars) 
      stop("LHS cannot appear in RHS")
  }
  else if (!is.numeric(rhs)) {
    stop("RHS must be expression, variable or number")
  }
  res <- list(lhs, rhs)
  attr(res, "lhs.is.target") <- lhs.is.target
  res
}

## Get likelihood profile from a likelihood function
get.likProf <- function(likfn, fitpars, var, range, n){
  parrange <- seq(range[1], range[2], length.out=n)
  sapply(parrange, function(x){pars <- fitpars; pars[var] <- x; likfn(pars)})
}

plot.likProf <- function(likfn, td, fitpars, var, range, ylim=NULL, n=100, add=FALSE, lcol = "black", ...){
  TH <- max(branching.times(td$phy))
  LP <- get.likProf(likfn, fitpars, var, range, n)
  if(is.null(ylim)){
    ylim = c(max(LP) - 20, max(LP)+1)
    ylim.min <- ylim[1]
  } else {ylim.min=ylim[1]}
  if(!add){
    plot(seq(TH, 0,length.out=n), mlik.prof, type="n", ylim=ylim, xlim=c(TH, 0), xaxt="n",ylab="Log Likelihood", xlab="Time before present (by)",...)
    axis(1, at=c(3.8, 3 ,2 ,1, 0), labels=c(3.8, 3, 2, 1, 0))
    text(2.32, ylim.min, "Great Oxygenation Event",cex=0.7)
    text(0.607, ylim.min, "Secondary Oxygenation Event",cex=0.7)
    rect(2.45, -100000, 1.85, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
    rect(0.663, -100000, 0.551, 10000, col=rgb(100,100,100, maxColorValue=255,alpha=25), border=NA)
  }
  lines(seq(TH, 0, length.out=n), LP, col=lcol)
}

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

binarify <- function(x){
  x[x=="0&1"] <- 0
  x[x %in% 2:10] <- 1
  as.numeric(as.character(x))
}

start.gen <- function(fn, td){
  TH <- max(branching.times(td$phy))
  arg <- argnames(fn)
  rs <- grep("r", arg)
  qs <- grep("q", arg)
  ts <- grep("t", arg)
  qs <- qs[!(qs %in% ts)]
  pars <- rep(0, length(arg))
  pars[rs] <- exp(runif(length(rs), -5, 5))
  pars[qs] <- exp(runif(length(qs), -5, 5))
  pars[ts] <- runif(length(ts), 0, TH)
  pars
}
lower.gen <- function(fn, td, ln=FALSE){
  TH <- max(branching.times(td$phy))
  arg <- argnames(fn)
  rs <- grep("r", arg)
  qs <- grep("q", arg)
  ts <- grep("t", arg)
  qs <- qs[!(qs %in% ts)]
  pars <- rep(0, length(arg))
  pars[rs] <- 0
  pars[qs] <- 0
  pars[ts] <- 0
  pars
}
upper.gen <- function(fn, td, ln=FALSE){
  TH <- max(branching.times(td$phy))
  arg <- argnames(fn)
  rs <- grep("r", arg)
  qs <- grep("q", arg)
  ts <- grep("t", arg)
  qs <- qs[!(qs %in% ts)]
  pars <- rep(0, length(arg))
  pars[rs] <- 10000
  pars[qs] <- 10000
  pars[ts] <- TH
  pars
}
defactor <- function(x){
  as.numeric(as.character(x))
}

start.genGO1 <- function(fn, td){
  TH <- max(branching.times(td$phy))
  arg <- argnames(fn)
  rs <- grep("r", arg)
  qs <- grep("q", arg)
  ts <- grep("t", arg)
  qs <- qs[!(qs %in% ts)]
  pars <- rep(0, length(arg))
  pars[rs] <- exp(runif(length(rs), -5, 5))
  pars[qs] <- exp(runif(length(qs), -5, 5))
  pars[ts] <- runif(length(ts), 1.85, 2.45)
  pars
}
start.genGO2 <- function(fn, td){
  TH <- max(branching.times(td$phy))
  arg <- argnames(fn)
  rs <- grep("r", arg)
  qs <- grep("q", arg)
  ts <- grep("t", arg)
  qs <- qs[!(qs %in% ts)]
  pars <- rep(0, length(arg))
  pars[rs] <- exp(runif(length(rs), -5, 5))
  pars[qs] <- exp(runif(length(qs), -5, 5))
  pars[ts] <- runif(length(ts), 0.551, 0.663)
  pars
}


profile.mle <- function(fns, sequence, starts, previous=NULL, extra=c("r1", "r2"), method="subplex", cores=1){
  registerDoMC(cores=cores)
  res <- foreach(i=1:length(fns)) %dopar% {
    tmp <- list()
    for(j in 1:length(sequence)) {
      fn <- fns[[i]]
      ft <<- sequence[j]
      fn <- constrain2(fn, q10.tc~ft, extra=extra)
      startx <- starts[[j]]
      if(is.null(previous)){
        tmp[[j]] <- find.mle(fn, x.init=startx, method="optim")
      } else{
        ftmp <- find.mle(fn, x.init=startx, method="optim")
        if(ftmp$lnL > previous[[i]][[j]]$lnL){
          tmp[[j]] <- ftmp
          print(paste("Found a better lnL for trait", i, "at", j, "bybp", sep=" "))
        } else {
          tmp[[j]] <- previous[[i]][[j]]
        }
      }
    }
    tmp
  }
  return(res)
}


