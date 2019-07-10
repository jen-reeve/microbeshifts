setwd("~/GitHub/microbeshifts/analysis")
library(phytools)
library(phangorn)
library(tidyverse)
library(truncdist)
library(plotly)
file_path<-"shiftfunctions.R"
source(file_path)
file_path<-"cyanofunctions.R"
source(file_path)

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

# make tree for this
maxtree1<-maxCladeCred(cyanosTrees[[1]])

# discrete
sal_bin<-setNames(dat$Freshwater_habitat,rownames(dat))

mtrees<-make.simmap(maxtree1,sal_bin[names(sal_bin)%in%cyanos.species],type="discrete",model="ARD",nsim=100)

obj<-densityMap(mtrees,states=levels(as.factor(sal_bin))[2:1],plot=FALSE)
plot(obj,fsize=c(0,1),show.tip.label=FALSE)

pd<-summary(mtrees)
cols<-setNames(c("#7570B3","#D95F02"),levels(as.factor(sal_bin)))
plot(pd,fsize=0.,ftype="i",colors=cols,ylim=c(-2,Ntip(maxtree1)))
add.simmap.legend(colors=cols[2:1],prompt=FALSE,x=0,y=20,vertical=TRUE)

# simple optimum
dat_simp <- dat %>% 
  mutate(sal_simp = case_when(Freshwater_habitat == 1  ~ 0, 
                              Freshwater_habitat == 0 ~ 35,
                              TRUE ~ 0))
row.names(dat_simp) <- row.names(dat)

simpsal<-setNames(dat_simp$sal_simp,rownames(dat_simp))

simpmap<-contMap(maxtree1,simpsal[names(simpsal)%in%maxtree1$tip.label])
plot(simpmap,fsize=c(0,1))

fresh<-0
marine <-35
count<-100

p<- plot_ly(x=fresh,y=count,
            type="bar") %>%
  add_trace(x=marine,
            y=count,
            type="bar") %>%
  layout(barmode="overlay",
         yaxis=list(showticklabels=FALSE,
                    showgrid=FALSE),
         xaxis=list(showticklabels=F),
         showlegend=FALSE,
         colorway=c("#D95F02","#7570B3"))
p

# optimum from distribution
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
distsal <- setNames(dat_dist$sal_dist,rownames(dat_dist))

distmap<-contMap(maxtree1,distsal[names(distsal)%in%maxtree1$tip.label])
plot(distmap,fsize=c(0,1))

# showing distributions used
fresh<-rtrunc(100,a=0,b=Inf,spec="norm",mean=0,sd=5)
marine <- rtrunc(100,a=-Inf,b=Inf,spec="norm",mean=35,sd=4)

p<- plot_ly(x=fresh,type="histogram") %>%
  add_trace(x=marine,
            type="histogram") %>%
  layout(barmode="overlay",
         yaxis=list(showticklabels=FALSE,
                    showgrid=FALSE),
         xaxis=list(showticklabels=F),
         showlegend=FALSE,
         colorway=c("#D95F02","#7570B3"))
p