# DOES THE DATA HAVE PHYLOGENY??!?

# get tree
source("./dataPrepOnly.R")
plot(tree)

# get data
cyanData<-clean.data()

# mash tree and data together using aRbor stuff
td <- make.treedata(tree, cyanData)

# now loop through all variables and test for signal
allTheVariables<-colnames(td$dat)

signalResults<-list()
for(i in 1:length(allTheVariables)) {
  theTrait<-allTheVariables[i]
  selectedData<-select_(td,theTrait)
  filteredData<-filter(selectedData, !is.na(selectedData$dat))

  signalResults[[i]]<-physigArbor(filteredData, charType="discrete")
}

# create table of results
res<-matrix(nrow=length(allTheVariables), ncol=3)
rownames(res)<-allTheVariables
for(i in 1:length(allTheVariables)) {
 res[i,]<-c(signalResults[[i]][[1]]$lambdaValue, signalResults[[i]][[1]]$chisqTestStat, signalResults[[i]][[1]]$chisqPVal) 
}  

res
