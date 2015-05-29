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

colnames(res) <- c("Lambda", "Chi-square Value", "p-value")
res <- data.frame("Trait"=rownames(res), res)
res$Trait <- as.character(res$Trait)
res$Trait[3] <- "Nonfreshwater_habitat"
write.csv(res, file="../output/phylogeneticSignalTable.csv")

prettyres <- res
prettyres$Lambda <- round(prettyres$Lambda, 3)
prettyres$Chi.square.Value <- round(prettyres$Chi.square.Value, 2)
prettyres$pvalue <- ifelse(prettyres$p.value > 0.05, "", "*")
prettyres$pvalue[prettyres$p.value < 0.01] <- "**"
prettyres$pvalue[prettyres$p.value < 0.001] <- "***"
results <- read.csv('../output/ResultsTable.csv')
fulltable <- merge(results, prettyres[,-3], by="Trait")
fulltable <- arrange(fulltable, desc(Shift.Cluster), desc(Max.Support), desc(Time.Estimate))
write.csv(fulltable, "../output/Table1.csv")
