## Analyze cyanobacteria dataset

## # Prepare the data table
require(aRbor)

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
      'Gas_vesicles', 'False_Branching', 'True_Branching', 'Fission_in_multiple_planes', 'Multiseriate_trichomes',
      'Baeocytes', 'Extracellular_sheath', 'Mucilage', 'celldiam_mean')
dat <- cyanodat
dat$Morphology[cyanodat$Morphology=="0&2"] <- 0
dat$Motility[cyanodat$Motility=="2"] <- 1
dat$Multiseriate_trichomes[cyanodat$Multiseriate_trichomes!=0] <- 1
dat$Mucilage[cyanodat$Mucilage!=0] <- 1
dat$Habit[cyanodat$Habit=="0&1"] <- 0
dat <- dat[,whichcol]

tree1 <- read.tree("../data/Tree1MrB.tre")
tree1$edge.length <- tree1$edge.length/1000
td1 <- make.treedata(tree1, dat, name_column=0)
#Check to see that all characters are binary (except for cell diameter)
all(sapply(td1$dat[,-ncol(td1$dat)], function(x) levels(factor(x))==c("0", "1")))

#Visualize the distribution of traits:
td1$dat <- td1$dat[,-ncol(td1$dat)]
asrDiscrete <- aceArbor(td1, charType="discrete", aceType="marginal", discreteModelType="ER", na.rm="bytrait")
attributes(asrDiscrete)$na.drop
par(mfrow=c(2,2))
plot(asrDiscrete, type="fan", show.tip.label=FALSE)

asrCont <- aceArbor(td1, charType="continuous")