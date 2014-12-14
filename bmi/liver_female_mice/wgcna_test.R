source("http://bioconductor.org/biocLite.R")
biocLite("impute")
install.packages("WGCNA")

library(WGCNA)

constDataDirectory <- "data"
constDataFile <- "LiverFemale3600.csv"
constTraitFile <- "ClinicalTraits.csv"
constPreprocessedFile <- "PreprocessedLiverData.RData"
constFemaleTOMDataFile <- "FemaleTOMData" # .RData is appended by the function itself
constNetworkConstructionDataFile <- "NetworkConstructionLiverData.RData"
constWorkingDirectory <- "~/Documents/projects/git-projects/r-projects/bmi/liver_female_mice"

setwd(constWorkingDirectory)
getwd()

options(stringsAsFactors=FALSE)

liverData <- read.csv(paste(constDataDirectory, constDataFile, sep="/"))
dim(liverData)
names(liverData)

# tse <- transposeSampleExpressions
tse <- as.data.frame(t(liverData[, -c(1:8)]))
colnames(tse) <- liverData$substanceBXH
rownames(tse) <- colnames(liverData)[-c(1:8)]

# gsg <- goodSampleSAndGenes
gsg <- goodSamplesGenes(tse)
gsg$allOK

if(!gsg$allOK){
  
  if(sum(!gsg$goodGenes) > 0){
    printFlush(paste("Removing genes: ", paste(names(tse)[!gsg$goodGenes], collapse=",")))
  }
  
  if(sum(!gsg$goodSamples) > 0){
    printFlush(paste("Removing samples: ", paste(rownames(tse)[!gsg$goodSamples], collapse=",")))
  }
  
  tse <- tse[gsg$goodSamples, gsg$goodGenes]
}

# st <- sampleTree
st <- flashClust(dist(tse), method="average")

sizeGrWindow(12, 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(st, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)
abline(h=15, col="red")

clust = cutreeStatic(st, cutHeight=15, minSize=10)
table(clust)

# gs <- goodSamples
gs = (clust == 1)
# transposeGoodSamplesExpression
tgse <- tse[gs, ]
nGenes <- ncol(tgse)
nSamples <- nrow(tgse)

# td <- traitData
td <- read.csv(paste(constDataDirectory, constTraitFile, sep="/"))
dim(td)
names(td)

allTraits <- td[, -c(31, 16)]
allTraits <- allTraits[, c(2, 11:36)]
dim(allTraits)
names(allTraits)

femaleSamples <- rownames(tgse)
traitRows <- match(femaleSamples, allTraits$Mice)
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]

collectGarbage()

# sampleTree2
st2 <- flashClust(dist(tgse), method="average")
traitColors <- numbers2colors(datTraits, signed=FALSE)
plotDendroAndColors(st2, traitColors, groupLabels=names(datTraits), main="Sample dendrogram and trait heatmap")

save(tgse, datTraits, file=paste(constDataDirectory, constPreprocessedFile, sep="/"))

# Part 2 - Network Construction and Module Detection
enableWGCNAThreads()
load(paste(constDataDirectory, constPreprocessedFile, sep="/"))

powers = c(c(1:10), seq(from=12, to=20, by=2))
sft = pickSoftThreshold(tgse, powerVector=powers, verbose=5)

sizeGrWindow(9, 5)
par(mfrow=c(1,2))
cex1=0.9

plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])* sft$fitIndices[, 2],
     xlab="Soft Threshold(power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n",
     main=paste("Scale Independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels=powers, cex=cex1, col="red")
abline(h=.90, col="red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab="Soft Threshold(power)",
     ylab="Mean Connectivity",
     type="n",
     main=paste("Mean Connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels=powers, cex=cex1, col="red")

net = blockwiseModules(tgse, power=6,
                       TOMType="unsigned",
                       minModuleSize=30,
                       reassignThreshold=0,
                       mergeCutHeight=.25,
                       numericLabels=TRUE,
                       pamRespectsDendro=FALSE,
                       saveTOMs=TRUE,
                       saveTOMFileBase=paste(constDataDirectory, constFemaleTOMDataFile, sep="/"),
                       verbose=3)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file=paste(constDataDirectory, constNetworkConstructionDataFile, sep="/"))

