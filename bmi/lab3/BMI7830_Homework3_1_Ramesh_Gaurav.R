###############
# Install Dependencies
###############

# source("http://bioconductor.org/biocLite.R")
# biocLite("impute")
# biocLite("GO.db")
# install.packages("WGCNA")

################
# Load required libraries
################

library(WGCNA)
library(hgu133plus2.db)
library(annotate)
library("GO.db")

###################################################
# Define constants(directory and file names for various intermediate in/output)
###################################################

constDataDirectory <- "data"
constDataFile <- "NSC_Lung_Cancer_Eset.tsv"
constAnnotationFile <- "GDS3627_SampleAnnotation.tsv"
constPreprocessedFile <- "Preprocessed_Data_Annotations.RData"
constTOMFile <- "NSC_Lung_Cancer_TOM"
constNetworkConstructionFile <- "NSC_Lung_Cancer_Network_Construction.RData"
constModTraitFile <- "ModTraitCorrelations.csv"
constGeneInfoFile <- "GeneInfo.csv"
constGeneEnrichmentFile <- "GeneEnrichmentData.csv"
constWorkingDirectory <- "~/Documents/projects/git-projects/r-projects/bmi/lab3"

setwd(constWorkingDirectory)
getwd()

options(stringsAsFactors=TRUE)
enableWGCNAThreads()

###################
# 1. Data Input and Cleaning
###################

inputData.original <- read.delim(paste(constDataDirectory, constDataFile, sep="/"), header=TRUE, sep="\t")

# re-initialize row names of the dataset and cleanup
rownames(inputData.original) <- inputData.original[, 1]

# probe ids are also read as columns.. remove them
# remove the first column(probe ids), to retain just the expression values
inputData.original <- inputData.original[, -1]

# transpose data frame for analysis: samples x genes(probes)
# WGCNA requires genes to be columns
inputData.trans <- as.data.frame(t(inputData.original))

# verbose = 3 => fully/strict verbose
goodSamplesAndGenes <- goodSamplesGenes(inputData.trans, verbose=3)
goodSamplesAndGenes$allOK

if(!goodSamplesAndGenes$allOK){
  
  if(sum(!goodSamplesAndGenes$goodGenes) > 0){
    printFlush(paste("Removing genes: ", paste(names(inputData.trans)[!goodSamplesAndGenes$goodGenes], collapse="\n")))
  }
  
  if(sum(!goodSamplesAndGenes$goodSamples) > 0){
    printFlush(paste("Removing samples: ", paste(rownames(inputData.trans)[!goodSamplesAndGenes$goodSamples], collapse="\n")))
  }
  
  # retain only good samples and good genes
  inputData.trans <- inputData.trans[goodSamplesAndGenes$goodSamples, goodSamplesAndGenes$goodGenes]
}

# prepare the label data frame(traits object)
traits <- read.delim(paste(constDataDirectory, constAnnotationFile, sep="/"), header=TRUE, sep="\t")
rownames(traits) <- traits[, 1]
traits <- subset(traits, select=c("Factors"))

# sample dendrogram to detect outliers(if any)
sampleTree <- flashClust(dist(inputData.trans), method="average")
sizeGrWindow(12, 9)
par(cex=.6)
par(mar=c(0, 4, 2, 0))
plot(sampleTree, cex.lab=1.5, cex.axis=1.5, cex.main=2)

# Relation of clinical trait to sample dendrogram
sampleTreeVsTraits <- flashClust(dist(inputData.trans), method="average")
colors <- numbers2colors(as.numeric(traits$Factors))
plotDendroAndColors(sampleTreeVsTraits, colors, groupLabels=names(traits))

# save cleaned input data and traits
save(inputData.trans, traits, file=paste(constDataDirectory, constPreprocessedFile, sep="/"))

load(paste(constDataDirectory, constPreprocessedFile, sep="/"))

################################
# 2. Network Construction and Module Detection
################################

# try out various soft thresholds
# powers <- seq(from=2, to=20, by=2)
# sft <- pickSoftThreshold(inputData.trans, powerVector=powers, verbose=5)

# above method gives 6 as a good soft threshold
softThreshold <- 6
bwnet <- blockwiseModules(inputData.trans, 
                          maxBlockSize=2000,
                          power=softThreshold, 
                          TOMType="unsigned", 
                          minModuleSize=30,
                          numericLabels=FALSE,
                          saveTOMs=TRUE,
                          saveTOMFileBase=paste(constDataDirectory, constTOMFile, sep="/"),
                          verbose=3)

# sizeGrWindow(12, 9)
# mergedColors = labels2colors(bwnet$colors)
# plotDendroAndColors(bwnet$dendrograms[[1]], 
#                     mergedColors[bwnet$blockGenes[[1]]], 
#                     dendroLabels=FALSE, 
#                     hang=.03, addGuide=TRUE, guideHang=.05)

bwModuleLabels <- bwnet$colors
bwModuleColors <- labels2colors(bwModuleLabels)
MEs <- bwnet$MEs
geneTree <- bwnet$dendrograms[[1]]
save(MEs, bwModuleLabels, bwModuleColors, geneTree, file=paste(constDataDirectory, constNetworkConstructionFile, sep="/"))

###############################################
# 3. Relating Modules to External Traits and Identifying Important Genes
###############################################

load(paste(constDataDirectory, constNetworkConstructionFile, sep="/"))

nGenes    = ncol(inputData.trans)
nSamples = nrow(inputData.trans)

MEs <- orderMEs(moduleEigengenes(inputData.trans, bwModuleColors)$eigengenes)

isFactor <- sapply(traits, is.factor)
factors2numbers <- traits
factors2numbers[isFactor] <- lapply(traits[isFactor], as.numeric)

moduleTraitCor <- cor(MEs, factors2numbers, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(12, 25)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep="")

dim(textMatrix) <- dim(moduleTraitCor)
par(mar=c(6, 8.5, 3, 3))
labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=names(traits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=FALSE,
               cex.text=0.2,
               cex.lab=0.5,
               zlim=c(-1, 1),
               main=paste("Module-Trait Relationships"))

traits <- as.data.frame(factors2numbers$Factors)
names(traits) <- "Tumor Type"

modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(inputData.trans, MEs, use="p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")

geneTraitSignificance <- as.data.frame(cor(inputData.trans, traits, use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) <- paste("GS.", names(traits), sep="")
names(GSPvalue) <- paste("p.GS", names(traits), sep="")

# sort the moduleTraitCorrelation based on absolute value of correlation
orderedModuleTraitCor <- abs(moduleTraitCor)
rowOrder <- order(orderedModuleTraitCor, decreasing=TRUE)
sortedCorrValues <- orderedModuleTraitCor[rowOrder]

# yellow, skyblue, thistle are found to have high correlation
sortedModValues  <- modNames[rowOrder]

#write the value to a file for future reference
cat(paste(sortedCorrValues, sortedModValues, sep=","), file=paste(constDataDirectory, constModTraitFile, sep="/"), sep="\n")

module="thistle"
column <- match(module, modNames)
moduleGenes <- bwModuleColors == module

sizeGrWindow(7, 7)
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab=paste("Module Membership in ", module, "module"),
                   ylab="Gene Significance for Traits",
                   main=paste("Module membership vs. gene significance\n"),
                   cex.main=1.2,
                   cex.lab=1.2, cex.axi= 1.2, col=module)

# most significant module was found to be yellow
#sModule <- significantModule
sModule <- "yellow"
# get the probe_ids for the genes in module yellow
sProbes <- names(inputData.trans)[bwModuleColors==sModule]
sProbeSymbols <- getSYMBOL(sProbes, "hgu133plus2.db")

# remove the probes with no gene id mapping
isNA <- is.na(sProbeSymbols)
sProbes <- sProbes[!isNA]
sProbeSymbols <- sProbeSymbols[!isNA]

# export relevant genes for enrichment using external tools, DAVID
sink(file="relevant_genes.txt")
cat(sProbeSymbols, sep="\n")
sink()

gene.symbols <- getSYMBOL(names(inputData.trans), "hgu133plus2.db")
gene.entrez_ids <- getEG(names(inputData.trans), "hgu133plus2.db")
gene.details <- data.frame(probe_ids=names(inputData.trans),
                           symbols=gene.symbols,
                           entrez_ids=gene.entrez_ids,
                           moduleColor=bwModuleColors,
                           geneTraitSignificance,
                           GSPvalue)

for (mod in 1:ncol(geneModuleMembership)){
  oldNames <- names(gene.details)
  gene.details <- data.frame(gene.details, geneModuleMembership[, rowOrder[mod]],
                             MMPvalue[, rowOrder[mod]])
  names(gene.details) <- c(oldNames, paste("MM.", modNames[rowOrder[mod]], sep=""),
                           paste("p.MM", modNames[rowOrder[mod]], sep=""))
}

geneOrder <- order(gene.details$moduleColor, -abs(gene.details$GS.Tumor.Type))
geneInfo <- gene.details[geneOrder, ]

write.csv(geneInfo, file=paste(constDataDirectory, constGeneInfoFile, sep="/"), row.names=FALSE)

##############################
# 4. Gene Ontology and Functional Enrichment
##############################

# enrichment using GO packages
GOEnr <- GOenrichmentAnalysis(bwModuleColors, geneInfo$entrez_ids, organism="human", nBestP=10)
tab <- GOEnr$bestPTerms[[4]]$enrichment

write.table(tab, file=paste(constDataDirectory, constGeneEnrichmentFile, sep="/"), sep=",", quote=TRUE, row.names=FALSE)

###############################
# 5. Network Visualization
###############################

tab <- read.table(paste(constDataDirectory, constGeneEnrichmentFile, sep="/"), header=TRUE, sep=",") 

intGenes <- tab$module == "yellow" | tab$module == "skyblue" | tab$module == "thistle"
intGeneEnrichment <- tab[intGenes, ]

# load TOM file
# load(paste(constDataDirectory, paste(constTOMFile, "block.1.RData", sep="-"), sep="/"))
nTOMFiles <- length(bwnet$TOMFiles)

# for all the ToM files saved, plot TOMplot
for(i in 1:nTOMFiles){
  TOMFileName <- bwnet$TOMFiles[i]
  TOMFileName
  load(TOMFileName)
  
  dissTOM <- 1 - TOM
  dissTOM <- as.matrix(dissTOM)
  plotTOM <- dissTOM^7
  
  png(filename=paste(TOMFileName, "img.png", sep="-"))
  TOMplot(plotTOM, bwnet$dendrograms[[i]], bwnet$colors[bwnet$blockGenes[[i]]], main=paste("Network Heatmap, Block", i, sep=" "))
  dev.off()
}

MET = orderMEs(cbind(MEs, factors2numbers))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)