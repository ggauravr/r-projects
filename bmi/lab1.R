# install the necessary bioconductor packages/libraries
source("http://bioconductor.org/biocLite.R")
biocLite(pkgs=c("RColorBrewer", "limma", "affy", "gcrma", 
                "Biobase", "affyPLM", "affydata", "GEOquery", "simpleaffy"))
biocLite("hgu133plus2.db")

# load the required packages
library(GEOquery)
library(simpleaffy)
library(Biobase)
library(RColorBrewer)
library(affyPLM)
library(limma)
library(hgu133plus2.db)
library(annotate)
library(ggplot2)

# define constants needed for the script
constDataset <- "GSE10245"
constDataDirectory <- "data"
constEsetFile <- "gcrmaExpressionSet.txt"
constWorkingDirectory <- "~/Documents/projects/git-projects/r-projects/bmi"

# set the working directory
# directory where the required files will be downloaded
setwd(constWorkingDirectory)

# download the dataset, extract the .tar file
# unzip the resulting gzipped files, to get the .cel files
getGEOSuppFiles(constDataset)
untar(paste(constDataset, "/", constDataset, "_RAW.tar", sep=""), exdir=constDataDirectory)
sapply(paste(constDataDirectory, list.files(constDataDirectory, pattern = "gz$") , sep="/"), gunzip)
  
# generate AffyBatch object from CEL and corresponding Annotation files
# normalize the expression data using GC-RMA method and store the result
# ExpressionSet object into a file for ease of future use

celAffyBatch <- read.affy(covdesc="GDS3627_SampleAnnotation.txt", path=constDataDirectory)
# gcrmaEset <- gcrma(celAffyBatch)
# write.exprs(gcrmaEset, file=constEsetFile)

# load expressionSet object from the file saved above
gcrmaEset <- readExpressionSet(constEsetFile)
pData(gcrmaEset) <- pData(celAffyBatch)
annotation(gcrmaEset) <- annotation(celAffyBatch)

# set colour palette
colors <- brewer.pal(8, "Set1")

# boxplot parameters:
#   las = 3 indicates vertical naming in the x-axis
# 	cex.axis adjusts the font of axis names, in terms of fraction of the default size

# boxplot of unnormalized intensity values
boxplot(celAffyBatch, col=colors, las=3, cex.axis=0.6)
# plot a boxplot of normalized intensity values
boxplot(gcrmaEset, col=colors, las=3, cex.axis=0.6)

# histogram of unnormalized values
hist(celAffyBatch, col=colors)
# histogram of normalized values
hist(gcrmaEset, col=colors)

celAffyBatch.filtered <- nsFilter(gcrmaEset, require.entrez=FALSE, remove.dupEntrez=FALSE)
# log shows the various features excluded and the reason of exclusion
# here, the only exclusions are from features where IQR < 0.5(default)
# and the Affymetrix Quality Control probesets, that is the probesets whose names
# begin with AFFY
# SEE: ?nsFilter to understand the various arguments and the output
celAffyBatch.filtered$filter.log
gcrmaEset.filtered <- celAffyBatch.filtered$eset

# prepare the design matrix and rename columns of the design as needed
samples <- factor(gcrmaEset$Factors)
design <- model.matrix(~0 + samples)
colnames(design) <- levels(samples)

# produce a linear model fit, for every row of the expression data
# try dim(fit) to check the rows(number of linear fits) vs cols(number of coefficients)
fit <- lmFit(gcrmaEset.filtered, design)
cont.matrix <- makeContrasts(ACvsSCC=adenocarcinoma - squamous_cell_carcinoma, levels=design)
eBayesFit <- eBayes(contrasts.fit(fit, cont.matrix))

# print the number of probesets with |logFC| > 2, 3, 4
# results: 481, 150 and 67 respectively
lapply(2:4, function(x) nrow(topTable(eBayesFit, coef=1, number=Inf, lfc=x)) )

# get the probeset_ids and the relevant statistics
# map the probeset_ids to their corresponding gene symbols, using the annotation db
# combine/map probeset_ids and their statistics to gene symbols
# write the results into a file
gene.list <- topTable(eBayesFit, coef=1, number=Inf, sort.by='logFC')
gene.symbols <- getSYMBOL(row.names(gene.list), "hgu133plus2.db")
gene.results <- cbind(gene.list, gene.symbols)
write.table(gene.results, file="gene_details.txt", sep="\t")

# find the top 10 up-regulated and top 10 down-regulated genes
# write them to a file for convenience
gene.up <- head(gene.results[gene.results$logFC > 0, ], n=10)
gene.down <- head(gene.results[gene.results$logFC < 0, ], n=10)
write.table(gene.up, file="gene_up.txt", sep="\t")
write.table(gene.down, file="gene_down.txt", sep="\t")

# convert NA values to symbols for consistent handling
gene.symbols[is.na(gene.symbols)] <- "????"
probe.matches <- gene.symbols == "EGFR" | gene.symbols == "SOX2" | gene.symbols == "TP53"
probe.not.positions <- which(!probe.matches)

# construct a better indicative vector for labels
# label format: "(probe_id, gene_symbol)"
probe.labels <- paste("(", row.names(gene.results), ", ", gene.symbols, ")", sep="")
# set unwanted probe labels to NA, for elimination
probe.labels[probe.not.positions] <- NA
write.table(probe.labels[!is.na(probe.labels)], file="gene_probe_map.txt")
ggplot(gene.list, 
  aes(
     label=probe.labels,
     x=gene.list$logFC, 
     y=-log10(gene.list$P.Value))) + 
  geom_point(aes(shape=probe.matches, 
                 colour=probe.matches)) + 
  geom_text(size=3) + 
  xlab("Log2 Fold-Change") + 
  ylab("-Log10 P-Value") + 
  scale_shape_discrete(name="Shape", breaks=c(TRUE, FALSE), labels=c("EGFR, TP53, SOX2", "Other")) +
  scale_colour_discrete(name="Color", breaks=c(TRUE, FALSE), labels=c("EGFR, TP53, SOX2", "Other")) +
  ggtitle("Volcano Plot(Fold-Change of Genes vs P-Value)")

