###################################
# CSE 5599/BMI 7830 Lab I
#
# Goal:
# Differential Gene Expression study of NSCLC(non-small cell lung cancer) sub-types: 
#   Adenocarcinoma(AC) and Squamous Cell Carcinoma(SCC)
#
###################################

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

############## 
# 1. Data Preparation
##############

funcGetDataFiles <- function(){
  # download the dataset, extract the .tar file
  # unzip the resulting gzipped files, to get the .cel files
  getGEOSuppFiles(constDataset)
  untar(paste(constDataset, "/", constDataset, "_RAW.tar", sep=""), exdir=constDataDirectory)
  sapply(paste(constDataDirectory, list.files(constDataDirectory, pattern = "gz$") , sep="/"), gunzip)
}
# funcGetDataFiles()

##############################
# 2. Annotation and AffyBatch object creation
##############################

# generate AffyBatch object from CEL and corresponding Annotation files
# normalize the expression data using GC-RMA method and store the result
# ExpressionSet object into a file for convience in future
celAffyBatch <- read.affy(covdesc="GDS3627_SampleAnnotation.txt", path=constDataDirectory)

#############
# 3. Normalization
#############

##########################################
# uncomment the below two lines, when run for the first time, 
# i.e. if the expression data is not stored in a file
##########################################

# gcrmaEset <- gcrma(celAffyBatch)
# write.exprs(gcrmaEset, file=constEsetFile)

# load expressionSet object from the file saved above
# set the pheno and annotation data from the AffyBatch object
gcrmaEset <- readExpressionSet(constEsetFile)
pData(gcrmaEset) <- pData(celAffyBatch)
annotation(gcrmaEset) <- annotation(celAffyBatch)

#################
# 4. Visualization and QC
#################

funcPreVisualize <- function(){
  # set colour palette
  colors <- brewer.pal(8, "Set1")

  #####################################################
  # boxplot parameters:
  #   las = 3 indicates vertical naming in the x-axis
  #   cex.axis adjusts the font of axis names, in terms of fraction of the default size
  #####################################################

  # boxplot of unnormalized intensity values
  boxplot(celAffyBatch, col=colors, las=3, cex.axis=0.6)
  # plot a boxplot of normalized intensity values
  boxplot(gcrmaEset, col=colors, las=3, cex.axis=0.6)

  # histogram of unnormalized values
  hist(celAffyBatch, col=colors)
  # histogram of normalized values
  hist(gcrmaEset, col=colors)  
}
# funcPreVisualize()

# probe-level metric calculations on the CEL files
celAffyBatch.qc <- fitPLM(celAffyBatch)
# Create an image of the first .CEL:
image(celAffyBatch.qc, which=1, add.legend=TRUE)

# GSM258551.CEL 258572.CEL 258577 258587 258605 appear to be shifted from 0
RLE(celAffyBatch.qc, main="RLE", las=3, cex.axis=0.6)
# We can also use NUSE (Normalised Unscaled Standard Errors).
# The median standard error should be 1 for most genes.
# GSM258553 258570 258580 258590 258603 258551 258572 258584.CEL appear to have medaian standard error above 1
NUSE(celAffyBatch.qc, main="NUSE", las=3, cex.axis=0.6)

#############################
# 5. Comparative Analysis and Model Fitting
#############################

#####################################################
# celAffyBatch.filtered$filter.log shows the various features excluded and the reason of exclusion
# here, the only exclusions are from features where IQR < 0.5(default)
# and the Affymetrix Quality Control probesets, that is the probesets whose names
# begin with AFFY
# SEE: ?nsFilter to understand the various arguments and the output
#####################################################

celAffyBatch.filtered <- nsFilter(gcrmaEset, require.entrez=FALSE, remove.dupEntrez=FALSE)
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

#######################
# 6. ProbeSet to Gene Mapping
#######################

# get the probeset_ids and the relevant statistics
# map the probeset_ids to their corresponding gene symbols, using the annotation db
# combine/map probeset_ids and their statistics to gene symbols
# write the results into a file
gene.summary <- topTable(eBayesFit, coef=1, number=Inf, sort.by='logFC')
gene.symbols <- getSYMBOL(row.names(gene.summary), "hgu133plus2.db")
gene.details <- cbind(gene.summary, gene.symbols)
write.table(gene.details, file="gene_details.txt", sep="\t")

#########################
# 7. Regulated Gene finding
#########################

# find the top 10 up-regulated and top 10 down-regulated genes
# write them to a file for convenience
gene.up <- head(gene.details[gene.details$logFC > 0, ], n=11)
gene.down <- head(gene.details[gene.details$logFC < 0, ], n=10)
# write.table(gene.up, file="gene_up.txt", sep="\t")
# write.table(gene.down, file="gene_down.txt", sep="\t")

# convert NA values to invalid symbols for consistent handling
gene.symbols[is.na(gene.symbols)] <- "????"

###################
# 8. Volcano Plots
###################

# find the boolean array, matching the given criteria(match EGFR | SOX2 | TP53)
probe.matches <- gene.symbols == "EGFR" | gene.symbols == "SOX2" | gene.symbols == "TP53"
probe.not.positions <- which(!probe.matches)

# construct a better indicative vector for labels
# label format: "(probe_id, gene_symbol)"
probe.labels <- paste("(", row.names(gene.details), ", ", gene.symbols, ")", sep="")
# set unwanted probe labels to NA, for elimination
probe.labels[probe.not.positions] <- NA
# write.table(probe.labels[!is.na(probe.labels)], file="gene_probe_map.txt")

# volcano plot for fold-change vs p-values
# genes(EGFR, TP53, SOX2) are indicated in a distinct color
ggplot(gene.summary, 
  aes(
     label=probe.labels,
     x=gene.summary$logFC, 
     y=-log10(gene.summary$P.Value))) + 
  
  # probe.matches has TRUE for the the three genes(EGFR, TP53, SOX2) and FALSE for all others
  # these aesthetics help color/shape the three genes differently than the others
  geom_point(aes(shape=probe.matches, 
                 colour=probe.matches)) + 
  geom_text(size=3) + 
  xlab("Log2 Fold-Change") + 
  ylab("-Log10 P-Value") + 
  scale_shape_discrete(name="Shape", breaks=c(TRUE, FALSE), labels=c("EGFR, TP53, SOX2", "Other")) +
  scale_colour_discrete(name="Color", breaks=c(TRUE, FALSE), labels=c("EGFR, TP53, SOX2", "Other")) +
  ggtitle("Volcano Plot(Fold-Change vs P-Value)")

# find the boolean array for, genes with abs(logFC) > 3 and p < 0.05
probe.matches <- abs(gene.details$logFC) > 3 & gene.details$P.Value < 0.05
probe.not.positions <- which(!probe.matches)

# volcano plot for fold-change vs p-values
# probes with |logFC|> 3 and p-value < 0.05 are indicated in a distinct color
ggplot(gene.summary, 
  aes(
     x=gene.summary$logFC, 
     y=-log10(gene.summary$P.Value))) + 
  
  # probe.matches has TRUE for the genes that match the above criteria and FALSE for all others
  geom_point(aes(shape=probe.matches, 
                 colour=probe.matches)) + 
  
  xlab("Log2 Fold-Change") + 
  ylab("-Log10 P-Value") + 
  scale_shape_discrete(name="Shape", breaks=c(TRUE, FALSE), labels=c("Statistically Significant Genes", "Other")) +
  scale_colour_discrete(name="Color", breaks=c(TRUE, FALSE), labels=c("Statistically Significant Genes", "Other")) +
  ggtitle("Volcano Plot(Fold-Change vs P-Value)") +
  
  # indicate the fold-change cut-off
  geom_vline(xintercept=3, linetype=3) +
  geom_vline(xintercept=-3, linetype=3)+
  # indicate the p-value cut-off
  geom_hline(yintercept=-log10(0.05), linetype=3)

##############################
# 9. Enrichment Analysis
##############################

# select genes with |logFC| > 2, for enrichment analysis
gene.train <- gene.details[abs(gene.details$logFC) > 2, ]
# filter out NA values
gene.train <- gene.train[!is.na(gene.train$gene.symbols), ]

# get the probe_ids of the selected genes
# convert factors to characters for convenient output to a file
gene.train.probe_ids <- row.names(gene.train)
gene.train.symbols <- as.character(gene.train$gene.symbols)

# change the source of output
# format the output with one gene symbol per line(to paste into ToppFun/ToppGene)
# reset the output to terminal
sink(file="gene_train_symbols.txt")
cat(gene.train.symbols, sep="\n")
sink()

sink(file="gene_train_probeids.txt")
cat(gene.train.probe_ids, sep="\n")
sink()
