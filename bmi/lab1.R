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
write.table(gene.up, file="genes_up_regulated.txt", sep="\t")
write.table(gene.down, file="genes_down_regulated.txt", sep="\t")

par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)
plot(
    # logFC on x- and -log p-value on y-axis
    gene.list$logFC, 
    -log10(gene.list$P.Value),
     
    # set x and y-axis limits
    xlim=c(-10, 10), 
    ylim=c(0, 15),

    # set x- and y-axis labels
    xlab="log2 fold change", ylab="-log10 p-value")


sum(abs(gene.list$logFC) > 2 & gene.list$P.Value < 0.001)

##no_of_genes=27,306
no_of_genes = dim(probeset.list)[1]
##Gives 693 genes
sum(abs(gene.list$logFC) > 2 & gene.list$P.Value < 0.05/no_of_genes)

install.packages("ggplot2")
require(ggplot2)



##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
gene.list$threshold = as.factor(abs(gene.list$logFC) > 2 & gene.list$P.Value < 0.05/no_of_genes)

##Construct the plot object
g = ggplot(data=gene.list, aes(x=lproogFC, y=-log10(P.Value), colour=threshold)) + 
  geom_point(alpha=0.4, size=1.75) + theme(legend.position = "none") + xlim(c(-10, 10)) + ylim(c(0, 15)) + xlab("log2 fold change") + ylab("-log10 p-value")
g

ind <- which(gene.list$gene.symbols == 'ESR1')
##Graph not shown
##g + geom_text(aes(x=gene.list[ind[1],]$logFC, y=-log10(gene.list[ind[1],]$P.Value), label=gene.list[ind[1],]$gene.symbols, size=1.2), colour="black")
g + geom_text(aes(x=gene.list[ind,]$logFC, y=-log10(gene.list[ind,]$P.Value), label=gene.list[ind,]$gene.symbols, size=1.2), colour="black")

gene.symbols[is.na(gene.symbols)] <- "XXXXX"
probe.matches <- gene.symbols == "EGFR" | gene.symbols == "SOX2" | gene.symbols == "TP53"
probe.not.positions <- which(!probe.matches)
probe.labels <- gene.symbols
probe.labels[probe.not.positions] <- NA
# probe.labels <-  gene.symbols[probe.positions]

ggplot(gene.list, 
        aes(
         label=probe.labels,
         x=gene.list$logFC, 
         y=-log10(gene.list$P.Value) )) + 
        geom_point(aes(shape=probe.matches, color=probe.matches)) + 
        geom_text(size=3) +
        xlab("log2 fold-change") + 
        ylab("-log10 p-value")


