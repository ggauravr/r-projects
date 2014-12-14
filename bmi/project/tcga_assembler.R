workingDirectory <- "~/Downloads/TCGA-Assembler"
dataDirectory <- "/media/ggauravr/Other/Downloads/RawData.TCGA-Assembler"
processedDirectory <- "/media/ggauravr/Other/Downloads/ProcessedData.TCGA-Assembler"

setwd(workingDirectory)

source("Module_A.r")
source("Module_B.r")

# DirectoryTraverseResult_Jul-08-2014
data <- DownloadRNASeqData(
  traverseResultFile = "./DirectoryTraverseResult_Jan-30-2014.rda", 
  saveFolderName = dataDirectory, 
  cancerType = "LUAD", 
  assayPlatform = "RNASeqV2", 
  dataType = c("rsem.genes.normalized_results"));

DownloadClinicalData(
  traverseResultFile = "./DirectoryTraverseResult_Jan-30-2014.rda", 
  saveFolderName=dataDirectory, 
  cancerType = "LUAD", 
  clinicalDataType=c("patient"));


GeneExpData = ProcessRNASeqData(
  inputFilePath = paste(dataDirectory, "LUAD__unc.edu__illuminahiseq_rnaseqv2__rsem.genes.normalized_results__Jan-30-2014.txt", sep="/"),
  outputFileName = "LUAD__illuminahiseq_rnaseqv2__GeneExp", 
  outputFileFolder = processedDirectory, 
  dataType = "GeneExp", 
  verType = "RNASeqV2")

getPBarcodeFromSBarcode <- function(sBarcode){
  pBarcode <- unlist(strsplit(sBarcode, "-"))
  pBarcode <- paste(pBarcode[1], pBarcode[2], pBarcode[3], sep="-")
}

sampleBarcodes <- colnames(Data)
matchPatientBarcodes <- unlist(lapply(sampleBarcodes, getPBarcodeFromSBarcode))
nRows <- dim(Data)[1]
nCols <- dim(Data)[2]
Data <- rbind(Data, matchPatientBarcodes)

cData <- read.table(cDataFile, header=T, sep="\t", stringsAsFactors=F)
cDataSub <- cData[, c("bcr_patient_barcode", "diagnosis", "pathologic_stage")]
