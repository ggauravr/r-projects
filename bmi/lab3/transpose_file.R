inputData.original <- read.delim("<tab-spaced or comma-separated input filename>", header=TRUE, sep="\t")

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

# save cleaned input data
# the below save command works for any R variable or Object that you wanna save
# to avoid recomputations and just load it back with the load command
save(inputData.trans, file="SavedInput.RData")

# just uncomment the below line to load the variable back into your R Environement
# load("SavedInput.RData")
