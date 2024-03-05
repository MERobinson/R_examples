library(tidyverse)
library(GEOquery)
library(limma)
library(fgsea)
library(affy)

###########################################################
##### Downloading and examining pre-processed data ########
###########################################################

#' microarray data in GEO can be downloaded with the getGEO function from GEOquery

# download data from GEO
geo <- getGEO("GSE23743") # this returns a list with all different platforms used in the study
geo <- geo[[1]] # select the platform you want to look at (theres only one here anyway)

#' this downloads the data as an ExpressionSet object - which contains 3 important parts
#' 1 - the actual assay data - i.e. the expression levels for each probe - can be accessed with exprs() function
#' 2 - the phenotype data - i.e. sample/column info - can be accessed with pData() function
#' 3 - the feature data - i.e. information about which genes the probes map to - accessed with fData() function 
head(exprs(geo))
head(pData(geo), 2)
head(fData(geo), 2)

#' NOTE: importantly, there is no standardization about exactly what the authors have to provide,
#' so you need to check each dataset to see what pData/fData has been provided, and whether the 
#' expression data has been a) logged, and b) normalized already.

# check if data is normalized/logged
boxplot(exprs(geo)) # you can see the majority of values are very low but with a long upper tail - indicating this isn't logged
exprs(geo) <- log2(exprs(geo)) # log the expression values
boxplot(exprs(geo)) # this data seems to have been pre-normalized (distributions are similar)

#' this example dataset was already RMA normalized and we have now logged the expression values
#' so we could go ahead with analysis and skip the next section, however the below section 
#' provides info about how to download the raw data and apply RMA normalization if needed

###########################################
##### Re-normalizing from raw data ########
###########################################

# can download the raw CEL files
getGEOSuppFiles("GSE23743",  makeDirectory = T) # download raw CEL files
untar("GSE23743/GSE23743_RAW.tar") # decompress
file_names <- list.celfiles("GSE23743/GSE23743_RAW") # list the CEL files
sample_names <- sub("^.*(GSM[0-9]+).+$", "\\1", file_names) # get the sample names (GEO sample IDs)

# sort out phenotype and feature information for creating the ExpressionSet object
head(pData(geo)) # look at what phenotype info was provided in the GEO record
pdat <- AnnotatedDataFrame(pData(geo)[,c("title","geo_accession","cell line:ch1")]) # take whatever columns are useful from the pdata
colnames(pdat) <- c("sample_id", "geo_accession", "cell_line") # rename columns
pdat$treatment <- ifelse(grepl("Imatinib-treated", pdat$sample_id),
                         "Imatinib", "control") # add treatment info - was missing in the GEO record
head(fData(geo)) # look at what info is provided in the feature data from GEO record
fdat <- fData(geo)[,c("ID","Gene Symbol","ENTREZ_GENE_ID", "Gene Title")] # take useful parts of the feature data
colnames(fdat) <- c("probe_id", "gene_symbol", "entrez_id", "gene_name") # rename

# load and normalize data into an ExpressionSet object
rma <- justRMA(celfile.path = ".", sampleNames = sample_names,
               phenoData = pdat, normalize = T)
fData(rma) <- fdat

###########################################
####### Differential expression  ##########
###########################################

#' note the columns names from pData need changing here if you did the above normalization step

# set sample info
treatment <- factor(sub(".+(Imatinib-treated|non-treated)","\\1",pData(geo)$title),
                    levels = c("non-treated","Imatinib-treated"),
                    labels = c("untreated","Imatinib"))
cell_line <- factor(gsub("\\-","",pData(geo)[,"cell line:ch1"]))
design <- model.matrix(~cell_line+treatment)
fit <- lmFit(exprs(geo), design=design)
bayes <- eBayes(fit)
deres <- topTable(bayes, sort = "none", n = Inf, coef="treatmentImatinib") %>%
  mutate(probe_id = fData(geo)[,"ID"],
         gene_symbol = fData(geo)[,"Gene Symbol"]) 
deres <- cbind(deres, exprs(geo)) %>%
  arrange(P.Value)
colnames(deres)[9:16] <- pData(geo)$title
write.csv(deres, "GSE23743_Imatinib_vs_untreated_Ph_DE_res.csv")
