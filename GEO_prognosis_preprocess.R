library(GEOquery)
library(tidyverse)
library(affy)
library(org.Hs.eg.db)
library(DBI)
library(RMariaDB)
orgdb <- org.Hs.eg.db

### DB connection
dbname <- "expr"
cnf <- list.files("/poolio/db", pattern = paste0(dbname, ".cnf"), full.names = T)
db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)

########################### load data ###########################

# get feature info and load in
library(hgu133plus2.db)
annodb <- hgu133plus2.db
finfo <- select(annodb, keys = keys(annodb, "PROBEID"),
                keytype = "PROBEID", 
                columns = c("SYMBOL","ALIAS","ENTREZID","ENSEMBL"))
finfo <- cbind(platform = "HG-U133_Plus_2",
               finfo)

# load sample info
si <- readRDS("GSE61804_sampleinfo_formatted.rds")

# download CEL
studylist <- unique(si$study_accession)
studylist <- studylist[!is.na(studylist)]
names(studylist) <- studylist
errors <- lapply(studylist, function(gse) {
  print(paste0("Processing ", gse))
  if (!grepl("^GSE", gse)) {
    print("Study does not appear to be a GEO series - skipping.")
    return(NULL)
  }
  gse_si <- si[which(si$study_accession == gse),]
  gse_dir <- file.path("cel", gse)
  dir.create(gse_dir, showWarnings = F)
  lapply(gse_si$sample_accession, function(gsm) {
    filename <- file.path(gse_dir, paste0(gsm, ".CEL.gz"))
    if (!file.exists(filename)) {
      getGEOSuppFiles(gsm, makeDirectory = F, baseDir = gse_dir)
      defaultname <- list.celfiles(gse_dir, full.names = T,
                                   pattern = as.character(gsm))
      tryCatch(file.rename(defaultname, filename),
               error = function(e) {
                 print(paste0("ERROR: ", gsm))
                 return(gsm) 
               })
      return(NULL)
    }
  }) %>% unlist()
})
errors[sapply(errors, is.null)] <- NULL

# filter out samples with issues
si <- si[!si$sample_accession %in% unlist(errors),]

# normalise data in batches by study
errors <- lapply(studylist, function(study) {
  print(paste0("Processing ", study))
  outfile <- file.path("es" , paste0(study, "_log2MAS5_ES.rds"))
  if (file.exists(outfile)) {
    print("Normalised ES already exists")
    return(NULL)
  } else {
    celdir <- file.path("cel", study)
    celfiles <- list.celfiles(celdir)
    names(celfiles) <- sub("^(.+).CEL.gz$", "\\1", celfiles)
    pdat <- as.data.frame(si[match(names(celfiles), si$sample_accession),])
    row.names(pdat) <- celfiles
    e <- tryCatch(affy <- ReadAffy(filenames = celfiles, 
                                   sampleNames = names(celfiles),
                                   phenoData = pdat, 
                                   celfile.path = celdir), 
                  error = function(e) {
                    print(paste0("Couldn't process study: ", study))
                    print(e)
                    return(e) })
    if (inherits(e, "error")) return(study)
    affylist <- suppressWarnings(split(affy, ceiling(seq_along(affy) / 500)))
    normlist <- lapply(affylist, function(affysub) {
      mas5(affysub, sc = 500) })
    norm <- Reduce(Biobase::combine, normlist)
    sampleNames(norm) <- norm$sample_accession
    exprs(norm) <- log2(exprs(norm))
    saveRDS(norm, outfile)
    return(NULL)
  }
})

# add to db
lapply(studylist, function(study) {
  print(paste0("Processing ", study))
  esfile <- file.path("es" , paste0(study, "_log2MAS5_ES.rds"))
  es <- readRDS(esfile)
  expr_long <- cbind(feature_id = featureNames(es), exprs(es)) %>%
    as.data.frame() %>%
    gather("sample_id", "expr_val", -1)
  exprlist <- split(expr_long, ceiling(1:nrow(expr_long) / 1e7))
  lapply(exprlist, function(sub) {
    print(".")
    write_delim(sub, "tmp.txt", delim = "\t", 
                col_names = F, na = '\\N')
    (query <- paste0("LOAD DATA LOCAL INFILE '", getwd(),
                     "/tmp.txt' INTO TABLE expr_data;"))
    queryres <- dbSendQuery(db, query)
    dbClearResult(queryres)
  })
})

# get list of all successfully added sample ids
query <- paste0("SELECT sample_accession FROM expr_data ",
                "WHERE feature_id = '1007_s_at';")
queryres <- dbSendQuery(db, query)
sample_acc <- dbFetch(queryres)
sample_acc <- unique(sample_acc$sample_accession)
dbClearResult(queryres)

# add processed si to table
si <- si[si$sample_accession %in% sample_acc,]
write_delim(si, "tmp.txt", delim = "\t", 
            col_names = F, na = '\\N')
(query <- paste0("LOAD DATA LOCAL INFILE '", getwd(),
                 "/tmp.txt' IGNORE INTO TABLE sample_info;"))
queryres <- dbSendQuery(db, query)
dbClearResult(queryres)
dbDisconnect(db)
