library(tidyverse)
library(survival)
library(survminer)

# list datasets
datlist <- list.files("../ES", ".rds", full.names = T)
names(datlist) <- sub("^.+/(.+).rds", "\\1", datlist)

# for each
lapply(names(datlist), function(datname) {
  
  if (file.exists(paste0("../res/", datname, "_res.rds"))) return(NULL)
  print(datname)
  
  # load ES
  es <- readRDS(datlist[[datname]]) 
  es <- es[, which(!is.na(es$follow_up_years_death) &
                     !is.na(es$censor_death))]
  
  # list all genes
  fdat <- fData(es) %>%
    separate_rows(SYMBOL)
  genelist <- unique(fdat$SYMBOL)
  genelist <- genelist[which(genelist != "NA" & !is.na(genelist))]
  names(genelist) <- genelist
  
  # list disease subtypes
  subtypes <- table(es$ds_subtype)
  print(subtypes)
  subtypes <- subtypes[subtypes > 10]
  
  # subtypes
  stdf <- lapply(names(subtypes), function(subtype) {
    print(subtype)
    # subset data
    essub <- es[, which(es$ds_subtype==subtype)]
    time <- essub$follow_up_years_death
    censor <- essub$censor_death
    # apply to each gene
    lapply(genelist, function(gene) {
      probes <- unique(fdat[which(fdat$SYMBOL == gene), ]$PROBEID)
      idx <- which(fData(essub)$PROBEID %in% probes)
      if (length(idx) > 1) { expr <- colMeans(exprs(essub)[idx, ]) 
      } else { expr <- exprs(essub)[idx, ] }
      fit <- coxph(Surv(time, censor) ~ expr)
      fitsum <- summary(fit)
      cbind(fitsum$coefficients, fitsum$conf.int[,3:4, drop=F]) %>%
        as.data.frame() %>%
        mutate(gene = gene, 
               n = fitsum$n,
               .before=1)
    }) %>% bind_rows() %>%
      mutate(ds_group = essub$ds_group[1],
             ds_type = essub$ds_type[1],
             ds_subtype = essub$ds_subtype[1],
             .before=1)
  }) %>% bind_rows()
  
  # overall (disease type)
  time <- es$follow_up_years_death
  censor <- es$censor_death
  oadf <- lapply(genelist, function(gene) {
    probes <- unique(fdat[which(fdat$SYMBOL == gene), ]$PROBEID)
    idx <- which(fData(es)$PROBEID %in% probes)
    if (length(idx) > 1) { expr <- colMeans(exprs(es)[idx, ])
    } else { expr <- exprs(es)[idx, ] }
    fit <- coxph(Surv(time, censor) ~ expr)
    fitsum <- summary(fit)
    cbind(fitsum$coefficients, fitsum$conf.int[,3:4, drop=F]) %>%
      as.data.frame() %>%
      mutate(gene = gene,
             n = fitsum$n,
             .before=1)
  }) %>% bind_rows() %>%
    mutate(ds_group = es$ds_group[1],
           ds_type = es$ds_type[1],
           ds_subtype = "overall",
           .before=1)
  
  # combine, add study info and return
  studydf <- stdf %>% mutate(dataset = datname, .before=1)
  studydf <- rbind(oadf, stdf) %>% mutate(dataset = datname, .before=1)
  saveRDS(studydf, paste0("../res/", datname, "_res.rds"))
}) 

# load res
resfiles <- list.files("../res", ".rds", full.names = T)
resfiles <- resfiles[!grepl("prognosis", resfiles)]
genedf <- lapply(resfiles, function(x) readRDS(x)) %>% 
  bind_rows() %>%
  distinct() %>%
  dplyr::filter(!gene %in% c("NA","3") & !is.na(gene))
colnames(genedf)[8:13] <- c("exp_coef","se_coef","z","pval","lower_95ci","upper_95ci")

# meta-effects
library(meta)
metadf <- lapply(unique(genedf$gene), function(genename) {
  genedat <- genedf %>% dplyr::filter(gene==genename)
  stmeta <- lapply(unique(genedat$ds_subtype), function(st) {
    stdat <- genedat[which(genedat$ds_subtype == st), ]
    dstypes <- unique(stdat$ds_type)
    lapply(dstypes, function(ds) {
      dsdat <- stdat[which(stdat$ds_type == ds), ]
      tryCatch({
        metares <- metagen(TE = dsdat$coef, seTE = dsdat$se_coef, studlab = dsdat$study,
                         fixed = FALSE, random = TRUE)
        data.frame(dataset = "meta", 
                   ds_group = dsdat$ds_group[1],
                   ds_type = dsdat$ds_type[1],
                   ds_subtype = st,
                   gene = genename,
                   n = sum(dsdat$n),
                   coef = metares$TE.random,
                   exp_coef = exp(metares$TE.random),
                   se_coef = metares$seTE.random,
                   z = metares$statistic.random,
                   pval = metares$pval.random,
                   lower_95ci = exp(metares$lower.random),
                   upper_95ci = exp(metares$upper.random))
      }, error = function(e) return(NULL))
    }) %>% bind_rows()
  }) %>% bind_rows()
  dsmeta <- lapply(unique(genedat$ds_type), function(ds) {
    dsdat <- genedat[which(genedat$ds_type == ds), ]
    tryCatch({
      metares <- metagen(TE = dsdat$coef, seTE = dsdat$se_coef, fixed = FALSE, random = TRUE)
      data.frame(dataset = "meta", 
                 ds_group = dsdat$ds_group[1],
                 ds_type = ds,
                 ds_subtype = NA,
                 gene = genename,
                 n = sum(dsdat$n),
                 coef = metares$TE.random,
                 exp_coef = exp(metares$TE.random),
                 se_coef = metares$seTE.random,
                 z = metares$statistic.random,
                 pval = metares$pval.random,
                 lower_95ci = exp(metares$lower.random),
                 upper_95ci = exp(metares$upper.random))
    }, error = function(e) return(NULL))
  }) %>% bind_rows()
  rbind(dsmeta, stmeta)
}) %>% bind_rows()
resdf <- rbind(genedf, metadf)
saveRDS(resdf, "../res/prognosis_gene_scores.rds")

# add solid tumor total
(groups <- unique(genedf$ds_group))
solid_groups <- groups[!groups %in% c("Lymphoid","Myeloid" )]
solidsub <- metadf %>%
  dplyr::filter(ds_group %in% solid_groups & is.na(ds_subtype))
solidmeta <- lapply(unique(solidsub$gene), function(genename) {
  genedat <- solidsub %>% dplyr::filter(gene==genename)
  tryCatch({
    metares <- metagen(TE = genedat$coef, seTE = genedat$se_coef, fixed = FALSE, random = TRUE)
    data.frame(dataset = "meta", 
               ds_group = "Solid tumor",
               ds_type = "Solid tumor",
               ds_subtype = NA,
               gene = genename,
               n = sum(genedat$n),
               coef = metares$TE.random,
               exp_coef = exp(metares$TE.random),
               se_coef = metares$seTE.random,
               z = metares$statistic.random,
               pval = metares$pval.random,
               lower_95ci = exp(metares$lower.random),
               upper_95ci = exp(metares$upper.random))
  }, error = function(e) return(NULL))
}) %>% bind_rows()
resdf <- rbind(resdf, solidmeta)
resdf <- resdf %>% distinct()
saveRDS(resdf, "../res/prognosis_gene_scores.rds")

