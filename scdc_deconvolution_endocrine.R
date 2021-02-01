# Dieses Skript dient zur Dekonvolution der Benchmark Datensets (bulk RNA-seq: healthy, panNEN, GEP-NEN)
# auf Basis der Trainings Datensets (scRNA-seq: healthy) mit SCDC und Ensemble 
# Trainings Datensets: Baron, Segerstolpe, Lawlor (SCDC+artdeco Quellen), Haber (GEO?)
# Benchmark Datensets: Califano, Missiaglia, Fadista, RepSet, Riemer, Sadanandam, Scarpa
#!! NOTE !!# 
# This will be endocrine only model (alpha, beta, gamma, delta)

#### load Metadata file ####
library(SCDC)
metadata <- read.table(file = "~/Praktikum/Data/Clinial_data_Meta_information.tsv", sep = "\t", header = T)



#### load Benchmark data sets ####
# -> need true proportions?

# Califano: no grading, panNEN (i think, GEP-NEN were taken out, check with metadata again)
alvarez <- read.table(file = "~/Praktikum/Data/Alverez.S105.tsv", sep = '\t', header = TRUE) 
alvarez_meta_idx <- which(metadata$Study == "Alvarez") 
alvarez_meta <- metadata[alvarez_meta_idx,]
alvarez_meta <- alvarez_meta[which(alvarez_meta$Location=="pancreas"),]
alvarez_meta <- alvarez_meta[which(alvarez_meta$Subtype!= "Outlier"),] #105!

# Fadista: all G0 (healthy), pancreatic
fadista <- read.table(file = "~/Praktikum/Data/Fadista.S89.tsv", sep = '\t', header = TRUE)
rownames(fadista) <- fadista[,1]
fadista <- fadista[,-1]
fadista_meta_idx <- which(metadata$Study == "Fadista")
fadista_meta <- metadata[fadista_meta_idx,]

# Missiaglia: grading, K-67
missiaglia <- read.table(file = "~/Praktikum/Data/Missaglia.S75.tsv", sep = '\t', header = TRUE)
missiaglia_meta_idx <- which(metadata$Study == "GSE73338")
missiaglia_meta <- metadata[missiaglia_meta_idx,] 
missiaglia_sample_idx <- which(missiaglia_meta$Name %in% colnames(missiaglia))
missiaglia_meta <- missiaglia_meta[missiaglia_sample_idx,]

# Riemer: grading (!!NOT RIGHT!!), NEC + NET (no Controls, only Cancer) panNEN + GEP-NEN
riemer <- read.table(file = "~/Praktikum/Data/Riemer.S40.tsv", sep = '\t', header = TRUE)
rownames(riemer) <- riemer[,1]
riemer <- riemer[,-1]
riemer_meta_idx <- which(metadata$Study == "Riemer") 
riemer_meta <- metadata[riemer_meta_idx,]
riemer_meta <- riemer_meta[which(riemer_meta$Subtype=="Cancer"),] #40!

# Sadanandam: grading (!!NOT RIGHT!!), panNEN (but others too?)
sadanandam <- read.table(file = "~/Praktikum/Data/Sadanandam.S29.tsv", sep = '\t', header = TRUE)
sad_meta_idx <- which(metadata$Study == "Sadanandam") 
sad_meta <- metadata[sad_meta_idx,]

# Scarpa: grading, panNEN, NET+NEC
scarpa <- read.table(file = "~/Praktikum/Data/Scarpa.S29.tsv", sep = '\t', header = TRUE)
scarpa_meta_idx <- which(metadata$Study == "Scarpa") 
scarpa_meta <- metadata[scarpa_meta_idx,]

# RepSet: Scarpa+Riemer (panNEN+GEP-NEN), grading (!!NOT RIGHT!!)
repset <- read.table(file = "~/Praktikum/Data/RepSet.S57.tsv", sep = '\t', header = TRUE)
repset_meta <- rbind(scarpa_meta, riemer_meta)
rownames(repset_meta) <- repset_meta$Name
rownames(repset_meta)[1:55] <- paste("X",rownames(repset_meta)[1:55], sep = "")
repset_meta$Name <- as.character(repset_meta$Name)
repset_meta$Name[1:55] <- paste("X", repset_meta$Name[1:55], sep = "")
repset_meta_idx <- which(repset_meta$Name %in% colnames(repset))
repset_meta <- repset_meta[repset_meta_idx,]



#### Turn Benchmark Data into ExpressionSet objects ####
fdata_alvarez <- rownames(alvarez)
pdata_alvarez <- cbind(cellname = colnames(alvarez), subjects = as.character(alvarez_meta$Subtype))
eset_alvarez <- getESET(alvarez, fdata = fdata_alvarez, pdata = pdata_alvarez)

fdata_fadista <- rownames(fadista)
pdata_fadista <- cbind(cellname = colnames(fadista), subjects = as.character(fadista_meta$Grading))
eset_fadista <- getESET(fadista, fdata = fdata_fadista, pdata = pdata_fadista)

fdata_missiaglia <- rownames(missiaglia)
pdata_missiaglia <- cbind(cellname = colnames(missiaglia), subjects = as.character(missiaglia_meta$Grading))
eset_missiaglia <- getESET(missiaglia, fdata = fdata_missiaglia, pdata = pdata_missiaglia)

fdata_riemer <- rownames(riemer)
pdata_riemer <- cbind(cellname = colnames(riemer), subjects = as.character(riemer_meta$Grading))
eset_riemer <- getESET(riemer, fdata = fdata_riemer, pdata = pdata_riemer)

fdata_sadanandam <- rownames(sadanandam)
pdata_sadanandam <- cbind(cellname = colnames(sadanandam), subjects = as.character(sad_meta$Grading))
eset_sadanandam <- getESET(sadanandam, fdata = fdata_sadanandam, pdata = pdata_sadanandam)

fdata_scarpa <- rownames(scarpa)
pdata_scarpa <- cbind(cellname = colnames(scarpa), subjects = as.character(scarpa_meta$Grading))
eset_scarpa <- getESET(scarpa, fdata = fdata_scarpa, pdata = pdata_scarpa)

fdata_repset <- rownames(repset)
pdata_repset <- cbind(cellname = colnames(repset), subjects = as.character(repset_meta$Grading))
eset_repset <- getESET(repset, fdata = fdata_repset, pdata = pdata_repset)



#### load Training data sets ####
library(data.table)
library(Biobase)
library(tidyverse)
library(stringr)

## baron ##
#baron <- readRDS(file = "~/SCDC/SCDC/vignettes/data/baron.rds")
baron_meta_cnt <- readRDS(file = "~/Praktikum/Data/Baron/Baron.RDS")

baron_meta <- baron_meta_cnt[,1:4]
rownames(baron_meta) <- baron_meta$X
colnames(baron_meta)[4] <- "cluster"
colnames(baron_meta)[1] <- "sample"

baron_cnt <- baron_meta_cnt[,5:ncol(baron_meta_cnt)]
rownames(baron_cnt) <- baron_meta$X
baron_cnt <- t(baron_cnt)
baron_eset <- getESET(exprs = baron_cnt, fdata = rownames(baron_cnt), pdata = baron_meta)


## segerstolpe ##
#segerstolpe <- readRDS(file = "~/SCDC/SCDC/vignettes/data/segerstolpe.rds")
segerstolpe_meta_cnt <- readRDS(file = "~/Praktikum/Data/Segerstolpe/segerstolpe_raw.RDS")

segerstolpe_cnt <- segerstolpe_meta_cnt[-1,]
segerstolpe_cnt <- segerstolpe_cnt[,-2]
segerstolpe_cnt <- distinct(segerstolpe_cnt)
segerstolpe_cnt <- segerstolpe_cnt[!duplicated(segerstolpe_cnt[,1]),]
rownames(segerstolpe_cnt) <- segerstolpe_cnt[,1]
segerstolpe_cnt <- segerstolpe_cnt[,-1]

segerstolpe_meta <- t(rbind(colnames(segerstolpe_cnt), segerstolpe_meta_cnt[1,-c(1,2)]))
#colnames(segerstolpe_meta) <- c("name", "cluster")
colnames(segerstolpe_meta) <- c("name", "subtype")
# remove T2D's
t2d_idx <- which(sapply(1:nrow(segerstolpe_meta), function(x) str_detect(rownames(segerstolpe_meta)[x], "T2D")))
segerstolpe_meta <- segerstolpe_meta[-t2d_idx,]
segerstolpe_cnt <- segerstolpe_cnt[,-t2d_idx]
# add column about sample from segerstolpe@phenoData@data[,34:36] -> just match and see if cluster is same
seger_scdc_pheno <- cbind(rownames(segerstolpe@phenoData@data), segerstolpe@phenoData@data[,34:36])
pheno_idx <- match(unname(segerstolpe_meta[,1]), seger_scdc_pheno[,1])
seger_scdc_pheno <- seger_scdc_pheno[pheno_idx,]
segerstolpe_meta <- cbind(segerstolpe_meta, seger_scdc_pheno)

seger_tmp <- apply(segerstolpe_cnt, 2, unlist)
seger_tmp <- apply(seger_tmp, 2, as.numeric)
rownames(seger_tmp) <- rownames(segerstolpe_cnt)
segerstolpe_cnt <- seger_tmp

segerstolpe_eset <- getESET(exprs = segerstolpe_cnt, fdata = rownames(segerstolpe_cnt), pdata = segerstolpe_meta)


## lawlor ###
dt = fread("~/artdeco/artdeco/data/Lawlor.csv.gz") # expression values
lawlor_dt <- as.matrix(dt[,-1])
rownames(lawlor_dt) <- dt$V1

lawlor_meta <- metadata[which(metadata$Study=="Lawlor"),] # meta data
lawlor_idx <- which(lawlor_meta$Name %in% colnames(lawlor_dt))
lawlor_meta <- lawlor_meta[lawlor_idx,]
rownames(lawlor_meta) <- lawlor_meta$Name
lawlor_meta$Subtype <- tolower(lawlor_meta$Subtype)
lawlor_meta <- lawlor_meta[colnames(lawlor_dt),]
lawlor_meta_diff <- lawlor_meta
colnames(lawlor_meta_diff)[1] <- "sample" 
colnames(lawlor_meta_diff)[11] <- "cluster"

# single cell: 
lawlor <- getESET(exprs = lawlor_dt, fdata = dt$V1, pdata = lawlor_meta_diff)


# lawlor_gz <- gzfile('~/artdeco/artdeco/data/Lawlor.csv.gz','rt')  
# lawlor_csv <- read.csv(lawlor_gz,header=T, sep = ';') 
# lawlor_meta <- metadata[which(metadata$Study=="Lawlor"),]
# lawlor_idx <- which(lawlor_meta$Name %in% colnames(lawlor_csv))
# lawlor_meta <- lawlor_meta[lawlor_idx,]
# rownames(lawlor_meta) <- lawlor_meta$Name
# lawlor_meta$Subtype <- tolower(lawlor_meta$Subtype)
# #write.csv(lawlor_meta, file = "~/Praktikum/Data/Lawlor_meta.csv", quote = FALSE, row.names = T, col.names = T)
# #lawlor_meta_test <- read.csv(file = "~/Praktikum/Data/Lawlor_meta.csv", header = T, sep = ",", row.names = 1)
# lawlor_subtype_idx <- which(lawlor_meta$Subtype %in% c("alpha", "beta", "gamma", "delta"))
# lawlor_subtype <- lawlor_meta[lawlor_subtype_idx,]
# lawlor_subtype_csv <- lawlor_csv[,lawlor_subtype_idx]
# 
# fdata <- rownames(lawlor_subtype_csv)
# pData <- cbind(cellname = colnames(lawlor_subtype_csv), subjects = lawlor_subtype$Subtype)
# lawlor <- getESET(lawlor_subtype_csv, fdata = fdata, pdata = pData)
#
# phenoData <- new("AnnotatedDataFrame", data=lawlor_subtype)
# featureData <- new("AnnotatedDataFrame", data=data.frame(rownames(lawlor_subtype_csv)))
# rownames(featureData@data) <- rownames(lawlor_subtype_csv)
#
# lawlor <- ExpressionSet(assayData=as.matrix(lawlor_subtype_csv), phenoData = phenoData, featureData = featureData)



#### SCDC QC on scRNA-seq data ####
library(SCDC)
qc_baron <- SCDC_qc(baron_eset, ct.varname = "cluster", sample = "sample", 
                    scsetname = "Baron", ct.sub = c("alpha","beta","delta","gamma"), 
                    qcthreshold = 0.7)
qc_segerstolpe <- SCDC_qc(segerstolpe_eset, ct.varname = "cluster", sample = "sample", 
                          scsetname = "Segerstolpe", ct.sub = c("alpha","beta","delta","gamma"), 
                          qcthreshold = 0.7)
qc_lawlor <- SCDC_qc_ONE(lawlor, ct.varname = "cluster", sample = "sample", scsetname = "Lawlor",
                         ct.sub = c("alpha","beta","delta","gamma"))



#### Ensemble deconvolution ####
# 7 times
pancreas.sc <- list(baron = qc_baron$sc.eset.qc,
                    seger = qc_segerstolpe$sc.eset.qc,
                    lawlor = qc_lawlor$sc.eset.qc)

ensemble_alvarez <- SCDC_ENSEMBLE(bulk.eset = eset_alvarez, sc.eset.list = pancreas.sc, ct.varname = "cluster",
                             sample = "sample", ct.sub =  c("alpha","beta","delta","gamma"), search.length = 0.01)
ensemble_alvarez$w_table


ensemble_fadista <- SCDC_ENSEMBLE(bulk.eset = eset_fadista, sc.eset.list = pancreas.sc, ct.varname = "cluster",
                                  sample = "sample", ct.sub =  c("alpha","beta","delta","gamma"), search.length = 0.01)
ensemble_fadista$w_table


ensemble_missiaglia <- SCDC_ENSEMBLE(bulk.eset = eset_missiaglia, sc.eset.list = pancreas.sc, ct.varname = "cluster",
                                  sample = "sample", ct.sub =  c("alpha","beta","delta","gamma"), search.length = 0.01)
ensemble_missiaglia$w_table


ensemble_riemer <- SCDC_ENSEMBLE(bulk.eset = eset_riemer, sc.eset.list = pancreas.sc, ct.varname = "cluster",
                                  sample = "sample", ct.sub =  c("alpha","beta","delta","gamma"), search.length = 0.01)
ensemble_riemer$w_table


ensemble_sadanandam <- SCDC_ENSEMBLE(bulk.eset = eset_sadanandam, sc.eset.list = pancreas.sc, ct.varname = "cluster",
                                  sample = "sample", ct.sub =  c("alpha","beta","delta","gamma"), search.length = 0.01)
ensemble_sadanandam$w_table


ensemble_scarpa <- SCDC_ENSEMBLE(bulk.eset = eset_scarpa, sc.eset.list = pancreas.sc, ct.varname = "cluster",
                                  sample = "sample", ct.sub =  c("alpha","beta","delta","gamma"), search.length = 0.01)
ensemble_scarpa$w_table


ensemble_repset <- SCDC_ENSEMBLE(bulk.eset = eset_repset, sc.eset.list = pancreas.sc, ct.varname = "cluster",
                                  sample = "sample", ct.sub =  c("alpha","beta","delta","gamma"), search.length = 0.01)
ensemble_repset$w_table


ensemble_endocrine_dec_res <- list(ensemble_alvarez, ensemble_fadista, ensemble_missiaglia, ensemble_repset, ensemble_riemer,
                                  ensemble_sadanandam, ensemble_scarpa)
names(ensemble_endocrine_dec_res) <- c("Alvarez", "Fadista", "Missiaglia", "RepSet", "Riemer", "Sadanandam", "Scarpa")
#saveRDS(ensemble_endocrine_dec_res, file = "~/Praktikum/decon_res/ensemble_endocrine_dec_res.RDS")
saveRDS(ensemble_endocrine_dec_res, file = "~/Praktikum/decon_res2/ensemble_endocrine_dec_res2.RDS")


#### SCDC deconvolution ####
## Baron 7 times ##
baron_alvarez_deco <- SCDC_prop(bulk.eset = eset_alvarez, sc.eset = qc_baron$sc.eset.qc, 
                                    ct.varname = "cluster", sample = "sample", 
                                    ct.sub = c("alpha","beta","delta","gamma"))

baron_fadista_deco <- SCDC_prop(bulk.eset = eset_fadista, sc.eset = qc_baron$sc.eset.qc, 
                                ct.varname = "cluster", sample = "sample", 
                                ct.sub = c("alpha","beta","delta","gamma"))

baron_missiaglia_deco <- SCDC_prop(bulk.eset = eset_missiaglia, sc.eset = qc_baron$sc.eset.qc, 
                                ct.varname = "cluster", sample = "sample", 
                                ct.sub = c("alpha","beta","delta","gamma"))

baron_riemer_deco <- SCDC_prop(bulk.eset = eset_riemer, sc.eset = qc_baron$sc.eset.qc, 
                                   ct.varname = "cluster", sample = "sample", 
                                   ct.sub = c("alpha","beta","delta","gamma"))

baron_sadanandam_deco <- SCDC_prop(bulk.eset = eset_sadanandam, sc.eset = qc_baron$sc.eset.qc, 
                                   ct.varname = "cluster", sample = "sample", 
                                   ct.sub = c("alpha","beta","delta","gamma"))

baron_scarpa_deco <- SCDC_prop(bulk.eset = eset_scarpa, sc.eset = qc_baron$sc.eset.qc, 
                                   ct.varname = "cluster", sample = "sample", 
                                   ct.sub = c("alpha","beta","delta","gamma"))

baron_repset_deco <- SCDC_prop(bulk.eset = eset_repset, sc.eset = qc_baron$sc.eset.qc, 
                                   ct.varname = "cluster", sample = "sample", 
                                   ct.sub = c("alpha","beta","delta","gamma"))

baron_endocrine_dec_res <- list(baron_alvarez_deco, baron_fadista_deco, baron_missiaglia_deco, baron_repset_deco, 
                                baron_riemer_deco, baron_sadanandam_deco, baron_scarpa_deco)
names(baron_endocrine_dec_res) <- c("Alvarez", "Fadista", "Missiaglia", "RepSet", "Riemer", "Sadanandam", "Scarpa")
#saveRDS(baron_endocrine_dec_res, file = "~/Praktikum/decon_res/baron_endocrine_dec_res.RDS")
saveRDS(baron_endocrine_dec_res, file = "~/Praktikum/decon_res2/baron_endocrine_dec_res2.RDS")



## Segerstolpe 7 times ##
seger_alvarez_deco <- SCDC_prop(bulk.eset = eset_alvarez, sc.eset = qc_segerstolpe$sc.eset.qc, 
                                ct.varname = "cluster", sample = "sample", 
                                ct.sub = c("alpha","beta","delta","gamma"))

seger_fadista_deco <- SCDC_prop(bulk.eset = eset_fadista, sc.eset = qc_segerstolpe$sc.eset.qc, 
                                ct.varname = "cluster", sample = "sample", 
                                ct.sub = c("alpha","beta","delta","gamma"))

seger_missiaglia_deco <- SCDC_prop(bulk.eset = eset_missiaglia, sc.eset = qc_segerstolpe$sc.eset.qc, 
                                   ct.varname = "cluster", sample = "sample", 
                                   ct.sub = c("alpha","beta","delta","gamma"))

seger_riemer_deco <- SCDC_prop(bulk.eset = eset_riemer, sc.eset = qc_segerstolpe$sc.eset.qc, 
                               ct.varname = "cluster", sample = "sample", 
                               ct.sub = c("alpha","beta","delta","gamma"))

seger_sadanandam_deco <- SCDC_prop(bulk.eset = eset_sadanandam, sc.eset = qc_segerstolpe$sc.eset.qc, 
                                   ct.varname = "cluster", sample = "sample", 
                                   ct.sub = c("alpha","beta","delta","gamma"))

seger_scarpa_deco <- SCDC_prop(bulk.eset = eset_scarpa, sc.eset = qc_segerstolpe$sc.eset.qc, 
                               ct.varname = "cluster", sample = "sample", 
                               ct.sub = c("alpha","beta","delta","gamma"))

seger_repset_deco <- SCDC_prop(bulk.eset = eset_repset, sc.eset = qc_segerstolpe$sc.eset.qc, 
                               ct.varname = "cluster", sample = "sample", 
                               ct.sub = c("alpha","beta","delta","gamma"))

segerstolpe_endocrine_dec_res <- list(seger_alvarez_deco, seger_fadista_deco, seger_missiaglia_deco, seger_repset_deco, 
                                seger_riemer_deco, seger_sadanandam_deco, seger_scarpa_deco)
names(segerstolpe_endocrine_dec_res) <- c("Alvarez", "Fadista", "Missiaglia", "RepSet", "Riemer", "Sadanandam", "Scarpa")
#saveRDS(segerstolpe_endocrine_dec_res, file = "~/Praktikum/decon_res/segerstolpe_endocrine_dec_res.RDS")
saveRDS(segerstolpe_endocrine_dec_res, file = "~/Praktikum/decon_res2/segerstolpe_endocrine_dec_res2.RDS")


## Lawlor 7 times ##
lawlor_alvarez_deco <- SCDC_prop_ONE(bulk.eset = eset_alvarez, sc.eset = qc_lawlor$sc.eset.qc, 
                                ct.varname = "cluster", sample = "sample", 
                                ct.sub = c("alpha","beta","delta","gamma"))

lawlor_fadista_deco <- SCDC_prop_ONE(bulk.eset = eset_fadista, sc.eset = qc_lawlor$sc.eset.qc, 
                                ct.varname = "cluster", sample = "sample", 
                                ct.sub = c("alpha","beta","delta","gamma"))

lawlor_missiaglia_deco <- SCDC_prop_ONE(bulk.eset = eset_missiaglia, sc.eset = qc_lawlor$sc.eset.qc, 
                                   ct.varname = "cluster", sample = "sample", 
                                   ct.sub = c("alpha","beta","delta","gamma"))

lawlor_riemer_deco <- SCDC_prop_ONE(bulk.eset = eset_riemer, sc.eset = qc_lawlor$sc.eset.qc, 
                               ct.varname = "cluster", sample = "sample", 
                               ct.sub = c("alpha","beta","delta","gamma"))

lawlor_sadanandam_deco <- SCDC_prop_ONE(bulk.eset = eset_sadanandam, sc.eset = qc_lawlor$sc.eset.qc, 
                                   ct.varname = "cluster", sample = "sample", 
                                   ct.sub = c("alpha","beta","delta","gamma"))

lawlor_scarpa_deco <- SCDC_prop_ONE(bulk.eset = eset_scarpa, sc.eset = qc_lawlor$sc.eset.qc, 
                               ct.varname = "cluster", sample = "sample", 
                               ct.sub = c("alpha","beta","delta","gamma"))

lawlor_repset_deco <- SCDC_prop_ONE(bulk.eset = eset_repset, sc.eset = qc_lawlor$sc.eset.qc, 
                               ct.varname = "cluster", sample = "sample", 
                               ct.sub = c("alpha","beta","delta","gamma"))


lawlor_endocrine_dec_res <- list(lawlor_alvarez_deco, lawlor_fadista_deco, lawlor_missiaglia_deco, lawlor_repset_deco, 
                                lawlor_riemer_deco, lawlor_sadanandam_deco, lawlor_scarpa_deco)
names(lawlor_endocrine_dec_res) <- c("Alvarez", "Fadista", "Missiaglia", "RepSet", "Riemer", "Sadanandam", "Scarpa")
#saveRDS(lawlor_endocrine_dec_res, file = "~/Praktikum/decon_res/lawlor_endocrine_dec_res.RDS")
# no change here
saveRDS(lawlor_endocrine_dec_res, file = "~/Praktikum/decon_res2/lawlor_endocrine_dec_res2.RDS")
