# Dieses Skript dient zur Berechnung der Pearson Correlation zwischen
# 1) predicted cell type proportions with Ki-67 count levels
# 2) predicted cell type proportions with histopathology-based grading
# Input:
# Meta Data with grading info : ~/Praktikum/Data/Clinial_data_Meta_information.tsv
# 7 Benchmark Data sets with Ki-67 values : ~/Praktikum/Data/*[0-9].tsv
# 8 (12) times 7 deco results (endo, endo+ex)x(scdc-baron, scdc-segerstolpe, scdc-lawlor, ensemble) in ~/Praktikum/decon_res/*



#### load Metadata file ####
metadata <- read.table(file = "~/Praktikum/Data/Clinial_data_Meta_information.tsv", sep = "\t", header = T)



#### load Benchmark data sets ####
# -> need true proportions?

# Califano
alvarez <- read.table(file = "~/Praktikum/Data/Alverez.S105.tsv", sep = '\t', header = TRUE) 
alvarez_meta_idx <- which(metadata$Study == "Alvarez") 
alvarez_meta <- metadata[alvarez_meta_idx,]
alvarez_meta <- alvarez_meta[which(alvarez_meta$Location=="pancreas"),]
alvarez_meta <- alvarez_meta[which(alvarez_meta$Subtype!= "Outlier"),] #105!

# Fadista
fadista <- read.table(file = "~/Praktikum/Data/Fadista.S89.tsv", sep = '\t', header = TRUE)
rownames(fadista) <- fadista[,1]
fadista <- fadista[,-1]
fadista_meta_idx <- which(metadata$Study == "Fadista")
fadista_meta <- metadata[fadista_meta_idx,]

# Missiaglia
missiaglia <- read.table(file = "~/Praktikum/Data/Missaglia.S75.tsv", sep = '\t', header = TRUE)
missiaglia_meta_idx <- which(metadata$Study == "GSE73338")
missiaglia_meta <- metadata[missiaglia_meta_idx,] 
missiaglia_sample_idx <- which(missiaglia_meta$Name %in% colnames(missiaglia))
missiaglia_meta <- missiaglia_meta[missiaglia_sample_idx,]

# Riemer
riemer <- read.table(file = "~/Praktikum/Data/Riemer.S40.tsv", sep = '\t', header = TRUE)
rownames(riemer) <- riemer[,1]
riemer <- riemer[,-1]
riemer_meta_idx <- which(metadata$Study == "Riemer") 
riemer_meta <- metadata[riemer_meta_idx,]
riemer_meta <- riemer_meta[which(riemer_meta$Subtype=="Cancer"),] #40!

# Sadanandam
sadanandam <- read.table(file = "~/Praktikum/Data/Sadanandam.S29.tsv", sep = '\t', header = TRUE)
sad_meta_idx <- which(metadata$Study == "Sadanandam") 
sad_meta <- metadata[sad_meta_idx,]

# Scarpa
scarpa <- read.table(file = "~/Praktikum/Data/Scarpa.S29.tsv", sep = '\t', header = TRUE)
scarpa_meta_idx <- which(metadata$Study == "Scarpa") 
scarpa_meta <- metadata[scarpa_meta_idx,]

# RepSet
repset <- read.table(file = "~/Praktikum/Data/RepSet.S57.tsv", sep = '\t', header = TRUE)
repset_meta <- rbind(scarpa_meta, riemer_meta)
rownames(repset_meta) <- repset_meta$Name
rownames(repset_meta)[1:55] <- paste("X",rownames(repset_meta)[1:55], sep = "")
repset_meta$Name <- as.character(repset_meta$Name)
repset_meta$Name[1:55] <- paste("X", repset_meta$Name[1:55], sep = "")
repset_meta_idx <- which(repset_meta$Name %in% colnames(repset))
repset_meta <- repset_meta[repset_meta_idx,]



#### MKi-67 values ####
ki67_alvarez <- alvarez[which(rownames(alvarez) == "MKI67"),]
ki67_fadista <- fadista[which(rownames(fadista) == "MKI67"),]
ki67_missiaglia <- missiaglia[which(rownames(missiaglia) == "MKI67"),]
ki67_repset <- repset[which(rownames(repset) == "MKI67"),]
ki67_riemer <- riemer[which(rownames(riemer) == "MKI67"),]
ki67_sadanandam <- sadanandam[which(rownames(sadanandam) == "MKI67"),]
ki67_scarpa <- scarpa[which(rownames(scarpa) == "MKI67"),]

ki67_values <- list(ki67_alvarez, ki67_fadista, ki67_missiaglia, ki67_repset, ki67_riemer, ki67_sadanandam, ki67_scarpa)
names(ki67_values) <- c("Alvarez", "Fadista", "Missiaglia", "RepSet", "Riemer", "Sadanandam", "Scarpa")

for (i in 1:length(ki67_values)) {
  ki67_values[[i]] <- log(as.double(ki67_values[[i]])+1)
}

log(as.double(transcriptome_data[ki_index[1],])+1)



#### Grading info ####
#grading_fadista <- fadista_meta$Grading
#names(grading_fadista) <- fadista_meta$Name

grading_missiaglia <- missiaglia_meta$Grading
names(grading_missiaglia) <- missiaglia_meta$Name

grading_repset <- repset_meta$Grading
names(grading_repset) <- repset_meta$Name

grading_riemer <- riemer_meta$Grading
names(grading_riemer) <- riemer_meta$Name

grading_sadanandam <- sad_meta$Grading
names(grading_sadanandam) <- sad_meta$Name

grading_scarpa <- scarpa_meta$Grading
names(grading_scarpa) <- scarpa_meta$Name

grading_info <- list(grading_missiaglia, grading_repset, grading_riemer, grading_sadanandam, grading_scarpa)
names(grading_info) <- c("Missiaglia", "RepSet", "Riemer", "Sadanandam", "Scarpa")

grading_info_numeric <- sapply(1:length(grading_info), function(x) as.vector(grading_info[[x]]))
grading_info_numeric <- sapply(1:length(grading_info_numeric), function(x) 
  as.integer(str_replace_all(grading_info_numeric[[x]],pattern ="G","")))



#### Endocrine only model ####

#### read in deconvolution results ####
# deconvolution results (cell type proportions per sample) are in
# baron_endocrine_dec_res[[1]]$prop.est.mvw and
# ensemble_endocrine_dec_res[[1]]$prop.only$baron

#baron_endocrine_dec_res <- readRDS("~/Praktikum/decon_res/baron_endocrine_dec_res.RDS")
baron_endocrine_dec_res <- readRDS("~/Praktikum/decon_res2/baron_endocrine_dec_res2.RDS")
baron_endocrine_cell_prop <- list()
for (i in 1:length(baron_endocrine_dec_res)) {
  baron_endocrine_cell_prop[[i]] <- baron_endocrine_dec_res[[i]]$prop.est.mvw
}
names(baron_endocrine_cell_prop) <- names(baron_endocrine_dec_res)


#ensemble_endocrine_dec_res <- readRDS("~/Praktikum/decon_res/ensemble_endocrine_dec_res.RDS")
ensemble_endocrine_dec_res <- readRDS("~/Praktikum/decon_res2/ensemble_endocrine_dec_res2.RDS")
ensemble_endocrine_ref_by_weight <- list()
for (i in 1:length(ensemble_endocrine_dec_res)) {
  weights <- ensemble_endocrine_dec_res[[i]]$w_table[1:5,1:3]
  mean_weights <- sapply(1:3, function(x) mean(weights[,x]))
  ensemble_endocrine_ref_by_weight[[i]] <- which(mean_weights == max(mean_weights))
}
ensemble_endocrine_ref_by_weight <- as.integer(ensemble_endocrine_ref_by_weight)
ensemble_endocrine_cell_prop <- list()
for (i in 1:length(ensemble_endocrine_ref_by_weight)) {
  ref_idx <- ensemble_endocrine_ref_by_weight[i]
  ensemble_endocrine_cell_prop[[i]] <- ensemble_endocrine_dec_res[[i]]$prop.only[ref_idx]
}
names(ensemble_endocrine_cell_prop) <- names(ensemble_endocrine_dec_res)
for (i in 1:length(ensemble_endocrine_cell_prop)){
  ensemble_endocrine_cell_prop[[i]] <- ensemble_endocrine_cell_prop[[i]][[1]]
}


#lawlor_endocrine_dec_res <- readRDS("~/Praktikum/decon_res/lawlor_endocrine_dec_res.RDS")
lawlor_endocrine_dec_res <- readRDS("~/Praktikum/decon_res2/lawlor_endocrine_dec_res2.RDS")
lawlor_endocrine_cell_prop <- list()
for (i in 1:length(lawlor_endocrine_dec_res)) {
  lawlor_endocrine_cell_prop[[i]] <- lawlor_endocrine_dec_res[[i]]$prop.est.mvw
}
names(lawlor_endocrine_cell_prop) <- names(lawlor_endocrine_dec_res)


#segerstolpe_endocrine_dec_res <- readRDS("~/Praktikum/decon_res/segerstolpe_endocrine_dec_res.RDS")
segerstolpe_endocrine_dec_res <- readRDS("~/Praktikum/decon_res2/segerstolpe_endocrine_dec_res2.RDS")
segerstolpe_endocrine_cell_prop <- list()
for (i in 1:length(segerstolpe_endocrine_dec_res)) {
  segerstolpe_endocrine_cell_prop[[i]] <- segerstolpe_endocrine_dec_res[[i]]$prop.est.mvw
}
names(segerstolpe_endocrine_cell_prop) <- names(segerstolpe_endocrine_dec_res)


cell_type_prop_endocrine <- list(ensemble_endocrine_cell_prop, baron_endocrine_cell_prop, 
                                 segerstolpe_endocrine_cell_prop, lawlor_endocrine_cell_prop)
names(cell_type_prop_endocrine) <- c("Ensemble", "Baron", "Segerstolpe", "Lawlor")



### cell type predictions ~ MKi67 ###
endocrine_corr_ki67 <- matrix(data = NA, nrow = length(ki67_values)*4, ncol = 10)
colnames(endocrine_corr_ki67) <- c("Bechmark_Data", "Algorithm_Ref", "alpha", "alpha_pvalue", 
                                   "beta", "beta_pvalue", "delta", "delta_pvalue", "gamma", "gamma_pvalue")
endocrine_corr_ki67[,1] <- as.character(sapply(1:length(ki67_values), function(x) rep(names(ki67_values)[x], 4)))
endocrine_corr_ki67[,2] <- rep(c("Ensemble", "Baron", "Segerstolpe", "Lawlor"), length(ki67_values))

for (i_bench in 1:length(ki67_values)) {
  for (i_alg in 1:length(cell_type_prop_endocrine)) {
    off_set <- rnorm(nrow(cell_type_prop_endocrine[[i_alg]][[i_bench]]),mean=0.001,sd=0.001)
    cor_all_celltypes <- cor(cell_type_prop_endocrine[[i_alg]][[i_bench]] + off_set, ki67_values[[i_bench]]) 
    cor_all_celltypes_pval <- sapply(1:4, function(x) cor.test(cell_type_prop_endocrine[[i_alg]][[i_bench]][,x]+off_set, 
                                                               ki67_values[[i_bench]])$p.value)
    
    row_matrix <- i_bench*4-(4-i_alg)
    row_matrix_value <- c(rbind(t(cor_all_celltypes), cor_all_celltypes_pval))
    endocrine_corr_ki67[row_matrix, 3:10] <- row_matrix_value
    
  }
}



### cell type predictions ~ Grading ###
cell_type_prop_endocrine_grading <- lapply(cell_type_prop_endocrine, function(x) x[-1])
cell_type_prop_endocrine_grading <- lapply(cell_type_prop_endocrine_grading, function(x) x[-1])

endocrine_corr_grading <- matrix(data = NA, nrow = length(grading_info)*4, ncol = 10)
colnames(endocrine_corr_grading) <- c("Bechmark_Data", "Algorithm_Ref", "alpha", "alpha_pvalue", 
                                   "beta", "beta_pvalue", "delta", "delta_pvalue", "gamma", "gamma_pvalue")
endocrine_corr_grading[,1] <- as.character(sapply(1:length(grading_info), function(x) rep(names(grading_info)[x], 4)))
endocrine_corr_grading[,2] <- rep(c("Ensemble", "Baron", "Segerstolpe", "Lawlor"), length(grading_info))

for (i_bench in 1:length(grading_info_numeric)) {
  for (i_alg in 1:length(cell_type_prop_endocrine_grading)) {
    off_set <- rnorm(nrow(cell_type_prop_endocrine_grading[[i_alg]][[i_bench]]),mean=0.001,sd=0.001)
    cor_all_celltypes <- cor(cell_type_prop_endocrine_grading[[i_alg]][[i_bench]] + off_set, grading_info_numeric[[i_bench]]) 
    cor_all_celltypes_pval <- sapply(1:4, function(x) cor.test(cell_type_prop_endocrine_grading[[i_alg]][[i_bench]][,x]+off_set, 
                                                               grading_info_numeric[[i_bench]])$p.value)
    
    row_matrix <- i_bench*4-(4-i_alg)
    row_matrix_value <- c(rbind(t(cor_all_celltypes), cor_all_celltypes_pval))
    endocrine_corr_grading[row_matrix, 3:10] <- row_matrix_value
    
  }
}








#### Endo+Exocrine model ####

## read in deconvolution results
# deconvolution results (cell type proportions per sample) are in
# baron_exocrine_dec_res[[1]]$prop.est.mvw and
# ensemble_exocrine_dec_res[[1]]$prop.only$baron

#baron_exocrine_dec_res <- readRDS("~/Praktikum/decon_res/baron_exocrine_dec_res.RDS")
baron_exocrine_dec_res <- readRDS("~/Praktikum/decon_res2/baron_exocrine_dec_res2.RDS")
baron_exocrine_cell_prop <- list()
for (i in 1:length(baron_exocrine_dec_res)) {
  baron_exocrine_cell_prop[[i]] <- baron_exocrine_dec_res[[i]]$prop.est.mvw
}
names(baron_exocrine_cell_prop) <- names(baron_exocrine_dec_res)


#ensemble_exocrine_dec_res <- readRDS("~/Praktikum/decon_res/ensemble_exocrine_dec_res.RDS")
ensemble_exocrine_dec_res <- readRDS("~/Praktikum/decon_res2/ensemble_exocrine_dec_res2.RDS")
ensemble_exocrine_ref_by_weight <- list()
for (i in 1:length(ensemble_exocrine_dec_res)) {
  weights <- ensemble_exocrine_dec_res[[i]]$w_table[1:5,1:3]
  mean_weights <- sapply(1:3, function(x) mean(weights[,x]))
  ensemble_exocrine_ref_by_weight[[i]] <- which(mean_weights == max(mean_weights))
}
ensemble_exocrine_ref_by_weight <- as.integer(ensemble_exocrine_ref_by_weight)
ensemble_exocrine_cell_prop <- list()
for (i in 1:length(ensemble_exocrine_ref_by_weight)) {
  ref_idx <- ensemble_exocrine_ref_by_weight[i]
  ensemble_exocrine_cell_prop[[i]] <- ensemble_exocrine_dec_res[[i]]$prop.only[ref_idx]
}
names(ensemble_exocrine_cell_prop) <- names(ensemble_exocrine_dec_res)
for (i in 1:length(ensemble_exocrine_cell_prop)){
  ensemble_exocrine_cell_prop[[i]] <- ensemble_exocrine_cell_prop[[i]][[1]]
}


#lawlor_exocrine_dec_res <- readRDS("~/Praktikum/decon_res/lawlor_exocrine_dec_res.RDS")
lawlor_exocrine_dec_res <- readRDS("~/Praktikum/decon_res2/lawlor_exocrine_dec_res2.RDS")
lawlor_exocrine_cell_prop <- list()
for (i in 1:length(lawlor_exocrine_dec_res)) {
  lawlor_exocrine_cell_prop[[i]] <- lawlor_exocrine_dec_res[[i]]$prop.est.mvw
}
names(lawlor_exocrine_cell_prop) <- names(lawlor_exocrine_dec_res)


#segerstolpe_exocrine_dec_res <- readRDS("~/Praktikum/decon_res/segerstolpe_exocrine_dec_res.RDS")
segerstolpe_exocrine_dec_res <- readRDS("~/Praktikum/decon_res2/segerstolpe_exocrine_dec_res2.RDS")
segerstolpe_exocrine_cell_prop <- list()
for (i in 1:length(segerstolpe_exocrine_dec_res)) {
  segerstolpe_exocrine_cell_prop[[i]] <- segerstolpe_exocrine_dec_res[[i]]$prop.est.mvw
}
names(segerstolpe_exocrine_cell_prop) <- names(segerstolpe_exocrine_dec_res)


cell_type_prop_exocrine <- list(ensemble_exocrine_cell_prop, baron_exocrine_cell_prop, 
                                 segerstolpe_exocrine_cell_prop, lawlor_exocrine_cell_prop)
names(cell_type_prop_exocrine) <- c("Ensemble", "Baron", "Segerstolpe", "Lawlor")



## cell type predictions ~ MKi67
exocrine_corr_ki67 <- matrix(data = NA, nrow = length(ki67_values)*4, ncol = 14)
colnames(exocrine_corr_ki67) <- c("Bechmark_Data", "Algorithm_Ref", "alpha", "alpha_pvalue", 
                                   "beta", "beta_pvalue", "delta", "delta_pvalue", "gamma", "gamma_pvalue",
                                   "acinar", "acinar_pvalue", "ductal", "ductal_pvalue")
exocrine_corr_ki67[,1] <- as.character(sapply(1:length(ki67_values), function(x) rep(names(ki67_values)[x], 4)))
exocrine_corr_ki67[,2] <- rep(c("Ensemble", "Baron", "Segerstolpe", "Lawlor"), length(ki67_values))

for (i_bench in 1:length(ki67_values)) {
  for (i_alg in 1:length(cell_type_prop_exocrine)) {
    off_set <- rnorm(nrow(cell_type_prop_exocrine[[i_alg]][[i_bench]]),mean=0.001,sd=0.001)
    cor_all_celltypes <- cor(cell_type_prop_exocrine[[i_alg]][[i_bench]] + off_set, ki67_values[[i_bench]]) 
    cor_all_celltypes_pval <- sapply(1:6, function(x) cor.test(cell_type_prop_exocrine[[i_alg]][[i_bench]][,x]+off_set, 
                                                               ki67_values[[i_bench]])$p.value)
    
    row_matrix <- i_bench*4-(4-i_alg)
    row_matrix_value <- c(rbind(t(cor_all_celltypes), cor_all_celltypes_pval))
    exocrine_corr_ki67[row_matrix, 3:14] <- row_matrix_value
    
  }
}



## cell type predictions ~ Grading
cell_type_prop_exocrine_grading <- lapply(cell_type_prop_exocrine, function(x) x[-1])
cell_type_prop_exocrine_grading <- lapply(cell_type_prop_exocrine_grading, function(x) x[-1])

exocrine_corr_grading <- matrix(data = NA, nrow = length(grading_info)*4, ncol = 14)
colnames(exocrine_corr_grading) <- c("Bechmark_Data", "Algorithm_Ref", "alpha", "alpha_pvalue", 
                                      "beta", "beta_pvalue", "delta", "delta_pvalue", "gamma", "gamma_pvalue",
                                      "acinar", "acinar_pvalue", "ductal", "ductal_pvalue")
exocrine_corr_grading[,1] <- as.character(sapply(1:length(grading_info), function(x) rep(names(grading_info)[x], 4)))
exocrine_corr_grading[,2] <- rep(c("Ensemble", "Baron", "Segerstolpe", "Lawlor"), length(grading_info))

for (i_bench in 1:length(grading_info_numeric)) {
  for (i_alg in 1:length(cell_type_prop_exocrine_grading)) {
    off_set <- rnorm(nrow(cell_type_prop_exocrine_grading[[i_alg]][[i_bench]]),mean=0.001,sd=0.001)
    cor_all_celltypes <- cor(cell_type_prop_exocrine_grading[[i_alg]][[i_bench]] + off_set, grading_info_numeric[[i_bench]]) 
    cor_all_celltypes_pval <- sapply(1:6, function(x) cor.test(cell_type_prop_exocrine_grading[[i_alg]][[i_bench]][,x]+off_set, 
                                                               grading_info_numeric[[i_bench]])$p.value)
    
    row_matrix <- i_bench*4-(4-i_alg)
    row_matrix_value <- c(rbind(t(cor_all_celltypes), cor_all_celltypes_pval))
    exocrine_corr_grading[row_matrix, 3:14] <- row_matrix_value
    
  }
}


##### export all 4 matrices (endo+exo) #####
saveRDS(endocrine_corr_ki67, file = "~/Praktikum/corr_res2/endocrine_cor_ki67_2.RDS")
saveRDS(endocrine_corr_grading, file = "~/Praktikum/corr_res2/endocrine_cor_grading2.RDS")
saveRDS(exocrine_corr_ki67, file = "~/Praktikum/corr_res2/exocrine_cor_ki67_2.RDS")
saveRDS(exocrine_corr_grading, file = "~/Praktikum/corr_res2/exocrine_cor_grading2.RDS")

write.csv(endocrine_corr_ki67, file = "~/Praktikum/corr_res2/endocrine_cor_ki67_2.csv", quote = F, sep = "\t",
          col.names = T)
write.csv(endocrine_corr_grading, file = "~/Praktikum/corr_res2/endocrine_cor_grading2.csv", quote = F, sep = "\t",
          col.names = T)
write.csv(exocrine_corr_ki67, file = "~/Praktikum/corr_res2/exocrine_cor_ki67_2.csv", quote = F, sep = "\t",
          col.names = T)
write.csv(exocrine_corr_grading, file = "~/Praktikum/corr_res2/exocrine_cor_grading2.csv", quote = F, sep = "\t",
          col.names = T)
