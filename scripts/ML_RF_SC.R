# Dieses Skript dient zum Erstellen eines Random Forest Machine Learning Modells
# In diesem Skript wird ebenfalls eine Survival-Kurve erstellt
# Features: reconstruction error of deconvolution + derived cell type proportions
# Classification objective: classify panNEN and non-pancreatic GEP-NEN with respect to their grading(, NEC or NET subtype)
# perform cross-validation

# for each model (algorithm+benchmark dataset) create one RF-model
# mostly concentrate on RepSet


#### Random Forest ####

library("dplyr")
library(pROC)
library("caret")
library("stringr")
library("e1071")
set.seed(1)

fitControl <- trainControl(
  method = "cv",
  number = 10,
  sampling = "down",
  savePred=T
)

pred_data = function(
  train_mat,
  truth_vec,
  model = ""
){
  if ( model == ""){
    model = caret::train(
      x = train_mat,
      y = truth_vec,
      method = "rf",
      norm.votes=T,
      #predict.all=FALSE,
      type = "Classification",
      metric= "Accuracy",
      ntree = 500,
      trControl = fitControl)
  }
  
  truth_vec = factor(truth_vec, levels = c("G1","G2","G3"))
  prediction_ml = predict(model, train_mat)
  con_mat = confusionMatrix(
    prediction_ml,
    truth_vec,
    positive = "G3"
  )
  
  return(con_mat)
}



#### load Metadata file ####
metadata <- read.table(file = "~/Praktikum/Data/Clinial_data_Meta_information.tsv", sep = "\t", header = T, stringsAsFactors = F)



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
#missiaglia_sample_idx <- which(missiaglia_meta$Name %in% colnames(missiaglia))
missiaglia_meta$Name <- as.character(missiaglia_meta$Name)
missiaglia_sample_idx <-  match(colnames(missiaglia), missiaglia_meta$Name)
missiaglia_meta <- missiaglia_meta[missiaglia_sample_idx,]

# Riemer
riemer <- read.table(file = "~/Praktikum/Data/Riemer.S40.tsv", sep = '\t', header = TRUE)
rownames(riemer) <- riemer[,1]
riemer <- riemer[,-1]
riemer_meta_idx <- which(metadata$Study == "Riemer") 
riemer_meta <- metadata[riemer_meta_idx,]
#riemer_meta <- riemer_meta[which(riemer_meta$Subtype=="Cancer"),] #40!
riemer_meta$Name <- as.character(riemer_meta$Name)
riemer_meta$Name[1:32] <- paste("X", riemer_meta$Name[1:32], sep = "")
riemer_meta <- riemer_meta[match(colnames(riemer), riemer_meta$Name),]

# Sadanandam
sadanandam <- read.table(file = "~/Praktikum/Data/Sadanandam.S29.tsv", sep = '\t', header = TRUE)
sad_meta_idx <- which(metadata$Study == "Sadanandam") 
sad_meta <- metadata[sad_meta_idx,]

# Scarpa
scarpa <- read.table(file = "~/Praktikum/Data/Scarpa.S29.tsv", sep = '\t', header = TRUE)
scarpa_meta_idx <- which(metadata$Study == "Scarpa") 
scarpa_meta <- metadata[scarpa_meta_idx,]
scarpa_meta$Name <- as.character(scarpa_meta$Name)
scarpa_meta$Name <- paste("X", scarpa_meta$Name, sep = "")
scarpa_meta <- scarpa_meta[match(colnames(scarpa), scarpa_meta$Name),]

# RepSet
repset <- read.table(file = "~/Praktikum/Data/RepSet.S57.tsv", sep = '\t', header = TRUE)
repset_meta <- rbind(scarpa_meta, riemer_meta)
#repset_meta$Name <- as.character(repset_meta$Name)
#repset_meta$Name[1:55] <- paste("X", repset_meta$Name[1:55], sep = "")
rownames(repset_meta) <- repset_meta$Name
repset_meta_idx <- match(colnames(repset), rownames(repset_meta))
repset_meta <- repset_meta[repset_meta_idx,]



#### Grading info ####
#grading_fadista <- fadista_meta$Grading
#names(grading_fadista) <- fadista_meta$Name

grading_missiaglia <- missiaglia_meta$Grading 
grading_repset <- repset_meta$Grading
grading_riemer <- riemer_meta$Grading
grading_sadanandam <- sad_meta$Grading
grading_scarpa <- scarpa_meta$Grading

grading_info <- list(grading_missiaglia, grading_repset, grading_riemer, grading_sadanandam, grading_scarpa)
names(grading_info) <- c("Missiaglia", "RepSet", "Riemer", "Sadanandam", "Scarpa")

#grading_info_numeric <- sapply(1:length(grading_info), function(x) as.vector(grading_info[[x]]))
grading_info_numeric <- sapply(1:length(grading_info), function(x) 
  as.integer(str_replace_all(grading_info_numeric[[x]],pattern ="G","")))
names(grading_info_numeric) <- names(grading_info)

truth_vec <- repset_meta[rownames(train_mat), "Grading"]



#### Endocrine only model ####

## read in deconvolution results ##
# deconvolution results (cell type proportions per sample) are in
# baron_endocrine_dec_res[[1]]$prop.est.mvw and
# ensemble_endocrine_dec_res[[1]]$prop.only$baron

#baron_endocrine_dec_res <- readRDS("~/Praktikum/decon_res/baron_endocrine_dec_res.RDS")
baron_endocrine_dec_res <- readRDS("~/Praktikum/decon_res3/baron_endocrine_dec_res3.RDS")
baron_endocrine_cell_prop <- list()
for (i in 1:length(baron_endocrine_dec_res)) {
  baron_endocrine_cell_prop[[i]] <- baron_endocrine_dec_res[[i]]$prop.est.mvw
}
names(baron_endocrine_cell_prop) <- names(baron_endocrine_dec_res)


#ensemble_endocrine_dec_res <- readRDS("~/Praktikum/decon_res/ensemble_endocrine_dec_res.RDS")
ensemble_endocrine_dec_res <- readRDS("~/Praktikum/decon_res3/ensemble_endocrine_dec_res3.RDS")
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
lawlor_endocrine_dec_res <- readRDS("~/Praktikum/decon_res3/lawlor_endocrine_dec_res3.RDS")
lawlor_endocrine_cell_prop <- list()
for (i in 1:length(lawlor_endocrine_dec_res)) {
  lawlor_endocrine_cell_prop[[i]] <- lawlor_endocrine_dec_res[[i]]$prop.est.mvw
}
names(lawlor_endocrine_cell_prop) <- names(lawlor_endocrine_dec_res)


#segerstolpe_endocrine_dec_res <- readRDS("~/Praktikum/decon_res/segerstolpe_endocrine_dec_res.RDS")
segerstolpe_endocrine_dec_res <- readRDS("~/Praktikum/decon_res3/segerstolpe_endocrine_dec_res3.RDS")
segerstolpe_endocrine_cell_prop <- list()
for (i in 1:length(segerstolpe_endocrine_dec_res)) {
  segerstolpe_endocrine_cell_prop[[i]] <- segerstolpe_endocrine_dec_res[[i]]$prop.est.mvw
}
names(segerstolpe_endocrine_cell_prop) <- names(segerstolpe_endocrine_dec_res)



cell_type_prop_endocrine <- list(ensemble_endocrine_cell_prop, baron_endocrine_cell_prop, 
                                 segerstolpe_endocrine_cell_prop, lawlor_endocrine_cell_prop)
names(cell_type_prop_endocrine) <- c("Ensemble", "Baron", "Segerstolpe", "Lawlor")



#### prepare features  for endocrine model ####
## start with RepSet of every algorithm
## add ensemble_endocrine_dec_res$RepSet$prop.list$lawlor$yeval$RMSDy.sample.table to cell type predictions

# only RepSet cell type proportion predictions
cell_type_prop_endocrine_repset <- list()
for (i_alg in 1:length(cell_type_prop_endocrine)) {
  cell_type_prop_endocrine_repset[[i_alg]] <- cell_type_prop_endocrine[[i_alg]]$RepSet
}
names(cell_type_prop_endocrine_repset) <- names(cell_type_prop_endocrine)

# add the RMSD values from Ensemble algorithm to cell type proportions
repset_weight_ref <- ensemble_endocrine_ref_by_weight[which(names(ensemble_endocrine_dec_res)=="RepSet")]
cell_type_prop_endocrine_repset$Ensemble <- cbind(cell_type_prop_endocrine_repset$Ensemble,
                                                  t(ensemble_endocrine_dec_res$RepSet$prop.list[[repset_weight_ref]]$
                                                      yeval$RMSDy.sample.table),
                                                  t(ensemble_endocrine_dec_res$RepSet$prop.list[[repset_weight_ref]]$
                                                      yeval$mADy.sample.table),
                                                  t(ensemble_endocrine_dec_res$RepSet$prop.list[[repset_weight_ref]]$
                                                      yeval$spearmany.sample.table))
cell_type_prop_endocrine_repset$Baron <- cbind(cell_type_prop_endocrine_repset$Baron,
                                               t(baron_endocrine_dec_res$RepSet$yeval$RMSDy.sample.table),
                                               t(baron_endocrine_dec_res$RepSet$yeval$mADy.sample.table),
                                               t(baron_endocrine_dec_res$RepSet$yeval$spearmany.sample.table))
cell_type_prop_endocrine_repset$Segerstolpe <- cbind(cell_type_prop_endocrine_repset$Segerstolpe,
                                               t(segerstolpe_endocrine_dec_res$RepSet$yeval$RMSDy.sample.table),
                                               t(segerstolpe_endocrine_dec_res$RepSet$yeval$mADy.sample.table),
                                               t(segerstolpe_endocrine_dec_res$RepSet$yeval$spearmany.sample.table))
cell_type_prop_endocrine_repset$Lawlor <- cbind(cell_type_prop_endocrine_repset$Lawlor,
                                               t(lawlor_endocrine_dec_res$RepSet$yeval$RMSDy.sample.table),
                                               t(lawlor_endocrine_dec_res$RepSet$yeval$mADy.sample.table),
                                               t(lawlor_endocrine_dec_res$RepSet$yeval$spearmany.sample.table))



#### perform RF  for endocrine model####
## repeat for all 4 algorithms
res_endo <- list()
d_endo <- list()

for (i_alg in 1:length(cell_type_prop_endocrine_repset)) {
  train_mat <- cell_type_prop_endocrine_repset[[i_alg]]
  colnames(train_mat) <- c("alpha", "beta", "delta", "gamma", "RMSD", "mAD", "spearman")
  res_endo[[i_alg]] <- pred_data(train_mat = train_mat, truth_vec = truth_vec)
  d_endo[[i_alg]] <- res_endo[[i_alg]]$byClass
}

names(res_endo) <- names(cell_type_prop_endocrine_repset)
names(d_endo) <- names(cell_type_prop_endocrine_repset)


endocrine_RF_performance <- matrix(data = NA, nrow = length(d_endo)*nrow(d_endo[[1]]), ncol = 2+ncol(d_endo[[1]]))
colnames(endocrine_RF_performance) <- c("Algorithm_Ref", "Class", colnames(d_endo[[1]]))
endocrine_RF_performance[,1] <- as.character(sapply(1:length(d_endo), function(x) rep(names(d_endo)[x], nrow(d_endo[[1]]))))
endocrine_RF_performance[,2] <- rep(c("G1", "G2", "G3"), length(d_endo))

j = 1
for (i in 1:length(d_endo)) {
  endocrine_RF_performance[j:(j+2), 3:ncol(endocrine_RF_performance)] <- d_endo[[i]]
  j = j+3
}
#endocrine_RF_performance[,-c(1,2)] <- apply(endocrine_RF_performance[,-c(1,2)], 2, as.numeric)


## for feature importance
train_mat <- cell_type_prop_endocrine_repset[[1]]
colnames(train_mat) <- c("alpha", "beta", "delta", "gamma", "RMSD", "mAD", "spearman")

model_endo = caret::train(
  x = train_mat,
  y = truth_vec,
  method = "rf",
  norm.votes=T,
  #predict.all=FALSE,
  type = "Classification",
  metric= "Accuracy",
  ntree = 500,
  trControl = fitControl)

gbmImp_endo <- varImp(model_endo, scale = F) 
plot(varImp(model_endo, scale = F))


# endocrine_RF_feature_importance <- matrix(data = NA, nrow = nrow(gbmImp_endo$importance), ncol = length(d_endo))
# colnames(endocrine_RF_feature_importance) <- names(d_endo)
# rownames(endocrine_RF_feature_importance) <- rownames(gbmImp_endo$importance)
# 
# for (k in 1:length(d_endo)) {
#   train_mat_fi <- cell_type_prop_endocrine_repset[[k]]
#   colnames(train_mat_fi) <- c("alpha", "beta", "delta", "gamma", "RMSD", "mAD", "spearman")
#   model_endo_fi = caret::train(
#     x = train_mat_fi,
#     y = truth_vec,
#     method = "rf",
#     norm.votes=T,
#     #predict.all=FALSE,
#     type = "Classification",
#     metric= "Accuracy",
#     ntree = 500,
#     trControl = fitControl)
#   
#   gbmImp_endo_fi <- varImp(model_endo_fi, scale = T)
#   fi <- gbmImp_endo_fi$importance
#   endocrine_RF_feature_importance[,k] <- fi$Overall
# }


##






#### Exocrine only model ####
## read in deconvolution results
# deconvolution results (cell type proportions per sample) are in
# baron_exocrine_dec_res[[1]]$prop.est.mvw and
# ensemble_exocrine_dec_res[[1]]$prop.only$baron

#baron_exocrine_dec_res <- readRDS("~/Praktikum/decon_res/baron_exocrine_dec_res.RDS")
baron_exocrine_dec_res <- readRDS("~/Praktikum/decon_res3/baron_exocrine_dec_res3.RDS")
baron_exocrine_cell_prop <- list()
for (i in 1:length(baron_exocrine_dec_res)) {
  baron_exocrine_cell_prop[[i]] <- baron_exocrine_dec_res[[i]]$prop.est.mvw
}
names(baron_exocrine_cell_prop) <- names(baron_exocrine_dec_res)


#ensemble_exocrine_dec_res <- readRDS("~/Praktikum/decon_res/ensemble_exocrine_dec_res.RDS")
ensemble_exocrine_dec_res <- readRDS("~/Praktikum/decon_res3/ensemble_exocrine_dec_res3.RDS")
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
lawlor_exocrine_dec_res <- readRDS("~/Praktikum/decon_res3/lawlor_exocrine_dec_res3.RDS")
lawlor_exocrine_cell_prop <- list()
for (i in 1:length(lawlor_exocrine_dec_res)) {
  lawlor_exocrine_cell_prop[[i]] <- lawlor_exocrine_dec_res[[i]]$prop.est.mvw
}
names(lawlor_exocrine_cell_prop) <- names(lawlor_exocrine_dec_res)


#segerstolpe_exocrine_dec_res <- readRDS("~/Praktikum/decon_res/segerstolpe_exocrine_dec_res.RDS")
segerstolpe_exocrine_dec_res <- readRDS("~/Praktikum/decon_res3/segerstolpe_exocrine_dec_res3.RDS")
segerstolpe_exocrine_cell_prop <- list()
for (i in 1:length(segerstolpe_exocrine_dec_res)) {
  segerstolpe_exocrine_cell_prop[[i]] <- segerstolpe_exocrine_dec_res[[i]]$prop.est.mvw
}
names(segerstolpe_exocrine_cell_prop) <- names(segerstolpe_exocrine_dec_res)


cell_type_prop_exocrine <- list(ensemble_exocrine_cell_prop, baron_exocrine_cell_prop, 
                                segerstolpe_exocrine_cell_prop, lawlor_exocrine_cell_prop)
names(cell_type_prop_exocrine) <- c("Ensemble", "Baron", "Segerstolpe", "Lawlor")



#### prepare features  for exocrine model ####
## start with RepSet of every algorithm
## add ensemble_exocrine_dec_res$RepSet$prop.list$lawlor$yeval$RMSDy.sample.table to cell type predictions

# only RepSet cell type proportion predictions
cell_type_prop_exocrine_repset <- list()
for (i_alg in 1:length(cell_type_prop_exocrine)) {
  cell_type_prop_exocrine_repset[[i_alg]] <- cell_type_prop_exocrine[[i_alg]]$RepSet
}
names(cell_type_prop_exocrine_repset) <- names(cell_type_prop_exocrine)

# add the RMSD values from Ensemble algorithm to cell type proportions
repset_weight_ref <- ensemble_exocrine_ref_by_weight[which(names(ensemble_exocrine_dec_res)=="RepSet")]
cell_type_prop_exocrine_repset$Ensemble <- cbind(cell_type_prop_exocrine_repset$Ensemble,
                                                  t(ensemble_exocrine_dec_res$RepSet$prop.list[[repset_weight_ref]]$
                                                      yeval$RMSDy.sample.table),
                                                  t(ensemble_exocrine_dec_res$RepSet$prop.list[[repset_weight_ref]]$
                                                      yeval$mADy.sample.table),
                                                 t(ensemble_exocrine_dec_res$RepSet$prop.list[[repset_weight_ref]]$
                                                     yeval$spearmany.sample.table))
cell_type_prop_exocrine_repset$Baron <- cbind(cell_type_prop_exocrine_repset$Baron,
                                              t(baron_exocrine_dec_res$RepSet$yeval$RMSDy.sample.table),
                                              t(baron_exocrine_dec_res$RepSet$yeval$mADy.sample.table), 
                                              t(baron_exocrine_dec_res$RepSet$yeval$spearmany.sample.table))
cell_type_prop_exocrine_repset$Segerstolpe <- cbind(cell_type_prop_exocrine_repset$Segerstolpe,
                                                    t(segerstolpe_exocrine_dec_res$RepSet$yeval$RMSDy.sample.table),
                                                    t(segerstolpe_exocrine_dec_res$RepSet$yeval$mADy.sample.table),
                                                    t(segerstolpe_exocrine_dec_res$RepSet$yeval$spearmany.sample.table))
cell_type_prop_exocrine_repset$Lawlor <- cbind(cell_type_prop_exocrine_repset$Lawlor,
                                              t(lawlor_exocrine_dec_res$RepSet$yeval$RMSDy.sample.table),
                                              t(lawlor_exocrine_dec_res$RepSet$yeval$mADy.sample.table),
                                              t(lawlor_exocrine_dec_res$RepSet$yeval$spearmany.sample.table))



#### perform RF  for exocrine model####
## repeat for all 4 algorithms
res_exo <- list()
d_exo <- list()

for (i_alg in 1:length(cell_type_prop_exocrine_repset)) {
  train_mat <- cell_type_prop_exocrine_repset[[i_alg]]
  colnames(train_mat) <- c("alpha", "beta", "delta", "gamma", "acinar", "ductal", "RMSD", "mAD", "spearman")
  res_exo[[i_alg]] <- pred_data(train_mat = train_mat, truth_vec = truth_vec)
  d_exo[[i_alg]] <- res_exo[[i_alg]]$byClass
}

names(res_exo) <- names(cell_type_prop_exocrine_repset)
names(d_exo) <- names(cell_type_prop_exocrine_repset)

exocrine_RF_performance <- matrix(data = NA, nrow = length(d_exo)*nrow(d_exo[[1]]), ncol = 2+ncol(d_exo[[1]]))
colnames(exocrine_RF_performance) <- c("Algorithm_Ref", "Class", colnames(d_exo[[1]]))
exocrine_RF_performance[,1] <- as.character(sapply(1:length(d_exo), function(x) rep(names(d_exo)[x], nrow(d_exo[[1]]))))
exocrine_RF_performance[,2] <- rep(c("G1", "G2", "G3"), length(d_exo))

j = 1
for (i in 1:length(d_exo)) {
  exocrine_RF_performance[j:(j+2), 3:ncol(exocrine_RF_performance)] <- d_exo[[i]]
  j = j+3
}



## for feature importance
train_mat <- cell_type_prop_exocrine_repset[[4]]
colnames(train_mat) <- c("alpha", "beta", "delta", "gamma", "acinar", "ductal", "RMSD", "mAD", "spearman")
model_exo = caret::train(
  #x = cell_type_prop_exocrine_repset[[1]],
  x = train_mat,
  y = truth_vec,
  method = "rf",
  norm.votes=T,
  #predict.all=FALSE,
  type = "Classification",
  metric= "Accuracy",
  ntree = 500,
  trControl = fitControl)

gbmImp_exo <- varImp(model_exo, scale = T) 
plot(varImp(model_exo))

##


##### export all 4 matrices (endo+exo) #####
saveRDS(endocrine_RF_performance, file = "~/Praktikum/RF_res/endocrine_RF_performance.RDS")
saveRDS(exocrine_RF_performance, file = "~/Praktikum/RF_res/exocrine_RF_performance.RDS")

write.csv(endocrine_RF_performance, file = "~/Praktikum/RF_res/endocrine_RF_performance.csv", quote = F, sep = "\t",
          col.names = T)
write.csv(exocrine_RF_performance, file = "~/Praktikum/RF_res/exocrine_RF_performance.csv", quote = F, sep = "\t",
          col.names = T)
