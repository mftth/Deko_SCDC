# Dieses Skript dient zum Erstellen von Grafiken, die im Manuskript zu finden sind
# Figure 3 (downloaded manuskript)
# Figure 4
# Figure 5

library(ggplot2)
library(tidyverse)
library(tidyr)
library(reshape)
library(ggpubr)

#### load Metadata file ####
metadata <- read.table(file = "~/Praktikum/Data/Clinial_data_Meta_information.tsv", sep = "\t", 
                       header = T, stringsAsFactors = F)



#### load Benchmark data sets ####
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

grading_repset <- repset_meta$Grading

#### Endocrine only model ####
#### read in deconvolution results ####
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


cell_type_prop_endocrine_repset <- list()
for (i_alg in 1:length(cell_type_prop_endocrine)) {
  cell_type_prop_endocrine_repset[[i_alg]] <- cell_type_prop_endocrine[[i_alg]]$RepSet
}
names(cell_type_prop_endocrine_repset) <- names(cell_type_prop_endocrine)



#### Endo+Exocrine model ####

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


cell_type_prop_exocrine_repset <- list()
for (i_alg in 1:length(cell_type_prop_exocrine)) {
  cell_type_prop_exocrine_repset[[i_alg]] <- cell_type_prop_exocrine[[i_alg]]$RepSet
}
names(cell_type_prop_exocrine_repset) <- names(cell_type_prop_exocrine)



###### Figure 3 ######
## Bar chart
## two next to each other (endo model + endo-exo model)
## do that for each reference data set (ensemble, baron, segerstolpe, lawlor)
## only RepSet is of interest

# add Gradings to celltype proportions
#for (i_alg in 1:length(cell_type_prop_endocrine_repset)) {
#  cell_type_prop_endocrine_repset[[i_alg]] <- cbind(cell_type_prop_endocrine_repset[[i_alg]], 
#                                                    as.factor(grading_repset))
#  cell_type_prop_exocrine_repset[[i_alg]] <- cbind(cell_type_prop_exocrine_repset[[i_alg]], 
#                                                  as.factor(grading_repset))
#}

require(gridExtra)

#melt_endo_ensemble <- melt(cell_type_prop_endocrine_repset$Ensemble)
endocrine_molten <- lapply(cell_type_prop_endocrine_repset, function(x) melt(x))
exocrine_molten <- lapply(cell_type_prop_exocrine_repset, function(x) melt(x))

endocrine_molten <- lapply(endocrine_molten, function(x) cbind(x, rep(grading_repset, 4)))
exocrine_molten <- lapply(exocrine_molten, function(x) cbind(x, rep(grading_repset, 6)))

for (i_alg in 1:length(endocrine_molten)) {
  colnames(endocrine_molten[[i_alg]]) <- c("Sample", "Cell_type", "Value", "Grading")
  colnames(exocrine_molten[[i_alg]]) <- c("Sample", "Cell_type", "Value", "Grading")
  
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}



barplotgrading_endo <- ggplot(endocrine_molten$Lawlor, aes(fill=Cell_type, y=Value, x=Grading)) +
  geom_bar(position = "fill", stat = "identity") +
  ggtitle("Endocrine model") + ylab("cell type proportion") +
  scale_fill_manual(values=c(alpha = "#48d3a9", beta = "#9f4fd7", delta = "#438c4d", gamma = 
                               "#4f7fd7")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.title.x = 
          element_blank(),
        axis.title.y = element_text(size = 13))

barplotgrading_exo <-ggplot(exocrine_molten$Lawlor, aes(fill=Cell_type, y=Value, x=Grading)) +
  geom_bar(position = "fill", stat = "identity") + ggtitle("Endocrine & Exocrine model") +
  scale_fill_manual(values=c(alpha = "#48d3a9", beta = "#9f4fd7", delta = "#438c4d", gamma = 
                               "#4f7fd7", 
                             acinar = "#f35b39", ductal = "#f3ce39"), name = "Cell type", ) +
  theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = 
          element_text(size = 12), axis.title.y = element_blank(), plot.title = 
          element_text(hjust = 0.5), axis.title.x =  element_blank())

legend_barplotgrading<-g_legend(barplotgrading_exo)

barplotgrading <-arrangeGrob(barplotgrading_endo, barplotgrading_exo + theme(legend.position = 
                                                                               "none"), nrow=1)
barplotgrading_xaxis <- arrangeGrob(text_grob("Grading per deconvolution model", hjust = 0.5, 
                                              vjust = 0.3, size = 13))


grid.arrange(barplotgrading, legend_barplotgrading, barplotgrading_xaxis,
             nrow=2, widths = c(9,2), heights=c(10,0.2))  



###### Figure 4 ######
## Bar chart of ML-RF results (accuracy, sensitivity, Specificity)
## do that for each reference data set (ensemble, baron, segerstolpe, lawlor)
## do that for endo & endo + exo models
## only RepSet is of interest

## endo model ##
endocrine_rf_perf <- readRDS("~/Praktikum/RF_res/endocrine_RF_performance2.RDS")
colnames(endocrine_rf_perf) <- c("Algorithm_Ref", "Class", "Sensitivity", "Specificity", 
                                 "Pos_Pred_Value", "Neg_Pred_Value", "Precision", "Recall", "F1", 
                                 "Prevalence", "Detection_Rate", "Detection_Prevalence", 
                                 "Balanced_Accuracy")

endo_rf_perf_data <- data.frame(endocrine_rf_perf[10:12,])
endo_rf_perf_data[,3:ncol(endo_rf_perf_data)] <- sapply(endo_rf_perf_data[,3:ncol(endo_rf_perf_data)], 
                                                        function(x) as.numeric(as.character(x))*100)


barplot_rf_endo_acc <- ggplot(data = endo_rf_perf_data, aes(x=Class, y=Balanced_Accuracy)) + 
  geom_bar(stat = "identity") + ggtitle("Accuracy") + ylab("percantage %") + ylim(0,100) +
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))

barplot_rf_endo_sens <- ggplot(data = endo_rf_perf_data, aes(x=Class, y=Sensitivity)) + 
  geom_bar(stat = "identity") + ggtitle("Sensitivity") + ylim(0,100) +
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

barplot_rf_endo_spec <- ggplot(data = endo_rf_perf_data, aes(x=Class, y=Specificity)) + 
  geom_bar(stat = "identity") + ggtitle("Specificity") + ylim(0,100) +
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

ggarrange(barplot_rf_endo_acc, barplot_rf_endo_sens, barplot_rf_endo_spec, nrow = 1)


## exo model ##
exocrine_rf_perf <- readRDS("~/Praktikum/RF_res/exocrine_RF_performance2.RDS")
colnames(exocrine_rf_perf) <- c("Model", "Class", "Sensitivity", "Specificity", 
                                "Pos_Pred_Value", "Neg_Pred_Value", "Precision", "Recall", "F1", 
                                "Prevalence", "Detection_Rate", "Detection_Prevalence", 
                                "Balanced_Accuracy")

exo_rf_perf_data <- data.frame(exocrine_rf_perf[1:3,])
exo_rf_perf_data[,3:ncol(exo_rf_perf_data)] <- sapply(exo_rf_perf_data[,3:ncol(exo_rf_perf_data)], 
                                                        function(x) as.numeric(as.character(x))*100)


barplot_rf_exo_acc <- ggplot(data = exo_rf_perf_data, aes(x=Class, y=Balanced_Accuracy)) + 
  geom_bar(stat = "identity") + ggtitle("Accuracy") + ylab("percantage %") + ylim(0,100) +
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))

barplot_rf_exo_sens <- ggplot(data = exo_rf_perf_data, aes(x=Class, y=Sensitivity)) + 
  geom_bar(stat = "identity") + ggtitle("Sensitivity") + ylim(0,100) +
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

barplot_rf_exo_spec <- ggplot(data = exo_rf_perf_data, aes(x=Class, y=Specificity)) + 
  geom_bar(stat = "identity") + ggtitle("Specificity") + ylim(0,100) +
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

ggarrange(barplot_rf_exo_acc, barplot_rf_exo_sens, barplot_rf_exo_spec, nrow = 1)


## add values from ki67 model from manuscript ##
ki67_acc <- c("G1"=94, "G2"=87, "G3"=92)
ki67_sens <- c("G1"=93, "G2"=86, "G3"=85)
ki67_spec <- c("G1"=95, "G2"=87, "G3"=100)

ki67_rf_perf_data <- data.frame(ki67_sens, ki67_spec, ki67_acc)
ki67_rf_perf_data <- cbind(rownames(ki67_rf_perf_data), ki67_rf_perf_data)
rownames(ki67_rf_perf_data) <- NULL
ki67_rf_perf_data <- cbind(rep("Ki67",3), ki67_rf_perf_data)
colnames(ki67_rf_perf_data) <- c("Model", "Class", "Sensitivity", "Specificity", 
                                 "Balanced_Accuracy")

# add other unknown values so i can rbind those with exo_rf_perf_data
ki67_unkown <- data.frame(matrix(NA, nrow = 3, ncol = 8))
colnames(ki67_unkown) <- c("Pos_Pred_Value", "Neg_Pred_Value", "Precision", "Recall", "F1", 
                           "Prevalence", "Detection_Rate", "Detection_Prevalence")

ki67_rf_perf_data <- cbind(ki67_rf_perf_data, ki67_unkown)
ki67_rf_perf_data <- ki67_rf_perf_data[,c("Model", "Class", "Sensitivity", "Specificity", 
                                          "Pos_Pred_Value", "Neg_Pred_Value", "Precision", 
                                          "Recall", "F1", "Prevalence", "Detection_Rate", 
                                          "Detection_Prevalence", "Balanced_Accuracy")]

combi_rf_perf_data <- rbind(exo_rf_perf_data, ki67_rf_perf_data)


barplot_rf_combi_acc <- ggplot(data = combi_rf_perf_data, aes(x=Class, y=Balanced_Accuracy, 
                                                              fill=Model)) + 
  geom_bar(stat = "identity", position = "dodge") + ggtitle("Accuracy") + ylab("percantage %") + 
  ylim(0,100) +
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))

barplot_rf_combi_sens <- ggplot(data = combi_rf_perf_data, aes(x=Class, y=Sensitivity, 
                                                               fill=Model)) +
  geom_bar(stat = "identity", position = "dodge") + ggtitle("Sensitivity") + ylim(0,100) +
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

barplot_rf_combi_spec <- ggplot(data = combi_rf_perf_data, aes(x=Class, y=Specificity, 
                                                               fill=Model)) +
  geom_bar(stat = "identity", position = "dodge") + ggtitle("Specificity") + ylim(0,100) +
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

ggarrange(barplot_rf_combi_acc, barplot_rf_combi_sens, barplot_rf_combi_spec, nrow = 1, 
          common.legend = T)
