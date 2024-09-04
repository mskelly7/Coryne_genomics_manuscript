# Corynebacterium WGS Analysis
# Matthew Kelly, MD, MPH 

remove(list=ls())
setwd("G:/My Drive/Manuscript Bots Microbiome #3 - Coryne Genomics") 
set.seed(1234)

version
library(tidyverse)
library(dplyr)
library(plyr)
library(data.table)
library(gridExtra)
library(httr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(DataCombine)
library(ggpubr)
library(readxl)
library(writexl)
library(Rcpp)
library(ggplot2)

excel_sheets("coryne_BGCs_071222.xlsx")
metadata <- read_excel("coryne_BGCs_071222.xlsx", sheet="Metadata", na="NA")
metadata <- arrange(metadata, by=strain_id)
metadata$unknown_strains <- NULL

pks <- read_excel("coryne_BGCs_071222.xlsx", sheet="PKS", na="NA")
pks_ids <- data.frame(table(pks$strain_id))
pks_ids[pks_ids$Freq > 1,]
# Strains contain a maximum of 2 PKS gene clusters
pks$pks_seq <- as.factor(pks$pks_seq)
pks <- arrange(pks, by=pks_seq)
pks <- pks %>% mutate(pks_unique=paste("pks",cumsum(!duplicated(pks_seq)), sep=""))
pks$pks_unique[is.na(pks$pks_seq)] <- NA
pks <- arrange(pks, by=strain_id)
pks <- pks[,c("strain_id", "pks_unique", "pks_type", "pks_seq")]
pks <- pks[order(pks$pks_type),]
sum(!is.na(pks$pks_seq))
length(na.omit(unique(pks$pks_unique)))
pks_final <- pks %>% distinct(strain_id, pks_unique, .keep_all = TRUE) %>%
  group_by(strain_id) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = row, values_from = c(pks_unique, pks_type, pks_seq))
pks_final <- pks_final[,c("strain_id", "pks_unique_1", "pks_type_1", "pks_seq_1", "pks_unique_2", "pks_type_2", "pks_seq_2")]

nrps <- read_excel("coryne_BGCs_071222.xlsx", sheet="NRPS", na="NA")
nrps_ids <- data.frame(table(nrps$strain_id))
nrps_ids[nrps_ids$Freq > 1,]
# Strains contain a maximum of 4 NRPS gene clusters
nrps$nrps_seq <- as.factor(nrps$nrps_seq)
nrps <- arrange(nrps, by=nrps_seq)
nrps <- nrps %>% mutate(nrps_unique=paste("nrps",cumsum(!duplicated(nrps_seq)), sep=""))
nrps$nrps_unique[is.na(nrps$nrps_seq)] <- NA
nrps <- arrange(nrps, by=strain_id)
nrps <- nrps[,c("strain_id", "nrps_unique", "nrps_seq")]
sum(!is.na(nrps$nrps_seq))
length(na.omit(unique(nrps$nrps_unique)))
nrps_final <- nrps %>% distinct(strain_id, nrps_unique, .keep_all = TRUE) %>%
  group_by(strain_id) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = row, values_from = c(nrps_unique, nrps_seq))
nrps_final <- nrps_final[,c("strain_id", "nrps_unique_1", "nrps_seq_1", "nrps_unique_2", "nrps_seq_2", 
                            "nrps_unique_3", "nrps_seq_3", "nrps_unique_4", "nrps_seq_4")]

nrps_like <- read_excel("coryne_BGCs_071222.xlsx", sheet="NRPS-like", na="NA")
nrps_like_ids <- data.frame(table(nrps_like$strain_id))
nrps_like_ids[nrps_like_ids$Freq > 1,]
# No strains contain >1 NRPS-like gene cluster
nrps_like$nrps_like_seq <- as.factor(nrps_like$nrps_like_seq)
nrps_like <- arrange(nrps_like, by=nrps_like_seq)
nrps_like <- nrps_like %>% mutate(nrps_like_unique=paste("nrps_like",cumsum(!duplicated(nrps_like_seq)), sep=""))
nrps_like$nrps_like_unique[is.na(nrps_like$nrps_like_seq)] <- NA
nrps_like <- arrange(nrps_like, by=strain_id)
nrps_like_final <- nrps_like[,c("strain_id", "nrps_like_unique", "nrps_like_seq")]
sum(!is.na(nrps_like_final$nrps_like_seq))
length(na.omit(unique(nrps_like_final$nrps_like_unique)))

napaa <- read_excel("coryne_BGCs_071222.xlsx", sheet="NAPAA", na="NA")
napaa_ids <- data.frame(table(napaa$strain_id))
napaa_ids[napaa_ids$Freq > 1,]
# Strains contain a maximum of 2 NAPAA gene clusters
napaa$napaa_seq <- as.factor(napaa$napaa_seq)
napaa <- arrange(napaa, by=napaa_seq)
napaa <- napaa %>% mutate(napaa_unique=paste("napaa",cumsum(!duplicated(napaa_seq)), sep=""))
napaa$napaa_unique[is.na(napaa$napaa_seq)] <- NA
napaa <- arrange(napaa, by=strain_id)
napaa <- napaa[,c("strain_id", "napaa_unique", "napaa_seq")]
sum(!is.na(napaa$napaa_seq))
length(na.omit(unique(napaa$napaa_unique)))
napaa_final <- napaa %>% distinct(strain_id, napaa_unique, .keep_all = TRUE) %>%
  group_by(strain_id) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = row, values_from = c(napaa_unique, napaa_seq))
napaa_final <- napaa_final[,c("strain_id", "napaa_unique_1", "napaa_seq_1", "napaa_unique_2", "napaa_seq_2")]

terpene <- read_excel("coryne_BGCs_071222.xlsx", sheet="Terpene", na="NA")
terpene_ids <- data.frame(table(terpene$strain_id))
terpene_ids[terpene_ids$Freq > 1,]
# Some strains contain 2 terpene gene clusters
terpene$terpene_seq <- as.factor(terpene$terpene_seq)
terpene <- arrange(terpene, by=terpene_seq)
terpene <- terpene %>% mutate(terpene_unique=paste("terp",cumsum(!duplicated(terpene_seq)), sep=""))
terpene$terpene_unique[is.na(terpene$terpene_seq)] <- NA
terpene <- arrange(terpene, by=strain_id)
terpene <- terpene[,c("strain_id", "terpene_unique", "terpene_seq")]
sum(!is.na(terpene$terpene_seq))
length(na.omit(unique(terpene$terpene_unique)))
terpene_final <- terpene %>% distinct(strain_id, terpene_unique, .keep_all = TRUE) %>%
  group_by(strain_id) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = row, values_from = c(terpene_unique, terpene_seq))
terpene_final <- terpene_final[,c("strain_id", "terpene_unique_1", "terpene_seq_1", "terpene_unique_2", "terpene_seq_2")]

siderophore <- read_excel("coryne_BGCs_071222.xlsx", sheet="Siderophore", na="NA")
siderophore_ids <- data.frame(table(siderophore$strain_id))
siderophore_ids[siderophore_ids$Freq > 1,]
# Some strains contain 2 siderophore gene clusters
siderophore$siderophore_seq <- as.factor(siderophore$siderophore_seq)
siderophore <- arrange(siderophore, by=siderophore_seq)
siderophore <- siderophore %>% mutate(siderophore_unique=paste("sid",cumsum(!duplicated(siderophore_seq)), sep=""))
siderophore$siderophore_unique[is.na(siderophore$siderophore_seq)] <- NA
siderophore <- arrange(siderophore, by=strain_id)
siderophore <- siderophore[,c("strain_id", "siderophore_unique", "siderophore_seq")]
sum(!is.na(siderophore$siderophore_seq))
length(na.omit(unique(siderophore$siderophore_unique)))
siderophore_final <- siderophore %>% distinct(strain_id, siderophore_unique, .keep_all = TRUE) %>%
  group_by(strain_id) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = row, values_from = c(siderophore_unique, siderophore_seq))
siderophore_final <- siderophore_final[,c("strain_id", "siderophore_unique_1", "siderophore_seq_1", "siderophore_unique_2", "siderophore_seq_2")]

ripp <- read_excel("coryne_BGCs_071222.xlsx", sheet="RiPP", na="NA")
ripp_ids <- data.frame(table(ripp$strain_id))
ripp_ids[ripp_ids$Freq > 1,]
# Some strains contain up to 3 bacteriocin gene clusters
ripp$ripp_seq <- as.factor(ripp$ripp_seq)
ripp <- arrange(ripp, by=ripp_seq)
ripp <- ripp %>% mutate(ripp_unique=paste("ripp",cumsum(!duplicated(ripp_seq)), sep=""))
ripp$ripp_unique[is.na(ripp$ripp_seq)] <- NA
ripp <- arrange(ripp, by=strain_id)
ripp <- ripp[,c("strain_id", "ripp_unique", "ripp_type", "ripp_seq_type", "ripp_seq")]
sum(!is.na(ripp$ripp_seq))
length(na.omit(unique(ripp$ripp_unique)))
ripp_final <- ripp %>% distinct(strain_id, ripp_unique, .keep_all = TRUE) %>%
  group_by(strain_id) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = row, values_from = c(ripp_unique, ripp_type, ripp_seq_type, ripp_seq))
ripp_final <- ripp_final[,c("strain_id", "ripp_unique_1", "ripp_type_1", "ripp_seq_type_1", "ripp_seq_1", "ripp_unique_2", "ripp_type_2", "ripp_seq_type_2", 
                            "ripp_seq_2", "ripp_unique_3", "ripp_type_3", "ripp_seq_type_3", "ripp_seq_3")]

other <- read_excel("coryne_BGCs_071222.xlsx", sheet="Other", na="NA")
other_ids <- data.frame(table(other$strain_id))
other_ids[other_ids$Freq > 1,]
# No strains contain >1 other gene cluster
other$other_seq <- as.factor(other$other_seq)
other <- arrange(other, by=other_seq)
other <- other %>% mutate(other_unique=paste("oth",cumsum(!duplicated(other_seq)), sep=""))
other$other_unique[is.na(other$other_seq)] <- NA
other <- arrange(other, by=strain_id)
other <- other[,c("strain_id", "other_type", "other_unique", "other_seq")]
table(other$other_type)
sum(!is.na(other$other_seq))
length(na.omit(unique(other$other_unique)))
betalactone <- subset(other, other_type=="Betalactone")
names(betalactone)[names(betalactone)=="other_unique"] <- "beta_unique"
betalactone$beta_unique <- gsub("oth", "beta", as.factor(betalactone$beta_unique))
names(betalactone)[names(betalactone)=="other_seq"] <- "beta_seq"
betalactone$other_type <- NULL
cytolysin <- subset(other, other_type=="Cytolysin")
names(cytolysin)[names(cytolysin)=="other_unique"] <- "cyto_unique"
cytolysin$cyto_unique <- gsub("oth", "cyto", as.factor(cytolysin$cyto_unique))
names(cytolysin)[names(cytolysin)=="other_seq"] <- "cyto_seq"
cytolysin$other_type <- NULL
ectoine <- subset(other, other_type=="Ectoine")
names(ectoine)[names(ectoine)=="other_unique"] <- "ect_unique"
ectoine$ect_unique <- gsub("oth", "ect", as.factor(ectoine$ect_unique))
names(ectoine)[names(ectoine)=="other_seq"] <- "ect_seq"
ectoine$other_type <- NULL
nucleoside <- subset(other, other_type=="Nucleoside")
names(nucleoside)[names(nucleoside)=="other_unique"] <- "nuc_unique"
nucleoside$nuc_unique <- gsub("oth", "nuc", as.factor(nucleoside$nuc_unique))
names(nucleoside)[names(nucleoside)=="other_seq"] <- "nuc_seq"
nucleoside$other_type <- NULL
phenazine <- subset(other, other_type=="Phenazine")
names(phenazine)[names(phenazine)=="other_unique"] <- "phen_unique"
phenazine$phen_unique <- gsub("oth", "phen", as.factor(phenazine$phen_unique))
names(phenazine)[names(phenazine)=="other_seq"] <- "phen_seq"
phenazine$other_type <- NULL
other_final <- merge(merge(merge(merge(betalactone, cytolysin, by="strain_id", all.x=TRUE, all.y=TRUE), ectoine, by="strain_id", all.x=TRUE, all.y=TRUE),
                           nucleoside, by="strain_id", all.x=TRUE, all.y=TRUE), phenazine, by="strain_id", all.x=TRUE, all.y=TRUE)
other_final$beta_unique[other_final$beta_unique=="beta7"] <- "beta1"
other_final$cyto_unique[other_final$cyto_unique=="cyto5"] <- "cyto1"
other_final$nuc_unique[other_final$nuc_unique=="nuc6"] <- "nuc1"
other_final$phen_unique[other_final$phen_unique=="phen8"] <- "phen1"

remove(pks_ids, nrps_ids, nrps_like_ids, napaa_ids, terpene_ids, siderophore_ids, ripp_ids, other_ids)

bgc_final <- merge(merge(merge(merge(merge(merge(merge(merge(metadata, pks_final, by="strain_id", all.x=FALSE), nrps_final, by="strain_id"), 
                                                             nrps_like_final, by="strain_id"), napaa_final, by="strain_id"), terpene_final, by="strain_id"), 
                                                             siderophore_final, by="strain_id"), ripp_final, by="strain_id"), 
                                                             other_final, by="strain_id", all.x=TRUE)

# *****************************************************
# FOCUS ONLY ON RESPIRATORY STRAINS OF CORYNEBACTERIUM
# *****************************************************

# Summary of all bacterial strains that have been sequenced through this project
nrow(metadata)
table(metadata$location)
# 78 strains from Botswana, 35 strains from Duke patients, 93 strains from the Lemon Laboratory

# Focus analyses on strains isolated from the respiratory tract and exclude strains corresponding to novel species 
metadata_resp <- subset(metadata, specimen_source=="Respiratory" & species!="unknown")
nrow(metadata_resp)
table(metadata_resp$country)
table(metadata_resp$location)
table(metadata_resp$species)

bgc_resp <- subset(bgc_final, specimen_source=="Respiratory" & species!="unknown")
num_pks_unique <- length(na.omit(unique(bgc_resp$pks_unique_1))) + length(na.omit(unique(bgc_resp$pks_unique_2)))
num_nrps_unique <- length(na.omit(unique(bgc_resp$nrps_unique_1))) + length(na.omit(unique(bgc_resp$nrps_unique_2))) + 
                                                                       length(na.omit(unique(bgc_resp$nrps_unique_3)))
num_nrps_like_unique <- length(na.omit(unique(bgc_resp$nrps_like_unique)))
num_napaa_unique <- length(na.omit(unique(bgc_resp$napaa_unique_1))) + length(na.omit(unique(bgc_resp$napaa_unique_2)))
num_terpene_unique <- length(na.omit(unique(bgc_resp$terpene_unique_1))) + length(na.omit(unique(bgc_resp$terpene_unique_2)))
num_siderophore_unique <- length(na.omit(unique(bgc_resp$siderophore_unique_1))) + length(na.omit(unique(bgc_resp$siderophore_unique_2)))
num_ripp_unique <- length(na.omit(unique(bgc_resp$ripp_unique_1))) + length(na.omit(unique(bgc_resp$ripp_unique_2)))
num_betalactone_unique <- length(na.omit(unique(bgc_resp$beta_unique)))
num_cytolysin <- length(na.omit(unique(bgc_resp$cyto_unique)))
num_ectoine_unique <- length(na.omit(unique(bgc_resp$ect_unique)))
num_nucleoside_unique <- length(na.omit(unique(bgc_resp$nuc_unique)))
num_phenazine_unique <- length(na.omit(unique(bgc_resp$phen_unique)))
total_num_bgc <- num_pks + num_nrps + num_nrps_like + num_napaa + num_terpene + num_siderophore + num_ripp + num_linaridin + num_betalactone +
  num_cytolysin + num_ectoine + num_nucleoside + num_phenazine

bgc_resp$pks1_num[!is.na(bgc_resp$pks_unique_1)] <- 1
bgc_resp$pks1_num[is.na(bgc_resp$pks_unique_1)] <- 0
bgc_resp$pks2_num[!is.na(bgc_resp$pks_unique_2)] <- 1
bgc_resp$pks2_num[is.na(bgc_resp$pks_unique_2)] <- 0
bgc_resp$num_pks <- bgc_resp$pks1_num + bgc_resp$pks2_num
bgc_resp$pks1_num <- NULL
bgc_resp$pks2_num <- NULL
table(bgc_resp$num_pks)
bgc_resp$nrps1_num[!is.na(bgc_resp$nrps_unique_1)] <- 1
bgc_resp$nrps1_num[is.na(bgc_resp$nrps_unique_1)] <- 0
bgc_resp$nrps2_num[!is.na(bgc_resp$nrps_unique_2)] <- 1
bgc_resp$nrps2_num[is.na(bgc_resp$nrps_unique_2)] <- 0
bgc_resp$nrps3_num[!is.na(bgc_resp$nrps_unique_3)] <- 1
bgc_resp$nrps3_num[is.na(bgc_resp$nrps_unique_3)] <- 0
bgc_resp$nrps4_num[!is.na(bgc_resp$nrps_unique_4)] <- 1
bgc_resp$nrps4_num[is.na(bgc_resp$nrps_unique_4)] <- 0
bgc_resp$num_nrps <- bgc_resp$nrps1_num + bgc_resp$nrps2_num + bgc_resp$nrps3_num + bgc_resp$nrps4_num
bgc_resp$nrps1_num <- NULL
bgc_resp$nrps2_num <- NULL
bgc_resp$nrps3_num <- NULL
bgc_resp$nrps4_num <- NULL
table(bgc_resp$num_nrps)
bgc_resp$num_nrps_like[!is.na(bgc_resp$nrps_like_unique)] <- 1
bgc_resp$num_nrps_like[is.na(bgc_resp$nrps_like_unique)] <- 0
table(bgc_resp$num_nrps_like)
bgc_resp$napaa1_num[!is.na(bgc_resp$napaa_unique_1)] <- 1
bgc_resp$napaa1_num[is.na(bgc_resp$napaa_unique_1)] <- 0
bgc_resp$napaa2_num[!is.na(bgc_resp$napaa_unique_2)] <- 1
bgc_resp$napaa2_num[is.na(bgc_resp$napaa_unique_2)] <- 0
bgc_resp$num_napaa <- bgc_resp$napaa1_num + bgc_resp$napaa2_num 
bgc_resp$napaa1_num <- NULL
bgc_resp$napaa2_num <- NULL
table(bgc_resp$num_napaa)
bgc_resp$terpene1_num[!is.na(bgc_resp$terpene_unique_1)] <- 1
bgc_resp$terpene1_num[is.na(bgc_resp$terpene_unique_1)] <- 0
bgc_resp$terpene2_num[!is.na(bgc_resp$terpene_unique_2)] <- 1
bgc_resp$terpene2_num[is.na(bgc_resp$terpene_unique_2)] <- 0
bgc_resp$num_terpene <- bgc_resp$terpene1_num + bgc_resp$terpene2_num
bgc_resp$terpene1_num <- NULL
bgc_resp$terpene2_num <- NULL
table(bgc_resp$num_terpene)
bgc_resp$siderophore1_num[!is.na(bgc_resp$siderophore_unique_1)] <- 1
bgc_resp$siderophore1_num[is.na(bgc_resp$siderophore_unique_1)] <- 0
bgc_resp$siderophore2_num[!is.na(bgc_resp$siderophore_unique_2)] <- 1
bgc_resp$siderophore2_num[is.na(bgc_resp$siderophore_unique_2)] <- 0
bgc_resp$num_siderophore <- bgc_resp$siderophore1_num + bgc_resp$siderophore2_num
bgc_resp$siderophore1_num <- NULL
bgc_resp$siderophore2_num <- NULL
table(bgc_resp$num_siderophore)
bgc_resp$ripp1_num[!is.na(bgc_resp$ripp_unique_1)] <- 1
bgc_resp$ripp1_num[is.na(bgc_resp$ripp_unique_1)] <- 0
bgc_resp$ripp2_num[!is.na(bgc_resp$ripp_unique_2)] <- 1
bgc_resp$ripp2_num[is.na(bgc_resp$ripp_unique_2)] <- 0
bgc_resp$ripp3_num[!is.na(bgc_resp$ripp_unique_3)] <- 1
bgc_resp$ripp3_num[is.na(bgc_resp$ripp_unique_3)] <- 0
bgc_resp$num_ripp <- bgc_resp$ripp1_num + bgc_resp$ripp2_num + bgc_resp$ripp3_num
bgc_resp$ripp1_num <- NULL
bgc_resp$ripp2_num <- NULL
bgc_resp$ripp3_num <- NULL
table(bgc_resp$num_ripp)
bgc_resp$num_betalactone[!is.na(bgc_resp$beta_unique)] <- 1
bgc_resp$num_betalactone[is.na(bgc_resp$beta_unique)] <- 0
table(bgc_resp$num_betalactone)
bgc_resp$num_cytolysin[!is.na(bgc_resp$cyto_unique)] <- 1
bgc_resp$num_cytolysin[is.na(bgc_resp$cyto_unique)] <- 0
table(bgc_resp$num_cytolysin)
bgc_resp$num_ectoine[!is.na(bgc_resp$ect_unique)] <- 1
bgc_resp$num_ectoine[is.na(bgc_resp$ect_unique)] <- 0
table(bgc_resp$num_ectoine)
bgc_resp$num_nucleoside[!is.na(bgc_resp$nuc_unique)] <- 1
bgc_resp$num_nucleoside[is.na(bgc_resp$nuc_unique)] <- 0
table(bgc_resp$num_nucleoside)
bgc_resp$num_phenazine[!is.na(bgc_resp$phen_unique)] <- 1
bgc_resp$num_phenazine[is.na(bgc_resp$phen_unique)] <- 0
table(bgc_resp$num_phenazine)

bgc_resp$num_bgc <- bgc_resp$num_pks + bgc_resp$num_nrps + bgc_resp$num_nrps_like + bgc_resp$num_napaa + bgc_resp$num_terpene +
  bgc_resp$num_siderophore + bgc_resp$num_ripp + bgc_resp$num_betalactone + bgc_resp$num_cytolysin + bgc_resp$num_ectoine + bgc_resp$num_nucleoside + 
  bgc_resp$num_phenazine
bgc_resp <- bgc_resp[ , colSums(is.na(bgc_resp)) < nrow(bgc_resp)]

# *****************************************************
# CREATE FIGURE SUMMARIZING BGCS BY SPECIES AND COUNTRY
# *****************************************************

# Summarize the respiratory strains and genomic data contained within this dataset
nrow(bgc_resp)
table(bgc_resp$location)
table(bgc_resp$species)
table(bgc_resp$location, bgc_resp$species)
summary(bgc_resp$ave_coverage, useNA="always")
summary(bgc_resp$completeness, useNA="always")
summary(bgc_resp$contamination, useNA="always")

# Summarize the BGCs among these strains
table(bgc_resp$species)
summary(bgc_resp$num_bgc)
tapply(bgc_resp$num_bgc, bgc_resp$species, mean)

# Summarize BGCs by location
tapply(bgc_resp$num_bgc, bgc_resp$country, mean)
wilcox.test(num_bgc ~ country, data=bgc_resp, exact=FALSE)

# Generate a figure summarizing the mean number of BGCs by Corynebacterium species
# Note that there are no gene clusters for betalactones in respiratory strains

bgc_fig1 <- bgc_resp[,c("species", "country", "num_cytolysin", "num_ectoine", "num_napaa", "num_nrps", "num_nrps_like", "num_nucleoside", 
                        "num_phenazine", "num_pks", "num_ripp", "num_siderophore", "num_terpene")]

bgc_fig1$species_fig <- bgc_fig1$species
bgc_fig1$species_fig[bgc_fig1$species=="accolens"] <- "C. accolens"
bgc_fig1$species_fig[bgc_fig1$species=="kefirresidentii"] <- "C. kefirresidentii"
bgc_fig1$species_fig[bgc_fig1$species=="propinquum"] <- "C. propinquum"
bgc_fig1$species_fig[bgc_fig1$species=="pseudodiphtheriticum"] <- "C. pseudodiphtheriticum"
bgc_fig1$species_fig[bgc_fig1$species=="tuberculostearicum"] <- "C. tuberculostearicum"
bgc_fig1$species_fig[bgc_fig1$species=="amycolatum" | bgc_fig1$species=="appendicis" | bgc_fig1$species=="aurimucosum" | bgc_fig1$species=="bovis" | 
                       bgc_fig1$species=="coyleae" | bgc_fig1$species=="freneyi" | bgc_fig1$species=="striatum" |
                       bgc_fig1$species=="mastitidis"] <- "Other"
table(bgc_fig1$species_fig, useNA="always")

bgc_fig1 <- bgc_fig1[,c("species_fig", "country", "num_cytolysin", "num_ectoine", "num_napaa", "num_nrps","num_nrps_like", "num_nucleoside", 
                        "num_phenazine", "num_pks", "num_ripp", "num_siderophore", "num_terpene")]

bgc_fig1 <- aggregate(bgc_fig1[,3:13], list(bgc_fig1$species_fig, bgc_fig1$country), mean)
bgc_fig1 <- gather(bgc_fig1, num, mean, num_cytolysin:num_terpene, factor_key = TRUE)
names(bgc_fig1)[names(bgc_fig1) == "Group.1"] <- "species_fig"
names(bgc_fig1)[names(bgc_fig1) == "Group.2"] <- "country"

colors_11 <- c("chartreuse", "dodgerblue1", "black", "yellow", "#6A3D9A", "green4", "#CAB2D6", "blue1", "#ff8c00", "brown", "#E31A1C")

bgc_fig1$num <- factor(bgc_fig1$num, levels=c("num_cytolysin", "num_ectoine", "num_napaa", "num_nrps","num_nrps_like", "num_nucleoside", 
                                              "num_phenazine", "num_pks", "num_ripp", "num_siderophore", "num_terpene"), 
                       labels=c("Cytolysin", "Ectoine", "NAPAA", "NRPS", "NRPS-like", "Nucleoside", "Phenazine", "PKS",
                                "RiPP", "Siderophore", "Terpene")) 
                                
bgc_fig1$species_fig <- factor(bgc_fig1$species_fig, levels = c("C. accolens", "C. kefirresidentii", "C. propinquum", 
                                                                "C. pseudodiphtheriticum", "C. tuberculostearicum", "Other"))

table(bgc_resp$country, bgc_resp$species)

fig1 <- ggplot(bgc_fig1, aes(fill=num, y=species_fig, x=mean)) + geom_bar(position="stack", stat="identity") + 
  geom_col(position = position_stack(reverse = TRUE)) +
  theme(legend.text=element_text(size=16), legend.title=element_text(size=18), 
        axis.title.y = element_text(angle=90, size=18), 
        axis.text.y = element_text(size=16, colour="black"), axis.title.x = element_text(size=18), 
        axis.text.x = element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust=1.4), 
        strip.text.x = element_text(size = 16), strip.background=element_rect(fill="white"), panel.background = element_blank()) + 
  scale_fill_manual(values=colors_11) + scale_y_discrete(limits=rev) + xlab("Mean number of BGCs") + ylab("") +
  guides(fill = guide_legend(ncol=1, byrow=TRUE, title=""))

png(file="R_Plots/BGCs_by_species.png", 
    width = 12, height = 8, units = 'in', res = 600)
print(fig1)
dev.off()

fig1_country <- ggplot(bgc_fig1, aes(fill=num, y=species_fig, x=mean)) + geom_bar(position="stack", stat="identity") + 
  geom_col(position = position_stack(reverse = TRUE)) +
  theme(legend.text=element_text(size=16), legend.title=element_text(size=18), 
        axis.title.y = element_text(angle=90, size=18), 
        axis.text.y = element_text(size=16, colour="black"), axis.title.x = element_text(size=18), 
        axis.text.x = element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust=1.4), 
        strip.text.x = element_text(size = 16), strip.background=element_rect(fill="white"), panel.background = element_blank()) + 
  scale_fill_manual(values=colors_11) + scale_y_discrete(limits=rev) + xlab("Mean number of BGCs") + ylab("") +
  guides(fill = guide_legend(ncol=1, byrow=TRUE, title="")) + facet_wrap("country")

png(file="R_Plots/BGCs_by_species_country.png", 
    width = 18, height = 10, units = 'in', res = 600)
print(fig1_country)
dev.off()

# *************************************************************
# EXPORT BGC SEQUENCES AS FASTA FILES FOR PHYLOGENETIC ANALYSES
# *************************************************************

# Creates a simple function that takes a dataframe that has a column "name" and "seq" and writes a fasta file from it

writeFasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))}
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)}

set.seed(1234)

pks_ids <- unique(append(bgc_resp$pks_unique_1, bgc_resp$pks_unique_2))
pks_ids <- pks_ids[!is.na(pks_ids)]
pks_fasta <- subset(pks, pks$pks_unique %in% pks_ids)
pks_fasta <- pks_fasta %>% group_by(pks_unique) %>% sample_n(1)
pks_fasta$name <- gsub(" ", "", paste(pks_fasta$pks_type, "_", pks_fasta$pks_unique, "_", pks_fasta$strain_id))
names(pks_fasta)[names(pks_fasta)=="pks_seq"] <- "seq"
pks_fasta$seq <- as.character(pks_fasta$seq)
pks_fasta <- pks_fasta[,c("name","seq")]
writeFasta(pks_fasta, "BGCs_PKS/pks.fasta")

ripp_ids <- unique(append(append(append(bgc_resp$ripp_unique_1, bgc_resp$ripp_unique_2), bgc_resp$ripp_unique_3), bgc_resp$ripp_unique_4))
ripp_ids <- ripp_ids[!is.na(ripp_ids)]
ripp_fasta <- subset(ripp, ripp$ripp_unique %in% ripp_ids)
ripp_fasta <- ripp_fasta %>% group_by(ripp_unique) %>% sample_n(1)
ripp_fasta$name <- gsub(" ", "", paste(ripp_fasta$ripp_unique, "_", ripp_fasta$strain_id))
names(ripp_fasta)[names(ripp_fasta)=="ripp_seq"] <- "seq"
ripp_fasta$seq <- as.character(ripp_fasta$seq)
ripp_fasta <- subset(ripp_fasta, ripp_type=="bacteriocin_IId")
ripp_fasta <- ripp_fasta[,c("name","seq")]
writeFasta(ripp_fasta, "BGCs_RiPPs/classIId.fasta")

remove(betalactone, ectoine, napaa, napaa_final, nrps, nrps_final, nrps_like, nrps_like_final, nucleoside, other, other_final, phenazine, pks, pks_final,
       ripp, ripp_final, siderophore, siderophore_final, terpene, terpene_final)
remove(pks_fasta, ripp_fasta, pks_ids, ripp_ids, writeFasta)

# *******************************************************
# CREATE PHYLOGENETIC TREES OF GENE CLUSTERS FOR ANALYSIS
# *******************************************************

library(ggtree)

pks_tree <- read.tree("BGCs_PKS/PKS_tree.nwk")
tree <- ggtree(pks_tree)
d <- tree$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d$label <- d$label*100

tree + geom_text(data=d, aes(label=label), hjust = 1.2, vjust = -0.4, size = 4) + theme_tree() + geom_tiplab(align=TRUE, linesize=0.2) + xlim(0,12)

# ************
# INTRODUCTION
# ************

nrow(bgc_resp)
sum(bgc_resp$num_bgc)

# *******
# METHODS
# *******

table(bgc_resp$location)
table(bgc_resp$species)

sum(bgc_resp$num_bgc) 
summary(bgc_resp$num_bgc)

sum(bgc_resp$num_pks) 
sum(bgc_resp$num_nrps) 
sum(bgc_resp$num_napaa) 
sum(bgc_resp$num_terpene) 
sum(bgc_resp$num_siderophore) 
sum(bgc_resp$num_ripp) 
ripps <- bgc_resp[,c("strain_id","species","ripp_type_1", "ripp_type_2")]
ripps <- subset(ripps, ripp_type_1!="" | ripp_type_2!="")
table(ripps$ripp_type_1)
table(ripps$ripp_type_2)

table(bgc_resp$pks_type_1)
table(bgc_resp$pks_type_2)

