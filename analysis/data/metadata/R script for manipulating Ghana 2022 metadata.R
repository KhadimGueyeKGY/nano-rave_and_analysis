########## R script for manipulating metadata for Ghana 2022 samples

## Set wd
setwd("~/Documents/science/malaria/ont/ghana_2023/data/metadata")

## Load packages
library("tidyverse")
library("readxl")

## Import data
meta <- read_xlsx("~/Documents/science/malaria/ont/ghana_2023/data/metadata/ghana22_metadata_2023.xlsx", sheet=1)

## Check out data
head(meta)
colnames(meta)

meta$DNA_conc_qubit_ngul <- as.numeric(meta$DNA_conc_qubit_ngul)
meta$DNA_yield_ng <- as.numeric(meta$DNA_yield_ng)
meta$parasitaemia_200WBC <- as.numeric(meta$parasitaemia_200WBC)


## Apply filters
meta_ed2 <- meta %>%
  filter(sample_location == "Navrongo",
         gDNA_sample_available == "TRUE",
         human_genetics_ethics == "TRUE")
#         sample_duplicate == "FALSE",
#         nanopore_sequence_bool == TRUE)

nrow(meta)
nrow(meta_ed2)

summary(meta_ed2$DNA_conc_qubit_ngul)
summary(meta_ed2$DNA_yield_ng)
summary(meta_ed2$parasitaemia_200WBC)

meta_ed3 <- meta_ed2 %>%
  filter(DNA_conc_qubit_ngul > 1) %>%
  filter(DNA_yield_ng > 100) %>%
  filter(parasitaemia_200WBC > 70) # Picked 70 because it gives N=82 samples. Plus + and - control = 84 samples.

nrow(meta_ed2)
nrow(meta_ed3)

# Remove duplicates
meta_ed4 <- meta_ed3 %>%
  arrange(desc(DNA_yield_ng)) %>%
  distinct(patient_id, .keep_all=TRUE)

nrow(meta_ed3)
nrow(meta_ed4)

msp1_sample_select = meta_ed4$sample_id


## Selecting samples for WGS
meta_ed5 <- meta_ed2 %>% # meta if want to include Accra samples
  filter(gDNA_sample_available == "TRUE") %>%
  filter(DNA_conc_qubit_ngul > 1) %>%
  filter(DNA_yield_ng > 500) %>%
  filter(parasitaemia_200WBC > 500)

# Remove duplicates
meta_ed6 <- meta_ed5 %>%
  arrange(desc(DNA_yield_ng)) %>%
  distinct(patient_id, .keep_all=TRUE)

nrow(meta)
nrow(meta_ed5)
nrow(meta_ed6)

# Randomly down-sample to N=8 samples
meta_ed7 <- sample_n(meta_ed6, 50)

meta_ed7
nrow(meta_ed7)

summary(meta_ed7$parasitaemia_200WBC)
summary(meta_ed7$DNA_yield_ng)

wgs_sample_select = meta_ed7$sample_id

## Construct dataframes selecting samples
msp1_select_df <- data.frame(msp1_sample_select) %>%
  rename(sample_id = msp1_sample_select)
msp1_select_df$msp1_select = "TRUE"
head(msp1_select_df)

wgs_select_df <- data.frame(wgs_sample_select) %>%
  rename(sample_id = wgs_sample_select)
wgs_select_df$wgs_select = "TRUE"
head(wgs_select_df)

# merge
meta_merg1 <- merge(meta, msp1_select_df,
                    by="sample_id",
                    all.x=T, all.y=T)
dim(meta_merg1)
dim(meta)

meta_merg2 <- merge(meta_merg1, wgs_select_df,
                    by="sample_id",
                    all.x=T, all.y=T)
dim(meta_merg1)
dim(meta_merg2)

# If it's not true, it's false (rather than NA)
meta_merg2$msp1_select <- as.character(meta_merg2$msp1_select)
meta_merg2$wgs_select <- as.character(meta_merg2$wgs_select)

vars_to_replace <- c("msp1_select", "wgs_select")
meta_temp <- meta_merg2[vars_to_replace]
meta_temp[is.na(meta_temp)] <- FALSE
meta_merg2[vars_to_replace] <- meta_temp

## Save file
write.csv(meta_merg2, "ghana22_metadata_sampleselect_2023.csv", row.names=F)





