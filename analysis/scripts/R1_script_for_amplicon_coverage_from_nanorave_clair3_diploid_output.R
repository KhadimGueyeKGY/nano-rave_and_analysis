##################### R script to read in BED files to check coverage - GAMP_22906B2
### From nano-rave output, Clair3 genotype. Example.

##### Set wd
path_setwd <- getwd()
setwd(path_setwd)

# ---------- Parse command line arguments ----------
multiplex_run = "B2"
#multiplex_seq_date = as.Date("2022-09-06")

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define argument names and default values
arg_names <- c("-m", "-g", "-d", "-a")
arg_values <- setNames(rep(NA, length(arg_names)), arg_names)

# Ensure sufficient arguments are provided
if (length(args) < length(arg_names) * 2) {
  stop("\n❌ Missing required arguments.\n\nUsage:\n  Rscript script.R -m MinION_run_name -g genotyping_method -d analysis_date -a analysis_dr\n\nArguments:\n  -m  MinION run name          (e.g. GAMP_22906B2)\n  -g  Genotyping method        (e.g. clair3_diploid)\n  -d  Analysis date            (format: YYYYMMDD, e.g. 20250511)\n  -a  Path to output directory (e.g. /mnt/c/.../output/)\n")
}

# Parse arguments using vectorized lookup
for (i in seq(1, length(args), by = 2)) {
  if (args[i] %in% arg_names) {
    arg_values[args[i]] <- args[i + 1]
  }
}

# Check if any required arguments are missing
if (any(is.na(arg_values))) {
  stop("\n❌ Missing required arguments.\n\nUsage:\n  Rscript script.R -m MinION_run_name -g genotyping_method -d analysis_date -a analysis_dr\n\nArguments:\n  -m  MinION run name          (e.g. GAMP_22906B2)\n  -g  Genotyping method        (e.g. clair3_diploid)\n  -d  Analysis date            (format: YYYYMMDD, e.g. 20250511)\n  -a  Path to output directory (e.g. /mnt/c/.../output/)\n")
}

# Assign parsed values to variables
MinION_run_name <- arg_values["-m"]
genotyping_method <- arg_values["-g"]
analysis_date <- arg_values["-d"]
analysis_dr <- arg_values["-a"]

# Display parsed values
cat("✅ Parsed arguments:\n")
cat("MinION Run Name:", MinION_run_name, "\n")
cat("Genotyping Method:", genotyping_method, "\n")
cat("Analysis Date:", analysis_date, "\n")
cat("Analysis Directory:", paste0(analysis_dr,'/',MinION_run_name,'/',genotyping_method), "\n")

##### Load packages
library("tidyverse")
library("readxl")

# ---------------------------------------------------

### Output file names
overall_cov_fn <- paste0(analysis_dr, MinION_run_name, "/", genotyping_method, "/", "coverage/", MinION_run_name, "_", genotyping_method, "_coverage_amplicon_run_summary_", analysis_date, ".csv")
cov_by_sample_fn <- paste0(analysis_dr, MinION_run_name, "/", genotyping_method, "/", "coverage/", MinION_run_name, "_", genotyping_method, "_coverage_by_run_sample", analysis_date, ".csv")


### Import data
## Identify file paths
# find BED files
bedfile_dir <- paste0(analysis_dr,MinION_run_name,"/", genotyping_method, "/genome_coverage")

# amplicon position data
amp_pos_raw_fn <- "../data/analysis_dependencies/amplicon_gene_positions.csv"

# Nanopore sequence run metadata
meta_fn <- "../data/metadata/metadata_ont_multiplex_barcodes_test.xlsx"

## Read in the data
amp_pos_raw <- read.csv(amp_pos_raw_fn) 
meta <- read_xlsx(meta_fn, sheet = paste0("multiplex_", multiplex_run)) 

### Set coverage threshold
min_cov_threshold <- 50

### Manipulations
## Prepare amplicon position data
amp_pos_crt <- amp_pos_raw %>%
  filter(gene=="crt")

amp_pos_dhfr <- amp_pos_raw %>%
  filter(gene=="dhfr")

amp_pos_dhps <- amp_pos_raw %>%
  filter(gene=="dhps")

amp_pos_mdr1 <- amp_pos_raw %>%
  filter(gene=="mdr1")

amp_pos_kelch13 <- amp_pos_raw %>%
  filter(gene=="kelch13")

amp_pos_csp <- amp_pos_raw %>%
  filter(gene=="csp")

#amp_pos_msp1 <- amp_pos_raw %>%
#  filter(gene=="msp1")

## ONT metadata
samp_names <- meta %>%
  select(ont_barcode, sample_id, patient_id, ont_multiplex_group, ont_seq_date)


########## CRT
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_crt.bedGraph", full.names=TRUE)

if (length(filenames) == 0) {
  stop("No files found matching the pattern '_crt.bedGraph' in the directory: ", bedfile_dir)
}

ldf <- lapply(filenames, read.table)
max_count <- length(ldf)

### Filter variants to only include those within the amplicon target region
list_crt <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_crt$amplicon_start) %>%
    filter(end < amp_pos_crt$amplicon_end)
  
  list_crt[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_crt_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_crt[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  
  ## Reproduce barcode ID
  id <- ifelse(nchar(j)==1, paste0('barcode0', as.character(j)),
               ifelse(nchar(j)==2, paste0('barcode', as.character(j)), NA)
  )
  
  summdf$ont_barcode <- id
  summdf$gene_target <- "crt"
  
  summdf2 <- merge(summdf, samp_names, by="ont_barcode",
                   all.x=T,all.y=F)
  
  colns <- c("ont_barcode", "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean",
             "coverage_Q3", "coverage_max", "gene_target", "sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date")
  
  colnames(summdf2) <- colns
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_crt_cov[[j]] <- summdf2
  
}

cov_df_crt <- data.frame(do.call("rbind", list_crt_cov))


# summary stats
summary_cov_crt <- summary(cov_df_crt$coverage_median)

median_median_cov_crt <- median(cov_df_crt$coverage_median)
median_Q1_cov_crt <- round(as.numeric(summary(cov_df_crt$coverage_median)[2]),0)
median_Q3_cov_crt <- round(as.numeric(summary(cov_df_crt$coverage_median)[5]),0)
median_IQR_cov_crt <- median_Q3_cov_crt - median_Q1_cov_crt
lower_cov_5pcnt_crt <- round(as.numeric(quantile(cov_df_crt$coverage_median, probs=.05)),0)
lower_cov_10pcnt_crt <- round(as.numeric(quantile(cov_df_crt$coverage_median, probs=.1)),0)



########## dhfr
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_dhfr.bedGraph", full.names=TRUE)
ldf <- lapply(filenames, read.table)

max_count <- length(ldf)

### Filter variants to only include those within the amplicon target region
list_dhfr <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_dhfr$amplicon_start) %>%
    filter(end < amp_pos_dhfr$amplicon_end)
  
  list_dhfr[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_dhfr_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_dhfr[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  
  ## Reproduce barcode ID
  id <- ifelse(nchar(j)==1, paste0('barcode0', as.character(j)),
               ifelse(nchar(j)==2, paste0('barcode', as.character(j)), NA)
  )
  
  summdf$ont_barcode <- id
  summdf$gene_target <- "dhfr"
  
  summdf2 <- merge(summdf, samp_names, by="ont_barcode",
                   all.x=T,all.y=F)
  
  colns <- c("ont_barcode", "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean",
             "coverage_Q3", "coverage_max", "gene_target", "sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date")
  
  colnames(summdf2) <- colns
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_dhfr_cov[[j]] <- summdf2
  
}

cov_df_dhfr <- data.frame(do.call("rbind", list_dhfr_cov))

# summary stats
summary_cov_dhfr <- summary(cov_df_dhfr$coverage_median)

median_median_cov_dhfr <- median(cov_df_dhfr$coverage_median)
median_Q1_cov_dhfr <- round(as.numeric(summary(cov_df_dhfr$coverage_median)[2]),0)
median_Q3_cov_dhfr <- round(as.numeric(summary(cov_df_dhfr$coverage_median)[5]),0)
median_IQR_cov_dhfr <- median_Q3_cov_dhfr - median_Q1_cov_dhfr
lower_cov_5pcnt_dhfr <- round(as.numeric(quantile(cov_df_dhfr$coverage_median, probs=.05)),0)
lower_cov_10pcnt_dhfr <- round(as.numeric(quantile(cov_df_dhfr$coverage_median, probs=.1)),0)



########## dhps
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_dhps.bedGraph", full.names=TRUE)
ldf <- lapply(filenames, read.table)

max_count <- length(ldf)

### Filter variants to only include those within the amplicon target region
list_dhps <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_dhps$amplicon_start) %>%
    filter(end < amp_pos_dhps$amplicon_end)
  
  list_dhps[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_dhps_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_dhps[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  
  ## Reproduce barcode ID
  id <- ifelse(nchar(j)==1, paste0('barcode0', as.character(j)),
               ifelse(nchar(j)==2, paste0('barcode', as.character(j)), NA)
  )
  
  summdf$ont_barcode <- id
  summdf$gene_target <- "dhps"
  
  summdf2 <- merge(summdf, samp_names, by="ont_barcode",
                   all.x=T,all.y=F)
  
  colns <- c("ont_barcode", "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean",
             "coverage_Q3", "coverage_max", "gene_target", "sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date")
  
  colnames(summdf2) <- colns
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_dhps_cov[[j]] <- summdf2
  
}

cov_df_dhps <- data.frame(do.call("rbind", list_dhps_cov))

# summary stats
summary_cov_dhps <- summary(cov_df_dhps$coverage_median)

median_median_cov_dhps <- median(cov_df_dhps$coverage_median)
median_Q1_cov_dhps <- round(as.numeric(summary(cov_df_dhps$coverage_median)[2]),0)
median_Q3_cov_dhps <- round(as.numeric(summary(cov_df_dhps$coverage_median)[5]),0)
median_IQR_cov_dhps <- median_Q3_cov_dhps - median_Q1_cov_dhps
lower_cov_5pcnt_dhps <- round(as.numeric(quantile(cov_df_dhps$coverage_median, probs=.05)),0)
lower_cov_10pcnt_dhps <- round(as.numeric(quantile(cov_df_dhps$coverage_median, probs=.1)),0)



########## mdr1
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_mdr1.bedGraph", full.names=TRUE)
ldf <- lapply(filenames, read.table)

max_count <- length(ldf)

### Filter variants to only include those within the amplicon target region
list_mdr1 <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_mdr1$amplicon_start) %>%
    filter(end < amp_pos_mdr1$amplicon_end)
  
  list_mdr1[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_mdr1_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_mdr1[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  
  ## Reproduce barcode ID
  id <- ifelse(nchar(j)==1, paste0('barcode0', as.character(j)),
               ifelse(nchar(j)==2, paste0('barcode', as.character(j)), NA)
  )
  
  summdf$ont_barcode <- id
  summdf$gene_target <- "mdr1"
  
  summdf2 <- merge(summdf, samp_names, by="ont_barcode",
                   all.x=T,all.y=F)
  
  colns <- c("ont_barcode", "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean",
             "coverage_Q3", "coverage_max", "gene_target", "sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date")
  
  colnames(summdf2) <- colns
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_mdr1_cov[[j]] <- summdf2
  
}

cov_df_mdr1 <- data.frame(do.call("rbind", list_mdr1_cov))

# summary stats
summary_cov_mdr1 <- summary(cov_df_mdr1$coverage_median)

median_median_cov_mdr1 <- median(cov_df_mdr1$coverage_median)
median_Q1_cov_mdr1 <- round(as.numeric(summary(cov_df_mdr1$coverage_median)[2]),0)
median_Q3_cov_mdr1 <- round(as.numeric(summary(cov_df_mdr1$coverage_median)[5]),0)
median_IQR_cov_mdr1 <- median_Q3_cov_mdr1 - median_Q1_cov_mdr1
lower_cov_5pcnt_mdr1 <- round(as.numeric(quantile(cov_df_mdr1$coverage_median, probs=.05)),0)
lower_cov_10pcnt_mdr1 <- round(as.numeric(quantile(cov_df_mdr1$coverage_median, probs=.1)),0)



########## kelch13
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_k13.bedGraph", full.names=TRUE)
ldf <- lapply(filenames, read.table)

max_count <- length(ldf)

### Filter variants to only include those within the amplicon target region
list_k13 <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_kelch13$amplicon_start) %>%
    filter(end < amp_pos_kelch13$amplicon_end)
  
  list_k13[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_k13_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_k13[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  
  ## Reproduce barcode ID
  id <- ifelse(nchar(j)==1, paste0('barcode0', as.character(j)),
               ifelse(nchar(j)==2, paste0('barcode', as.character(j)), NA)
  )
  
  summdf$ont_barcode <- id
  summdf$gene_target <- "kelch13"
  
  summdf2 <- merge(summdf, samp_names, by="ont_barcode",
                   all.x=T,all.y=F)
  
  colns <- c("ont_barcode", "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean",
             "coverage_Q3", "coverage_max", "gene_target", "sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date")
  
  colnames(summdf2) <- colns
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_k13_cov[[j]] <- summdf2
  
}

cov_df_k13 <- data.frame(do.call("rbind", list_k13_cov))

# summary stats
summary_cov_k13 <- summary(cov_df_k13$coverage_median)

median_median_cov_k13 <- median(cov_df_k13$coverage_median)
median_Q1_cov_k13 <- round(as.numeric(summary(cov_df_k13$coverage_median)[2]),0)
median_Q3_cov_k13 <- round(as.numeric(summary(cov_df_k13$coverage_median)[5]),0)
median_IQR_cov_k13 <- median_Q3_cov_k13 - median_Q1_cov_k13
lower_cov_5pcnt_k13 <- round(as.numeric(quantile(cov_df_k13$coverage_median, probs=.05)),0)
lower_cov_10pcnt_k13 <- round(as.numeric(quantile(cov_df_k13$coverage_median, probs=.1)),0)



########## csp
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_csp.bedGraph", full.names=TRUE)
ldf <- lapply(filenames, read.table)

max_count <- length(ldf)

### Filter variants to only include those within the amplicon target region
list_csp <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_csp$amplicon_start) %>%
    filter(end < amp_pos_csp$amplicon_end)
  
  list_csp[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_csp_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_csp[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  
  ## Reproduce barcode ID
  id <- ifelse(nchar(j)==1, paste0('barcode0', as.character(j)),
               ifelse(nchar(j)==2, paste0('barcode', as.character(j)), NA)
  )
  
  summdf$ont_barcode <- id
  summdf$gene_target <- "csp"
  
  summdf2 <- merge(summdf, samp_names, by="ont_barcode",
                   all.x=T,all.y=F)
  
  colns <- c("ont_barcode", "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean",
             "coverage_Q3", "coverage_max", "gene_target", "sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date")
  
  colnames(summdf2) <- colns
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_csp_cov[[j]] <- summdf2
  
}

cov_df_csp <- data.frame(do.call("rbind", list_csp_cov))

# summary stats
summary_cov_csp <- summary(cov_df_csp$coverage_median)

median_median_cov_csp <- median(cov_df_csp$coverage_median)
median_Q1_cov_csp <- round(as.numeric(summary(cov_df_csp$coverage_median)[2]),0)
median_Q3_cov_csp <- round(as.numeric(summary(cov_df_csp$coverage_median)[5]),0)
median_IQR_cov_csp <- median_Q3_cov_k13 - median_Q1_cov_csp
lower_cov_5pcnt_csp <- round(as.numeric(quantile(cov_df_csp$coverage_median, probs=.05)),0)
lower_cov_10pcnt_csp <- round(as.numeric(quantile(cov_df_csp$coverage_median, probs=.1)),0)



########## msp1
# read in bed file
#filenames <- list.files(bedfile_dir, pattern="_msp1.bedGraph", full.names=TRUE)
#ldf <- lapply(filenames, read.table)
#
#max_count <- length(ldf)
#
### Filter variants to only include those within the amplicon target region
#list_msp1 <- list()
#for (i in 1:max_count) {
#  
#  ldf[[i]] <- ldf[[i]] %>%
#    rename(gene_id = V1,
#           start = V2,
#           end = V3,
#           depth = V4)
#  
#  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
#  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
#  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
#  
#  ldf[[i]] <- ldf[[i]] %>%
#    filter(start > amp_pos_msp1$amplicon_start) %>%
#    filter(end < amp_pos_msp1$amplicon_end)
#  
#  list_msp1[[i]] <- ldf[[i]]
#  
#}
#
### Now extract coverage statistics for each amplicon
#list_msp1_cov <- list()
#
#for (j in 1:max_count) {
#  
#  summ <- summary(list_msp1[[j]]$depth)
#  summdf <- data.frame(summ=matrix(summ, ncol=6))
#  colnames(summdf) <- names(summ)
#  
#  ## Reproduce barcode ID
#  id <- ifelse(nchar(j)==1, paste0('barcode0', as.character(j)),
#               ifelse(nchar(j)==2, paste0('barcode', as.character(j)), NA)
#  )
#  
#  summdf$ont_barcode <- id
#  summdf$gene_target <- "msp1"
#  
#  summdf2 <- merge(summdf, samp_names, by="ont_barcode",
#                   all.x=T,all.y=F)
#  
#  colns <- c("ont_barcode", "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean",
#             "coverage_Q3", "coverage_max", "gene_target", "sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date")
#  
#  colnames(summdf2) <- colns
#  
#  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
#               "gene_target",
#               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
#  
#  summdf2 <- summdf2 %>%
#    select(all_of(col_ord))
#  
#  list_msp1_cov[[j]] <- summdf2
#  
#}
#
#cov_df_msp1 <- data.frame(do.call("rbind", list_msp1_cov))
#
## summary stats
#summary_cov_msp1 <- summary(cov_df_msp1$coverage_median)
#
#median_median_cov_msp1 <- median(cov_df_msp1$coverage_median)
#median_Q1_cov_msp1 <- round(as.numeric(summary(cov_df_msp1$coverage_median)[2]),0)
#median_Q3_cov_msp1 <- round(as.numeric(summary(cov_df_msp1$coverage_median)[5]),0)
#median_IQR_cov_msp1 <- median_Q3_cov_msp1 - median_Q1_cov_msp1
#lower_cov_5pcnt_msp1 <- round(as.numeric(quantile(cov_df_msp1$coverage_median, probs=.05)),0)
#lower_cov_10pcnt_msp1 <- round(as.numeric(quantile(cov_df_msp1$coverage_median, probs=.1)),0)
#


### Summary statistics for the whole run
gene_target <- c("crt", "dhfr", "dhps", "mdr1", "kelch13", "csp") #msp1

coverage_median <- c(median_median_cov_crt, median_median_cov_dhfr, median_median_cov_dhps, median_median_cov_mdr1, median_median_cov_k13, median_median_cov_csp) #median_median_cov_msp1
coverage_Q1 <- c(median_Q1_cov_crt, median_Q1_cov_dhfr, median_Q1_cov_dhps, median_Q1_cov_mdr1, median_Q1_cov_k13, median_Q1_cov_csp) #median_Q1_cov_msp1
coverage_Q3 <- c(median_Q3_cov_crt, median_Q3_cov_dhfr, median_Q3_cov_dhps, median_Q3_cov_mdr1, median_Q3_cov_k13, median_Q3_cov_csp) #median_Q3_cov_msp1
coverage_IQR <- c(median_IQR_cov_crt, median_IQR_cov_dhfr, median_IQR_cov_dhps, median_IQR_cov_mdr1, median_IQR_cov_k13, median_IQR_cov_csp) #median_IQR_cov_msp1
coverage_lower_5pcnt <- c(lower_cov_5pcnt_crt, lower_cov_5pcnt_dhfr, lower_cov_5pcnt_dhps, lower_cov_5pcnt_mdr1, lower_cov_5pcnt_k13, lower_cov_5pcnt_csp) #lower_cov_5pcnt_msp1
coverage_lower_10pcnt <- c(lower_cov_10pcnt_crt, lower_cov_10pcnt_dhfr, lower_cov_10pcnt_dhps, lower_cov_10pcnt_mdr1, lower_cov_10pcnt_k13, lower_cov_10pcnt_csp) #lower_cov_10pcnt_msp1

median_cov_stats <- data.frame(gene_target, coverage_median, coverage_Q1, coverage_Q3, coverage_IQR, coverage_lower_5pcnt, coverage_lower_10pcnt)

median_cov_stats



#### Identify low coverage samples - RELATIVE lower 5% or 10%
cov_df_crt$cov_below_5pcnt <-
  ifelse(cov_df_crt$coverage_median < lower_cov_5pcnt_crt,
         TRUE, FALSE)
cov_df_crt$cov_below_10pcnt <-
  ifelse(cov_df_crt$coverage_median < lower_cov_10pcnt_crt,
         TRUE, FALSE)

cov_df_dhfr$cov_below_5pcnt <-
  ifelse(cov_df_dhfr$coverage_median < lower_cov_5pcnt_dhfr,
         TRUE, FALSE)
cov_df_dhfr$cov_below_10pcnt <-
  ifelse(cov_df_dhfr$coverage_median < lower_cov_10pcnt_dhfr,
         TRUE, FALSE)

cov_df_dhps$cov_below_5pcnt <-
  ifelse(cov_df_dhps$coverage_median < lower_cov_5pcnt_dhps,
         TRUE, FALSE)
cov_df_dhps$cov_below_10pcnt <-
  ifelse(cov_df_dhps$coverage_median < lower_cov_10pcnt_dhps,
         TRUE, FALSE)

cov_df_mdr1$cov_below_5pcnt <-
  ifelse(cov_df_mdr1$coverage_median < lower_cov_5pcnt_mdr1,
         TRUE, FALSE)
cov_df_mdr1$cov_below_10pcnt <-
  ifelse(cov_df_mdr1$coverage_median < lower_cov_10pcnt_mdr1,
         TRUE, FALSE)

cov_df_k13$cov_below_5pcnt <-
  ifelse(cov_df_k13$coverage_median < lower_cov_5pcnt_k13,
         TRUE, FALSE)
cov_df_k13$cov_below_10pcnt <-
  ifelse(cov_df_k13$coverage_median < lower_cov_10pcnt_k13,
         TRUE, FALSE)

cov_df_csp$cov_below_5pcnt <-
  ifelse(cov_df_csp$coverage_median < lower_cov_5pcnt_csp,
         TRUE, FALSE)
cov_df_csp$cov_below_10pcnt <-
  ifelse(cov_df_csp$coverage_median < lower_cov_10pcnt_csp,
         TRUE, FALSE)

#cov_df_msp1$cov_below_5pcnt <-
#  ifelse(cov_df_msp1$coverage_median < lower_cov_5pcnt_msp1,
#         TRUE, FALSE)
#cov_df_msp1$cov_below_10pcnt <-
#  ifelse(cov_df_msp1$coverage_median < lower_cov_10pcnt_msp1,
#         TRUE, FALSE)



#### Identify low coverage samples - ABSOLUTE coverage cutoff (defined above eg 50x)
cov_df_crt$coverage_above_threshold <-
  ifelse(cov_df_crt$coverage_median > min_cov_threshold,
         "TRUE", "FALSE")

cov_df_dhfr$coverage_above_threshold <-
  ifelse(cov_df_dhfr$coverage_median > min_cov_threshold,
         "TRUE", "FALSE")

cov_df_dhps$coverage_above_threshold <-
  ifelse(cov_df_dhps$coverage_median > min_cov_threshold,
         "TRUE", "FALSE")

cov_df_mdr1$coverage_above_threshold <-
  ifelse(cov_df_mdr1$coverage_median > min_cov_threshold,
         "TRUE", "FALSE")

cov_df_k13$coverage_above_threshold <-
  ifelse(cov_df_k13$coverage_median > min_cov_threshold,
         "TRUE", "FALSE")

cov_df_csp$coverage_above_threshold <-
  ifelse(cov_df_csp$coverage_median > min_cov_threshold,
         "TRUE", "FALSE")

#cov_df_msp1$coverage_above_threshold <-
#  ifelse(cov_df_msp1$coverage_median > min_cov_threshold,
#         "TRUE", "FALSE")



#### Merge all the data
cov_df_sub_crt <- cov_df_crt %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(crt_coverage_median = coverage_median,
         crt_coverage_below_5pcnt = cov_below_5pcnt,
         crt_coverage_below_10pcnt = cov_below_10pcnt,
         crt_coverage_above_threshold = coverage_above_threshold)

cov_df_sub_dhfr <- cov_df_dhfr %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(dhfr_coverage_median = coverage_median,
         dhfr_coverage_below_5pcnt = cov_below_5pcnt,
         dhfr_coverage_below_10pcnt = cov_below_10pcnt,
         dhfr_coverage_above_threshold = coverage_above_threshold)

cov_df_sub_dhps <- cov_df_dhps %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(dhps_coverage_median = coverage_median,
         dhps_coverage_below_5pcnt = cov_below_5pcnt,
         dhps_coverage_below_10pcnt = cov_below_10pcnt,
         dhps_coverage_above_threshold = coverage_above_threshold)

cov_df_sub_mdr1 <- cov_df_mdr1 %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(mdr1_coverage_median = coverage_median,
         mdr1_coverage_below_5pcnt = cov_below_5pcnt,
         mdr1_coverage_below_10pcnt = cov_below_10pcnt,
         mdr1_coverage_above_threshold = coverage_above_threshold)

cov_df_sub_k13 <- cov_df_k13 %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(k13_coverage_median = coverage_median,
         k13_coverage_below_5pcnt = cov_below_5pcnt,
         k13_coverage_below_10pcnt = cov_below_10pcnt,
         k13_coverage_above_threshold = coverage_above_threshold)

cov_df_sub_csp <- cov_df_csp %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(csp_coverage_median = coverage_median,
         csp_coverage_below_5pcnt = cov_below_5pcnt,
         csp_coverage_below_10pcnt = cov_below_10pcnt,
         csp_coverage_above_threshold = coverage_above_threshold)

#cov_df_sub_msp1 <- cov_df_msp1 %>%
#  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
#  rename(msp1_coverage_median = coverage_median,
#         msp1_coverage_below_5pcnt = cov_below_5pcnt,
#         msp1_coverage_below_10pcnt = cov_below_10pcnt,
#         msp1_coverage_above_threshold = coverage_above_threshold)

# merge
cov_merged_df1 <- merge(cov_df_sub_crt,cov_df_sub_dhfr,
                        by="sample_id",
                        all.x=T, all.y=T)

cov_merged_df2 <- merge(cov_merged_df1,cov_df_sub_dhps,
                        by="sample_id",
                        all.x=T, all.y=T)

cov_merged_df3 <- merge(cov_merged_df2,cov_df_sub_mdr1,
                        by="sample_id",
                        all.x=T, all.y=T)

cov_merged_df4 <- merge(cov_merged_df3,cov_df_sub_k13,
                        by="sample_id",
                        all.x=T, all.y=T)

cov_merged_df_final <- merge(cov_merged_df4,cov_df_sub_csp,
                             by="sample_id",
                             all.x=T, all.y=T)

#cov_merged_df_final <- merge(cov_merged_df5,cov_df_sub_msp1,
#                        by="sample_id",
#                        all.x=T, all.y=T)

head(cov_merged_df_final)
median_cov_stats


##### Save outputs
# Ensure the output directory exists
dir.create(dirname(overall_cov_fn), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cov_by_sample_fn), recursive = TRUE, showWarnings = FALSE)

write.csv(median_cov_stats, file=overall_cov_fn, row.names=FALSE)
write.csv(cov_merged_df_final, file=cov_by_sample_fn, row.names=FALSE)


