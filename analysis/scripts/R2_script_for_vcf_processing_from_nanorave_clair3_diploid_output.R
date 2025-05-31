########### R script for processing VCFs from Nanopore run GAMP_22906B2 (example)
### Clair3 diploid genotyper option with nano-rave pipeline

##### Set wd
path_setwd <- getwd()
setwd(path_setwd)


##### assign sample and amplicon names + genotyping method

multiplex_run = "B2"
#multiplex_seq_date = as.Date("2022-09-06")

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define argument names and default values
arg_names <- c("-m", "-g", "-a")
arg_values <- setNames(rep(NA, length(arg_names)), arg_names)

# Ensure sufficient arguments are provided
if (length(args) < length(arg_names) * 2) {
  stop("\n❌ Missing required arguments.\n\nUsage:\n  Rscript script.R -m MinION_run_name -g genotyping_method -a analysis_dr\n\nArguments:\n  -m  MinION run name          (e.g. GAMP_22906B2)\n  -g  Genotyping method        (e.g. clair3_diploid)\n  -a  Path to output directory (e.g. /mnt/c/.../output/)\n")
}

# Parse arguments using vectorized lookup
for (i in seq(1, length(args), by = 2)) {
  if (args[i] %in% arg_names) {
    arg_values[args[i]] <- args[i + 1]
  }
}

# Check if any required arguments are missing
if (any(is.na(arg_values))) {
  stop("\n❌ Missing required arguments.\n\nUsage:\n  Rscript script.R -m MinION_run_name -g genotyping_method -a analysis_dr\n\nArguments:\n  -m  MinION run name          (e.g. GAMP_22906B2)\n  -g  Genotyping method        (e.g. clair3_diploid)\n  -a  Path to output directory (e.g. /mnt/c/.../output/)\n")
}

# Assign parsed values to variables
MinION_run_name <- arg_values["-m"]
genotyping_method <- arg_values["-g"]
analysis_dr <- arg_values["-a"]

# Display parsed values
cat("✅ Parsed arguments:\n")
cat("MinION Run Name:", MinION_run_name, "\n")
cat("Genotyping Method:", genotyping_method, "\n")
cat("Analysis Directory:", paste0(analysis_dr,'/',MinION_run_name,'/',genotyping_method), "\n")


##### Load packages
library("tidyverse")
library("vcfR")
library("readxl")


out_variant_info_fn <- paste0(path_setwd,"/",analysis_dr, MinION_run_name, "/", genotyping_method, "/", "genotypes/", MinION_run_name, "_variants_linked_to_metadata_intermediate_file", ".csv") # for vcf_varcount_df4
out_variant_info2_fn <- paste0(path_setwd,"/",analysis_dr, MinION_run_name, "/", genotyping_method, "/", "genotypes/", MinION_run_name,"_variants_linked_to_metadata_intermediate_file2", ".csv") # for vcf_tbl_all_varinfo2
out_keysnp_allsamp_genotypes_fn <- paste0(path_setwd,"/",analysis_dr, MinION_run_name, "/", genotyping_method, "/", "genotypes/", MinION_run_name, "_genotype_calls_samplecols", ".csv") # for keysnps_allsamps_gt
out_keysnp_allsamp_genotypes_samprows_fn <- paste0(path_setwd,"/",analysis_dr, MinION_run_name, "/", genotyping_method, "/", "genotypes/", MinION_run_name, "_genotype_calls_haplotypes_samplerows_v2", ".csv") # for snp_calls_t3
out_k13_calls_info_fn <- paste0(path_setwd,"/",analysis_dr, MinION_run_name, "/", genotyping_method, "/", "genotypes/", MinION_run_name, "_k13_variant_info", ".csv") # for k13_snps_allsamps_df
out_csp_genotypes_long_fn <- paste0(path_setwd,"/",analysis_dr, MinION_run_name, "/", genotyping_method, "/", "genotypes/", MinION_run_name, "_csp_nref_genotype_calls_long", ".csv") # for vcf_tbl_all_csp
out_csp_genotypes_samps_fn <- paste0(path_setwd,"/",analysis_dr, MinION_run_name, "/", genotyping_method, "/", "genotypes/", MinION_run_name, "_csp_genotype_calls_samplecols", ".csv") # for csp_samp_calls
out_msp1_genotypes_long_fn <- paste0(path_setwd,"/",analysis_dr, MinION_run_name, "/", genotyping_method, "/", "genotypes/", MinION_run_name, "_msp1_nref_genotype_calls_long", ".csv") # for vcf_tbl_all_msp1
out_plot_msp1_genetic_distance_fn <- paste0(path_setwd,"/",analysis_dr, MinION_run_name, "/", genotyping_method, "/", "genotypes/", MinION_run_name, "_plot_msp1_genetic_distance_NJT", ".png") # for msp1 NJT



##### Import data

## list of variants to call
snp_data <- read_excel("../data/analysis_dependencies/DR_variant_info_v2.xlsx", sheet = 1)

## Multiplex sample metadata
meta <- read_excel("../data/metadata/metadata_ont_multiplex_barcodes_test.xlsx", sheet =paste0("multiplex_",multiplex_run))

## kelch13 info
k13_annotated <- read.csv("../data/analysis_dependencies/k13_seq_annotated.csv")
genetic_code <- read.csv("../data/analysis_dependencies/genetic_code.csv")

### Import the vcf files
# construct directory for importing vcf files
nflow_results_dir <- paste0(analysis_dr,"/", MinION_run_name, "/", genotyping_method, "/variant_calling_unzip/")


## Extract the file names
setwd(nflow_results_dir)
dataFiles <- list() # creates a list
listvcf <- dir(pattern = "*.vcf") # creates the list of all the vcf files in the directory
for (k in 1:length(listvcf)){
  dataFiles[[k]] <- read.vcfR(listvcf[k])
}

# Alternative method for getting all the files
#dataFiles <- lapply(Sys.glob(paste(nflow_results_dir,MinION_run_name,"_barcode*",".vcf",sep="")), read.vcfR)

## Check it out
length(dataFiles)
dataFiles[[12]]
dataFiles[[2]]



### Check whether any vcf files don't have any variants

n <- length(dataFiles)

for (i in 1:n) {
  
  kn <- nrow(extract.gt(dataFiles[[i]]))
  
  result <- ifelse(kn>0,
                   print(paste("PASS! Vcf file no. ", i, " contains >=1 variant", sep="")),
                   print(paste("** FAIL! Vcf file no. ", i, " contains zero variants! **", sep="")))
}



######## Create dataframe of variant counts

mylist <- list()

n <- length(dataFiles)

for (i in 1:n) {
  
  # Count number of variants in the vcf file
  kn <- nrow(extract.gt(dataFiles[[i]]))

  mylist[[i]] <- kn
  
}

vcf_varcount_df <- data.frame(do.call("rbind", mylist)) %>%
  rename("variant_count" = "do.call..rbind...mylist.")
vcf_varcount_df$vcf_fn <- listvcf
vcf_varcount_df

vcf_varcount_df$zero_variants <- ifelse(vcf_varcount_df$variant_count<1, TRUE, FALSE)



## amplicon sequence names
vcf_varcount_df2 <- within(vcf_varcount_df, {
  
  # Convert to character
  vcf_fn <- as.character(vcf_fn)
  
  # Convert to integer for the numerical value of unique ID
  amplicon_id <- as.character(as.integer(factor(vcf_fn, levels = unique(vcf_fn))))
  
  # Rename with study format
  amplicon_id <- ifelse(nchar(amplicon_id) == 1, paste0('GH22', multiplex_run, '_SEQ', '00', amplicon_id),
                        ifelse(nchar(amplicon_id) == 2, paste0('GH22', multiplex_run, '_SEQ', '0', amplicon_id),
                               ifelse(nchar(amplicon_id) == 3, paste0('GH22', multiplex_run, '_SEQ', amplicon_id), NA)
                        )
  )
})


### Barcodes
## Extract the barcodes
vcf_varcount_df2$ont_barcode <- vcf_varcount_df2$vcf_fn
vcf_varcount_df2$ont_barcode <- substr(vcf_varcount_df2$ont_barcode, 14, 22) #### Change lengths if ONT run name is different length

#vcf_varcount_df2$ont_multiplex_run <- multiplex_run ### This will get imported with the metadata
#vcf_varcount_df2$ont_seq_date <- multiplex_seq_date

## Extract amplicon gene names
vcf_varcount_df2$amplicon_gene_target <- vcf_varcount_df2$vcf_fn
vcf_varcount_df2$amplicon_gene_target <- str_sub(vcf_varcount_df2$amplicon_gene_target, 24, -5)
vcf_varcount_df2$amplicon_gene_target

head(vcf_varcount_df2)


## Add in the sample IDs from metadata file
colnames(meta)
meta_ed <- meta %>% select(-any_of(c("notes", "patient_id")))


vcf_varcount_df3 <- merge(vcf_varcount_df2,meta_ed,
                          by="ont_barcode",
                          all.x=T, all.y=F)

col_ord1 <- c("sample_id",
              "collection_site",
              "ont_seq_date",
              "ont_multiplex_group",
              "ont_barcode",
              "amplicon_id",
              "amplicon_gene_target",
              "vcf_fn",
              "variant_count",
              "zero_variants")

vcf_varcount_df3 <- vcf_varcount_df3 %>% select(all_of(col_ord1))

vcf_varcount_df3 <- vcf_varcount_df3 %>%
  mutate(sample_id = ifelse(is.na(sample_id), paste0("Sample_", gsub("\\.v$", "", ont_barcode)), sample_id))

vcf_varcount_df3 <- vcf_varcount_df3 %>%
  mutate(sample_id = gsub("\\.v.*$", "", sample_id))

head(vcf_varcount_df3,7)
nrow(vcf_varcount_df3)



######## Remove vcf files from list if they have zero variants

empty_vcfs <- which(vcf_varcount_df3$zero_variants == TRUE)

dataFiles_noempty <- dataFiles[-c(empty_vcfs)];

length(dataFiles)
length(dataFiles_noempty)

# keep track of sample order
vcf_varcount_df4 <- vcf_varcount_df3 %>%
  filter(zero_variants==FALSE) %>%
  arrange(vcf_fn)

head(vcf_varcount_df4,7)
nrow(vcf_varcount_df4)



###############

### Testing with 1 sample
dataFiles_noempty[3]
# Extract genotypes + variant info
temp <- vcfR2tidy(dataFiles_noempty[[length(dataFiles_noempty)]], format_fields = c("GT", "DP", "AF"))
fix <- temp$fix %>% select(-ChromKey)
gt <- temp$gt %>% select(-ChromKey)
meta <- temp$meta
# Combine into 1 table
fixgt <- left_join(fix, gt, by = "POS")



############## loooooooooooop

n2 <- length(dataFiles_noempty)

mylist2 <- list()

for (j in 1:n2) {
  
  ##### Turning vcf into formatted tibble for each gene amplicon
  # Extract genotypes + variant info
  temp <- vcfR2tidy(dataFiles_noempty[[j]], format_fields = c("GT", "DP", "AF"))
  fix <- temp$fix %>% select(-ChromKey)
  gt <- temp$gt %>% select(-ChromKey)
  meta <- temp$meta
  # Combine into 1 table
  fixgt <- left_join(fix, gt, by = "POS")
  # Select columns with anticipated new names
  gtcols <- c("Gene", "Pos", "Ref", "Alt", "DP", "Qual", "Filter", "Pileup_call", "Alignment_call", "GT", "AF")
  # Fiddle with formatting
  ampx_tbl_fix_sampx <- fixgt %>%
    rename(Gene=CHROM,
           Pos=POS,
           Ref=REF,
           Alt=ALT,
           Qual=QUAL,
           Filter=FILTER,
           Pileup_call=P,
           Alignment_call="F",
           GT=gt_GT,
           DP=gt_DP,
           AF=gt_AF) %>%  # rename things
    select(all_of(gtcols))  # select desired vcf components

  ampx_tbl_fix_sampx$seq_count_temp = j
  
  mylist2[[j]] <- ampx_tbl_fix_sampx
  
}

# This may have emptied some tibbles of any content, if there were no variants that make it through filtering
xtemp <- length(dataFiles_noempty)
ytemp <- length(mylist2)

ifelse(xtemp==ytemp,
       "No change in number of vcfs included in the list - all good :-D",
       "Number of vcfs has changed! Danger!")

vcf_tbl_all <- data.frame(do.call("rbind", mylist2)) %>%
  arrange(seq_count_temp) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_DR1_crt", "crt")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_AG1_csp", "csp")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_DR1_dhfr", "dhfr")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_DR1_dhps", "dhps")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_DR1_k13", "k13")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_DR1_mdr1", "mdr1")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_AG1_msp1", "msp1"))

# Add SNP names
vcf_tbl_all2 <- within(vcf_tbl_all, {
  snp_id = paste0(Gene,"_",Pos,"_",Ref,Alt)
})
  


###### Merge with sample metadata

### key step - need to ensure order matches ###

xtemp = n_distinct(vcf_tbl_all2$seq_count_temp)
ytemp = nrow(vcf_varcount_df4)

ifelse(xtemp==ytemp,
       "Unique seq counts matches variant count - all good :-D",
       "Number of variants has changed! Danger!")

vcf_varcount_df4$seq_count_temp <- row.names(vcf_varcount_df4)
vcf_varcount_df4$seq_count_temp <- as.numeric(vcf_varcount_df4$seq_count_temp)

vcf_tbl_all2$seq_count_temp <- as.numeric(vcf_tbl_all2$seq_count_temp)

vcf_tbl_all_varinfo <- merge(vcf_tbl_all2, vcf_varcount_df4,
                             by="seq_count_temp",
                             all.x=T, all.y=T) %>%
  arrange(seq_count_temp)

head(vcf_tbl_all_varinfo,15)

nrow(vcf_tbl_all2)
nrow(vcf_tbl_all_varinfo)

xtemp = nrow(vcf_tbl_all2)
ytemp = nrow(vcf_tbl_all_varinfo)

ifelse(xtemp==ytemp,
       "Number of variants is the same as before merging with metadata - all good :-D",
       "Number of variants has changed! Danger!")


# Check the variant gene matches the vcf call gene names - should be identical
vcf_tbl_all_varinfo$gene_check_temp <- ifelse(vcf_tbl_all_varinfo$Gene==vcf_tbl_all_varinfo$amplicon_gene_target,
       0, 1)
ifelse(sum(vcf_tbl_all_varinfo$gene_check_temp)==0,
       "Gene names from vcf calls and variant metadata all match up - all good :-D",
       "Gene names don't match! Danger!")

### Now filter SNPs by quality score, only include biallelic pass, and DP>50

vcf_tbl_all_varinfo2 <- vcf_tbl_all_varinfo %>%
  filter(Filter=="PASS" & # quality pass
           Qual >= 1 & # Quality score >1
           nchar(Ref) < 2 & # single SNP ref
           nchar(Alt) < 2 & # single SNP alt (ie biallelic)
           DP >= 50) # >50x coverage depth

dim(vcf_tbl_all_varinfo)
dim(vcf_tbl_all_varinfo2)


## Clean up columns
col_ord2 <- c("sample_id",
              "collection_site",
              "ont_seq_date",
              "ont_multiplex_group",
              "ont_barcode",
              "amplicon_id",
              "vcf_fn",
              "Gene",
              "Pos",
              "Ref",
              "Alt",
              "Qual",
              "Filter",
              "Pileup_call",
              "Alignment_call",
              "DP",
              "AF",
              "GT",
              "snp_id")
              
vcf_tbl_all_varinfo2 <- vcf_tbl_all_varinfo2 %>% select(all_of(col_ord2))

head(vcf_tbl_all_varinfo2,15)
dim(vcf_tbl_all_varinfo2)


### Select the majority allele & turn genotypes into haplotypes
AF_majority_threshold = 0.51 # Set the threshold to define a majority call for hets, eg 51% nref-AF

vcf_tbl_all_varinfo2$GT_majority <- ifelse(vcf_tbl_all_varinfo2$AF > AF_majority_threshold,
                                           1,
                                           0)
vcf_tbl_all_varinfo2$GT
vcf_tbl_all_varinfo2$GT_majority

## Add a column whether likely het, defined by a cutoff
het_cutoff = 0.8 # e.g. if using 80%, then non-ref AF >80%==hom-nref or <20%==hom-ref; in between == het
vcf_tbl_all_varinfo2$is_het <- ifelse( (vcf_tbl_all_varinfo2$AF < het_cutoff) &
                                         (vcf_tbl_all_varinfo2$AF > (1-het_cutoff) ),
                                       TRUE, FALSE)

summary(vcf_tbl_all_varinfo2$is_het) # Tells you how many hets

vcf_tbl_all_varinfo2 <- vcf_tbl_all_varinfo2 %>% # Rename GT_majority to genotype
  rename(genotype=GT_majority)

# Note - The AF (allelic frequency) is the percentage of the alternative allele.
# 1-AF is the reference allele. If 1-AF >0.5, the reference allele is the majority.



### Add in SNP data for the samples - note this is ignoring kelch13 which we deal with separately (below)
colnames(snp_data)

snp_data_cols <- c("snp_id",
                   "chrom",
                   "chrom_pos",
                   "pos_in_gene_coding",
                   "pos_in_codon",
                   "snp_ref",
                   "snp_alt",
                   "codon_pos_in_gene",
                   "codon_ref",
                   "codon_alt",
                   "aa_ref",
                   "aa_alt",
                   "aa_mut_name",
                   "key_snp")

snp_data_v2 <- snp_data %>%
  select(all_of(snp_data_cols)) %>%
  filter(key_snp == TRUE) %>%
  select(-key_snp)

colnames(vcf_tbl_all_varinfo2)
sample_snp_calls_col <- c("sample_id",
                          "snp_id",
                          "genotype")
sample_snp_calls <- vcf_tbl_all_varinfo2 %>%
  select(all_of(sample_snp_calls_col))

### Make separate genotype calls for each sample in the batch

# Create integer for each group
sample_snp_calls <- sample_snp_calls %>%
  mutate(sample_integer = cumsum(!duplicated(sample_id)))


# loop through each sample integer to make list of genotype calls to add each one as a column to the key SNPs

maxn <- max(sample_snp_calls$sample_integer)

mylist_sample_vars <- list()

for (k in 1:maxn) {
  
  sub <- sample_snp_calls %>%
    filter(sample_integer==k) %>%
    select(c("snp_id", "sample_id", "genotype"))
  
  ncolname <- sub[1,"sample_id"]
  
  colnames(sub) <- c("snp_id","sample_id",ncolname)
  
  sub <- sub %>% select(all_of(c("snp_id", ncolname)))
  
  mylist_sample_vars[[k]] <- sub
  
}

#mylist_sample_vars[[3]]
length(mylist_sample_vars)

# get the list of SNPs to be included
snp_list <- snp_data_v2 %>%
  select("snp_id") %>%
  arrange("SNP_id")
head(snp_list)


# Add the genotype calls to the snp_list for each sample, to standardise row number
  
for(m in 1:maxn) {
  mylist_sample_vars[[m]] <- left_join(snp_list, mylist_sample_vars[[m]], by="snp_id") %>%
    arrange(snp_id) # keep snp_id for joining later
}

# Merge all sample columns by snp_id
keysnps_allsamps_gt <- reduce(mylist_sample_vars, full_join, by = "snp_id") %>%
  arrange(snp_id)

## Replace NA with wild-type genotype call
keysnps_allsamps_gt <- keysnps_allsamps_gt %>% replace(is.na(.), 0)


############ Check out the drug resistance (minus kelch13) genotype data!
keysnps_allsamps_gt





########### Deal with kelch13, csp and msp1 if we included msp1
#thing <- rep(1:1000, each=3)
#write.csv(thing, file="codon_num_blank.csv", row.names=F)

#### Get the k13 variants

colnames(vcf_tbl_all_varinfo2)
col_ord3 <- c("sample_id",
              "collection_site",
              "amplicon_id",
              "snp_id",
              "Gene",
              "Pos",
              "Ref",
              "Alt",
              "Qual",
              "DP",
              "genotype")

# k13
vcf_tbl_all_k13 <- vcf_tbl_all_varinfo2 %>%
  filter(Gene=="k13") %>%
  select(all_of(col_ord3)) %>%
  arrange(sample_id, snp_id)

# csp
vcf_tbl_all_csp <- vcf_tbl_all_varinfo2 %>%
  filter(Gene=="csp") %>%
  select(all_of(col_ord3)) %>%
  arrange(sample_id, snp_id)

# msp1
vcf_tbl_all_msp1 <- vcf_tbl_all_varinfo2 %>%
  filter(Gene=="msp1") %>%
  select(all_of(col_ord3)) %>%
  arrange(sample_id, snp_id)


##### Remove variants that are before or after my actual primers - unlikely to be real!!

vcf_tbl_all_k13$Pos <- as.numeric(vcf_tbl_all_k13$Pos)
vcf_tbl_all_csp$Pos <- as.numeric(vcf_tbl_all_csp$Pos)
vcf_tbl_all_msp1$Pos <- as.numeric(vcf_tbl_all_msp1$Pos)

vcf_tbl_all_k13 <- vcf_tbl_all_k13 %>%
  filter((Pos >= 1277) & (Pos <= 2145))
if (nrow(vcf_tbl_all_k13) == 0) {
  vcf_tbl_all_k13 <- tibble(
    sample_id = character(),
    collection_site = character(),
    amplicon_id = character(),
    snp_id = character(),
    Gene = character(),
    Pos = integer(),
    Ref = character(),
    Alt = character(),
    Qual = numeric(),
    DP = numeric(),
    genotype = character()
  )
}

vcf_tbl_all_csp <- vcf_tbl_all_csp %>%
  filter((Pos >= 168) & (Pos <= 1143))  # Actual end in 3D7 is 1143

vcf_tbl_all_msp1 <- vcf_tbl_all_msp1 %>%
  filter((Pos >= 105) & (Pos <= 5099)) # Actual end in 3D7 is 5099




########### kelch13 - special analysis
### Is anything left?? Otherwise, will get errors:
ifelse(nrow(vcf_tbl_all_k13)==0,
       "No kelch13 mutations detected, no need to continue analysis",
       "At least 1 kelch13 mutation - proceed with analysis")

######### Check whether there are multiple kelch mutations affecting the same codon
# script I've written deals with each variant individually, would struggle with 2 SNPs in same codon

## Check whether there are 2 SNPs in same codon across whole dataset
k13_pos_subs <- vcf_tbl_all_k13 %>%
  select(snp_id, Pos) %>%
  rename(pos = Pos)
k13_pos_subs

# Add to the k13 seq data
k13_seq_var <- left_join(k13_annotated, k13_pos_subs, "pos") %>%
  arrange(snp_id) %>%
  drop_na %>%
  distinct(snp_id, .keep_all=TRUE)

temptest <- n_distinct(k13_seq_var$codon_num)

ifelse(temptest==nrow(k13_seq_var),
       "All of the SNPs in this dataset occur in unique codons, so it's fine",
       "Danger! Multiple SNPs within same codon")



###### If there are any kelch13 variants then proceed with this bit
if(nrow(vcf_tbl_all_k13)>0) {
  
  # Create integer for each sample
  vcf_tbl_all_k13 <- vcf_tbl_all_k13 %>%
    mutate(sample_integer = cumsum(!duplicated(sample_id)))
  
  
  k13nmax <- max(vcf_tbl_all_k13$sample_integer)
  
  mylist_k13_samps <- list()
  
  for (q in 1:k13nmax) {
    
    # subset for each sample
    vcf_tbl_k13_q <- vcf_tbl_all_k13 %>% filter(sample_integer==q)
    
    # Create integer for each SNP
    vcf_tbl_k13_q <- vcf_tbl_k13_q %>%
      mutate(snp_integer = cumsum(!duplicated(snp_id)))
    
    # start nested loop for each SNP
    k13nmaxsub <- max(vcf_tbl_k13_q$snp_integer)
    
    mylist_k13_snps <- list()
    
    for (r in 1:k13nmaxsub) {
      
      vcf_tbl_k13_q$Pos <- as.numeric(vcf_tbl_k13_q$Pos)
      vcf_tbl_k13_q$Ref <- as.character(vcf_tbl_k13_q$Ref)
      vcf_tbl_k13_q$Alt <- as.character(vcf_tbl_k13_q$Alt)
      
      ## Define position and nucleotide of each SNP in the loop
      query_pos = vcf_tbl_k13_q[r,'Pos']
      query_alt_nt = vcf_tbl_k13_q[r,'Alt']
      
      ## Extract original codon and mutant codon
      nt_ref <- k13_annotated[query_pos,3]
      nt_codon_pos <- k13_annotated[query_pos,4]
      
      ## Define the ref codon
      ref_codon <- ifelse(nt_codon_pos==1, paste(nt_ref,k13_annotated[query_pos+1,3],k13_annotated[query_pos+2,3],sep=""),
                          ifelse(nt_codon_pos==2, paste(k13_annotated[query_pos-1,3],nt_ref,k13_annotated[query_pos+1,3],sep=""),
                                 ifelse(nt_codon_pos==3, paste(k13_annotated[query_pos-2,3],k13_annotated[query_pos-1,3],nt_ref,sep=""),
                                        paste("Something is wrong")
                                 )
                          )
      )
      ## Define the alt codon
      alt_codon <- ifelse(nt_codon_pos==1, paste(query_alt_nt,k13_annotated[query_pos+1,3],k13_annotated[query_pos+2,3],sep=""),
                          ifelse(nt_codon_pos==2, paste(k13_annotated[query_pos-1,3],query_alt_nt,k13_annotated[query_pos+1,3],sep=""),
                                 ifelse(nt_codon_pos==3, paste(k13_annotated[query_pos-2,3],k13_annotated[query_pos-1,3],query_alt_nt,sep=""),
                                        paste("Something is wrong")
                                 )
                          )
      )
      
      ## Translate codons into amino acids
      ref_codon_translation_info <- genetic_code %>%
        filter(codon==ref_codon) %>%
        rename(ref_codon = codon, ref_aa_letter = aa_letter, ref_aa_code = aa_code, ref_aa_name = aa_name)
      
      alt_codon_translation_info <- genetic_code %>%
        filter(codon==alt_codon) %>%
        rename(alt_codon = codon, alt_aa_letter = aa_letter, alt_aa_code = aa_code, alt_aa_name = aa_name)
      
      ## combine codon info
      aa_change_info <- cbind(ref_codon_translation_info,alt_codon_translation_info)
      
      ## Define mutation type
      aa_change_info$mut_type <-
        ifelse(aa_change_info$ref_aa_letter==aa_change_info$alt_aa_letter, "Synonymous",
               ifelse( ((aa_change_info$ref_aa_letter!=aa_change_info$alt_aa_letter) &
                          (aa_change_info$alt_aa_letter!="Stop")), "Missense",
                       ifelse( ((aa_change_info$ref_aa_letter!=aa_change_info$alt_aa_letter) &
                                  (aa_change_info$alt_aa_letter=="Stop")), "Nonsense",
                               "Something has gone wrong")
               )
        )
      
      ## Define the codon number 
      aa_change_info$codon_count <- ifelse(query_pos %% 3 == 0, query_pos/3, round((query_pos/3)+0.5,0) )
      
      ## Populate other variant fields
      aa_change_info$gene <- "k13"
      aa_change_info$pos <- query_pos
      aa_change_info$ref_nt <- nt_ref
      aa_change_info$alt_nt <- query_alt_nt
      aa_change_info$pos_in_codon <- nt_codon_pos
      
      aa_change_info$snp_id <- paste(aa_change_info$gene, "_", aa_change_info$pos, "_", aa_change_info$ref_nt, aa_change_info$alt_nt, sep="")
      
      aa_change_info$aa_mut <- paste(aa_change_info$ref_aa_letter, aa_change_info$codon_count, aa_change_info$alt_aa_letter, sep="")
      aa_change_info$aa_mut_id <- paste(aa_change_info$gene, "_", aa_change_info$ref_aa_letter, aa_change_info$codon_count, aa_change_info$alt_aa_letter, sep="")
      
      aa_change_info$codon_count <- as.numeric(aa_change_info$codon_count)
      
      ## Is the mutation within the propeller domain/ WHO region to worry about?
      aa_change_info$WHO_ARTR_region <- ifelse(aa_change_info$codon_count >= 349, TRUE, FALSE)
      
      ## Has the mutation previously associated with ART-R? (list from MalariaGEN GenRe technical documents, accessed 2022-08-07)
      k13_VOC <- c("P441L", "F446I", "G449A", "D452E", "N458Y", "C469Y", "C469F", "M476I", "K479I", "A481V", 
                   "Y493H", "R515K", "S522C", "P527L", "N537I", "N537D", "G538V", "R539T", "I543T", "P553L", 
                   "R561H", "V568G", "P574L", "R575K", "M579I", "C580Y", "D584V", "P667T", "F673I", "A675V", "H719N")
      
      aa_change_info$mut_ARTR_phenotype <- aa_change_info$aa_mut %in% k13_VOC
      
      ## Keep track of the sample
      aa_change_info$sample_id <- vcf_tbl_k13_q[1,"sample_id"]
      
      ## Keep track of SNP quality and DP
      aa_change_info$Qual <- vcf_tbl_k13_q[r,"Qual"]
      aa_change_info$DP <- vcf_tbl_k13_q[r,"DP"]
      
      ## Organise the columns
      cols_order <- c("sample_id",
                      "snp_id",
                      "Qual",
                      "DP",
                      "gene",
                      "pos",
                      "ref_nt",
                      "alt_nt",
                      "pos_in_codon",
                      "ref_codon",
                      "alt_codon",
                      "codon_count",
                      "WHO_ARTR_region",
                      "mut_type",
                      "ref_aa_letter",
                      "alt_aa_letter",
                      "ref_aa_code",
                      "alt_aa_code",
                      "ref_aa_name",
                      "alt_aa_name",
                      "aa_mut",
                      "aa_mut_id",
                      "mut_ARTR_phenotype")
      
      mutation_info <- aa_change_info %>% select(all_of(cols_order))
      
      mylist_k13_snps[[r]] <- mutation_info # Populate the list
      
    }
    
    k13_snps_df <- data.frame(do.call("rbind", mylist_k13_snps)) # Bind the list into a dataframe of SNPs for each sample
    
    mylist_k13_samps[[q]] <- k13_snps_df
    
  }
  
  k13_snps_allsamps_df <- data.frame(do.call("rbind", mylist_k13_samps)) # Bind the list into a df of samples
  
  k13_snps_allsamps_df
  
}




########################## High level interpretations: data manipulations

####### Check out the data to be used
keysnps_allsamps_gt

####### transposed version
snp_calls <- keysnps_allsamps_gt[2:ncol(keysnps_allsamps_gt)]
var_cols <- c(keysnps_allsamps_gt$snp_id)
snp_calls_t <- snp_calls %>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) 
names(snp_calls_t) <- c('sample_id', var_cols)


####### Create composite haplotype calls

## dhfr
# Options = IRNI (51I*-59R*-108N*-164I); or NRNI (51N-59R*-108N*-164I) or ICNI (51I*-59C-108N*-164I), or Other
# WT = NCSI
# Ensure column names are unique to avoid within() assignment errors
colnames(snp_calls_t) <- make.unique(colnames(snp_calls_t), sep = "_")

snp_calls_t2 <- within(snp_calls_t, {
  dhfr_haplotype <-
    ifelse( (dhfr_152_AT==1 &
               dhfr_175_TC==1 & 
               dhfr_323_GA==1 &
               dhfr_490_AT==0 ), "IRNI", 
            ifelse((dhfr_152_AT==0 &
                      dhfr_175_TC==1 & 
                      dhfr_323_GA==1 &
                      dhfr_490_AT==0 ), "NRNI",
                   ifelse((dhfr_152_AT==1 &
                             dhfr_175_TC==0 & 
                             dhfr_323_GA==1 &
                             dhfr_490_AT==0 ), "ICNI",
                          ifelse((dhfr_152_AT==0 &
                                    dhfr_175_TC==0 & 
                                    dhfr_323_GA==1 &
                                    dhfr_490_AT==0 ), "NCNI",
                                 ifelse((dhfr_152_AT==1 &
                                           dhfr_175_TC==1 & 
                                           dhfr_323_GA==0 &
                                           dhfr_490_AT==0 ), "IRSI",
                                        ifelse((dhfr_152_AT==0 &
                                                  dhfr_175_TC==1 & 
                                                  dhfr_323_GA==0 &
                                                  dhfr_490_AT==0 ), "NRSI",
                                               ifelse((dhfr_152_AT==1 &
                                                         dhfr_175_TC==0 & 
                                                         dhfr_323_GA==1 &
                                                         dhfr_490_AT==0 ), "ICRI",
                                                      ifelse((dhfr_152_AT==0 &
                                                                dhfr_175_TC==0 & 
                                                                dhfr_323_GA==1 &
                                                                dhfr_490_AT==0 ), "NCRI",
                                                             ifelse((dhfr_152_AT==0 &
                                                                       dhfr_175_TC==0 & 
                                                                       dhfr_323_GA==0 &
                                                                       dhfr_490_AT==0 ), "NCSI",
                                                                    ifelse((dhfr_152_AT==1 &
                                                                              dhfr_175_TC==1 & 
                                                                              dhfr_323_GA==1 &
                                                                              dhfr_490_AT==1 ), "IRNL", "Other")
                                                             )
                                                      )
                                               )
                                        )
                                 )
                          )
                   )
            )
    )
})

snp_calls_t2 %>% print(n=5, width=Inf)
snp_calls_t2 %>% count(dhfr_haplotype)
snp_calls_t2 %>% filter(dhfr_haplotype=="Other")


## dhps
# Options = AGKAA (436A*, 437G*, 540K, 581A, 613A); or AAKAA (436A*, 437A, 540K, 581A, 613A); or SGKAA (436S, 437G*, 540K, 581A, 613A); or SGEAA (436S, 437G*, 540E*, 581A, 613A)
# WT = SAKAA; 3D7 = SGKAA
snp_calls_t2 <- within(snp_calls_t2, {
  dhps_haplotype <-
    ifelse( (dhps_1306_TG==1 &
               dhps_1310_GC==0 & 
               dhps_1618_AG==0 &
               dhps_1742_CG==0 &
               dhps_1837_GT==0), "AGKAA", 
            ifelse( (dhps_1306_TG==1 &
                       dhps_1310_GC==1 & 
                       dhps_1618_AG==0 &
                       dhps_1742_CG==0 &
                       dhps_1837_GT==0), "AAKAA",
                    ifelse( (dhps_1306_TG==0 &
                               dhps_1310_GC==0 & 
                               dhps_1618_AG==0 &
                               dhps_1742_CG==0 &
                               dhps_1837_GT==0), "SGKAA", 
                            ifelse( (dhps_1306_TG==0 &
                                       dhps_1310_GC==1 & 
                                       dhps_1618_AG==0 &
                                       dhps_1742_CG==0 &
                                       dhps_1837_GT==0), "SAKAA",
                                    ifelse( (dhps_1306_TG==0 &
                                               dhps_1310_GC==1 & 
                                               dhps_1618_AG==0 &
                                               dhps_1742_CG==0 &
                                               dhps_1837_GT==1), "SAKAS",
                                            ifelse( (dhps_1306_TG==0 &
                                                       dhps_1310_GC==0 & 
                                                       dhps_1618_AG==0 &
                                                       dhps_1742_CG==0 &
                                                       dhps_1837_GT==1), "SGKAS",
                                                    ifelse( (dhps_1306_TG==1 &
                                                               dhps_1310_GC==0 & 
                                                               dhps_1618_AG==0 &
                                                               dhps_1742_CG==0 &
                                                               dhps_1837_GT==1), "AGKAS",
                                                            ifelse( (dhps_1306_TG==1 &
                                                                       dhps_1310_GC==1 & 
                                                                       dhps_1618_AG==0 &
                                                                       dhps_1742_CG==0 &
                                                                       dhps_1837_GT==1), "AAKAS",
                                                                    ifelse( (dhps_1306_TG==0 &
                                                                               dhps_1310_GC==0 & 
                                                                               dhps_1618_AG==0 &
                                                                               dhps_1742_CG==1 &
                                                                               dhps_1837_GT==0), "SGKGS",
                                                                            ifelse( (dhps_1306_TG==0 &
                                                                                       dhps_1310_GC==0 & 
                                                                                       dhps_1618_AG==1 &
                                                                                       dhps_1742_CG==0 &
                                                                                       dhps_1837_GT==0), "SGEAA", "Other")
                                                                    )
                                                            )
                                                    )
                                            )
                                    )
                            )
                    )
            )
    )
})
snp_calls_t2 %>% print(n=5, width=Inf)
snp_calls_t2 %>% count(dhps_haplotype)
snp_calls_t2 %>% filter(dhps_haplotype=="Other")


## dhfr + dhps
# Options = dhfr-IRNI + dhps-AGKAA (436A*, 437G*, 540K, 581A, 613A); or AAKAA (436A*, 437A, 540K, 581A, 613A); or SGKAA (436S, 437G*, 540K, 581A, 613A)
snp_calls_t2 %>% count(dhfr_haplotype, dhps_haplotype)

snp_calls_t2 <- within(snp_calls_t2, {
  dhfr_dhps_haplotype <-
    paste0( "dhfr-", dhfr_haplotype, ", dhps-", dhps_haplotype)
})

snp_calls_t2 %>% print(n=23, width=Inf)
snp_calls_t2 %>% count(dhfr_dhps_haplotype) %>% arrange(desc(n))


########### GENETIC REPORT CARD

# Artemisinin = ART
# Chloroquine = CQ
# Pyrimethamine = PYR
# Sulfadoxine = SX
# https://www.wwarn.org/sites/default/files/drug-abbreviations.pdf

snp_calls_t3 <- within(snp_calls_t2, {
  
  CQ <- ifelse(crt_227_AC==1,"R",
               ifelse(crt_227_AC==0,"S","Undetermined"))
  
  PYR <- ifelse(dhfr_323_GA==1,"R",
                ifelse(dhfr_323_GA==0,"S","Undetermined"))
  
  SX <- ifelse(dhps_1310_GC==0,"R",
               ifelse(dhps_1310_GC==1,"S","Undetermined"))
  
  SP.Rx <- ifelse( ((dhfr_152_AT==1) & (dhfr_175_TC==1) & (dhfr_323_GA==1)) ,"R",
                   ifelse(((dhfr_152_AT==0) | (dhfr_175_TC==0) | (dhfr_323_GA==0)) ,"S",
                          "Undetermined")
  )
  
  SP.IPTp <- ifelse( ((dhfr_152_AT==1) & (dhfr_175_TC==1) & (dhfr_323_GA==1)) &
                       ((dhps_1310_GC==0) & (dhps_1618_AG==1)) &
                       ((dhfr_490_AT==1) | (dhps_1742_CG==1) | (dhps_1837_GT==1) | (dhps_1837_GA==1)), "R",
                     "S")
})

# Sort out ART-R
if (!exists("k13_snps_allsamps_df")) {
  k13_snps_allsamps_df <- tibble::tibble(
    sample_id = character(),
    snp_id = character(),
    Gene = character(),
    Pos = integer(),
    Ref = character(),
    Alt = character(),
    Qual = numeric(),
    DP = numeric(),
    genotype = character(),
    mut_ARTR_phenotype = logical()
  )
}


if (nrow(k13_snps_allsamps_df) > 0) {
  k13_snps_who_artr <- k13_snps_allsamps_df %>%
    filter(mut_ARTR_phenotype %in% "TRUE") %>%
    filter(sample_id != "Control_KH2")
} else {
  k13_snps_who_artr <- tibble::tibble(
    sample_id = character(),
    snp_id = character(),
    Gene = character(),
    Pos = integer(),
    Ref = character(),
    Alt = character(),
    Qual = numeric(),
    DP = numeric(),
    genotype = character(),
    mut_ARTR_phenotype = logical()
  )
}

snp_calls_t3$ART <- ifelse( (nrow(vcf_tbl_all_k13)==0) | (nrow(k13_snps_who_artr)==0), "S", ">=1 mutation detected")

# Check it out
snp_calls_t3 %>% print(n=5, width=Inf)





############### create genotype format for csp

### Limit to CTD
vcf_tbl_csp_ctd <- vcf_tbl_all_csp %>%
  filter((Pos >= 816) & (Pos <= 1143))  # end of central repeat region to primer


# Create integer for each group
vcf_tbl_csp_ctd <- vcf_tbl_csp_ctd %>%
  mutate(sample_integer = cumsum(!duplicated(sample_id)))


# loop through each sample integer to make list of genotype calls to add each one as a column to the key SNPs

maxn <- max(vcf_tbl_csp_ctd$sample_integer)

mylist_csp_vars <- list()

for (t in 1:maxn) {
  
  sub <- vcf_tbl_csp_ctd %>%
    filter(sample_integer==t) %>%
    select(c("snp_id", "sample_id", "genotype"))
  
  ncolname <- sub[1,"sample_id"]
  
  colnames(sub) <- c("snp_id","sample_id",ncolname)
  
  sub <- sub %>% select(c("snp_id", ncolname))
  
  mylist_csp_vars[[t]] <- sub
  
}

## Join everything together

if (length(mylist_csp_vars) > 0) {
  csp_samp_calls <- mylist_csp_vars[[1]]  # Start with the first element
  if (length(mylist_csp_vars) > 1) {
    for (i in 2:length(mylist_csp_vars)) {  # Loop through the remaining elements
      csp_samp_calls <- full_join(csp_samp_calls, mylist_csp_vars[[i]], by = "snp_id")
    }
  }
} else {
  csp_samp_calls <- tibble::tibble(snp_id = character())  # Create an empty tibble
}


## replace NA with 0
csp_samp_calls <- csp_samp_calls %>%
  replace(is.na(.), 0) %>%
  arrange(snp_id)

head(csp_samp_calls)



############### Save outputs
outdir <- file.path(path_setwd,analysis_dr, MinION_run_name, genotyping_method, "genotypes")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

print("Output directory:")
print(outdir)


# Core drug resistance SNP info in various formats
write.csv(vcf_varcount_df4, file=out_variant_info_fn, row.names=FALSE)
write.csv(vcf_tbl_all_varinfo2, file=out_variant_info2_fn, row.names=FALSE)
write.csv(keysnps_allsamps_gt, file=out_keysnp_allsamp_genotypes_fn, row.names=FALSE)
write.csv(snp_calls_t3, file=out_keysnp_allsamp_genotypes_samprows_fn, row.names=FALSE)

# kelch13 variant info
if(nrow(vcf_tbl_all_k13)>0) {
  write.csv(k13_snps_allsamps_df, file=out_k13_calls_info_fn, row.names=FALSE)
}

# csp nref genotypes
write.csv(vcf_tbl_all_csp, file=out_csp_genotypes_long_fn, row.names=FALSE)
write.csv(csp_samp_calls, file=out_csp_genotypes_samps_fn, row.names=FALSE)

