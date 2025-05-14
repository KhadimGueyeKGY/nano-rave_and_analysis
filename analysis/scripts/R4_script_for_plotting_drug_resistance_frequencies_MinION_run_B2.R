############### R script for drug resistance genotype frequencies

### Set wd
path_setwd <- getwd()
setwd(path_setwd)


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


### Load packages
library("tidyverse")
library("readxl")


# Output path
analysis_dr <- paste0(analysis_dr,'/',MinION_run_name,'/',genotyping_method)

### Output filename
plot_dhfr_dhps_haplo_fn <- paste0(analysis_dr, "/plots", "/",MinION_run_name,"_",genotyping_method,"_dhfr_dhps_haplotype_counts_summary_qc_", analysis_date,".png")
plot_GRC_fn <- paste0(analysis_dr, "/plots", "/",MinION_run_name,"_",genotyping_method,"_drug_resistance_counts_summary_qc_", analysis_date,".png")


### Import data
# define paths
run_cov_fn <- paste0(analysis_dr,"/coverage/",MinION_run_name,"_",genotyping_method,"_coverage_by_run_sample",analysis_date,".csv")

# import data
cov_data <- read.csv(run_cov_fn)

## Genotypes for each sample within each MinION run
# define paths
run_geno_fn <- paste0(analysis_dr,"/genotypes/",MinION_run_name,"_genotype_calls_haplotypes_samplerows_v2.csv")

# import data
geno_raw <- read.csv(run_geno_fn)


### Step 1 - Remove all variants below 50x coverage

# Merge the coverage data into the genotype data
geno_cov <- merge(geno_raw, cov_data,
                    by="sample_id",
                    all.x=T,all.y=T)
geno_cov
dim(geno_cov)


### Replace genotypes with NA if coverage too low for that amplicon for that run

## Run
# crt
geno_cov_v2 <- within(geno_cov, {
  crt_227_AC = ifelse(crt_coverage_above_threshold=="FALSE",NA,crt_227_AC)

# dhfr
  dhfr_152_AT = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_152_AT)
  dhfr_175_TC = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_175_TC)
  dhfr_323_GA = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_323_GA)
  dhfr_323_GC = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_323_GC)
  dhfr_490_AT = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_490_AT)
  dhfr_492_AG = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_492_AG)

# dhps
  dhps_1306_TG = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1306_TG)
  dhps_1307_CT = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1307_CT)
  dhps_1618_AG = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1618_AG)
  dhps_1620_AT = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1620_AT)
  dhps_1742_CG = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1742_CG)
  dhps_1837_GT = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1837_GT)
  dhps_1837_GA = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1837_GA)

# mdr1
  mdr1_256_AT = ifelse(mdr1_coverage_above_threshold=="FALSE",NA,mdr1_256_AT)
  mdr1_551_AT = ifelse(mdr1_coverage_above_threshold=="FALSE",NA,mdr1_551_AT)

# dhfr_haplotype
  dhfr_haplotype = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_haplotype)

# dhps_haplotype
  dhps_haplotype = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_haplotype)

# dhps-dhfr haplotype
  dhfr_dhps_haplotype = ifelse((dhfr_coverage_above_threshold=="FALSE" | 
                                 dhps_coverage_above_threshold=="FALSE"),NA,dhfr_dhps_haplotype)

# GRC
  CQ = ifelse(crt_coverage_above_threshold=="FALSE",NA,CQ)
  PYR = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,PYR)
  SX = ifelse(dhps_coverage_above_threshold=="FALSE",NA,SX)
  SP.Rx = ifelse((dhfr_coverage_above_threshold=="FALSE" |
                      dhps_coverage_above_threshold=="FALSE"),NA,SP.Rx)
  SP.IPTp = ifelse((dhfr_coverage_above_threshold=="FALSE" |
                      dhps_coverage_above_threshold=="FALSE"),NA,SP.IPTp)
  ART = ifelse(k13_coverage_above_threshold=="FALSE",NA,ART)
})


n_distinct(geno_cov_v2$sample_id)
  

## Shed all the coverage stats
colnames(geno_cov_v2)
cols <- c("sample_id",
          "crt_227_AC",
          "dhfr_152_AT",
          "dhfr_175_TC",
          "dhfr_323_GA",
          "dhfr_323_GC",
          "dhfr_490_AT",
          "dhfr_492_AG",
          "dhps_1306_TG",
          "dhps_1307_CT",
          "dhps_1618_AG",
          "dhps_1620_AT",
          "dhps_1742_CG",
          "dhps_1837_GT",
          "dhps_1837_GA",
          "mdr1_256_AT",
          "mdr1_551_AT",
          "dhfr_haplotype",
          "dhps_haplotype",
          "dhfr_dhps_haplotype",
          "CQ",
          "PYR",
          "SX",
          "SP.Rx",
          "SP.IPTp",
          "ART")
geno_cov_v3 <- geno_cov_v2 %>% select(all_of(cols))

head(geno_cov_v3)


### Highlight controls
geno_cov_v4 <- within(geno_cov_v3, {
  control <- ifelse((sample_id=="Control_HB3" |
                       sample_id=="Control_Dd2" | 
                       sample_id=="Control_KH2"), "TRUE", "FALSE")
})
geno_cov_v4 %>% count(control)


geno_qc_controls <- geno_cov_v4 %>%
  filter(control=="TRUE")

geno_qc_nocontrols <- geno_cov_v4 %>%
  filter(control=="FALSE")

dim(geno_qc_nocontrols)



############ DATA VIZ

##### Plot 1 - dhfr and dhps genotypes
### Data manipulations

colnames(geno_qc_nocontrols)

geno_qc_nocontrols$crt_227_AC <- as.character(geno_qc_nocontrols$crt_227_AC)
geno_qc_nocontrols$mdr1_256_AT <- as.character(geno_qc_nocontrols$mdr1_256_AT)
geno_qc_nocontrols$mdr1_551_AT <- as.character(geno_qc_nocontrols$mdr1_551_AT)

cols_geno_plot <- c("sample_id",
                    "dhfr_haplotype",
                    "dhps_haplotype",
                    "dhfr_dhps_haplotype")

cols_GRC_plot <- c("sample_id",
                   "CQ",
                   "SX",
                   "PYR",
                   "SP.Rx",
                   "SP.IPTp",
                   "ART")

data_geno_plot <- geno_qc_nocontrols %>%
  select(all_of(cols_geno_plot)) %>%
  rename(DHFR = dhfr_haplotype,
         DHPS = dhps_haplotype,
         DHFR.DHPS = dhfr_dhps_haplotype) %>%
  mutate(DHFR.DHPS = str_replace(DHFR.DHPS, "dhfr-IRNI dhps-AGKAA", "IRNI.AGKAA")) %>%
  mutate(DHFR.DHPS = str_replace(DHFR.DHPS, "dhfr-IRNI dhps-SGEAA", "IRNI.SGEAA")) %>%
  mutate(DHFR.DHPS = str_replace(DHFR.DHPS, "dhfr-IRNI dhps-AAKAA", "IRNI.AAKAA")) %>%
  mutate(DHFR.DHPS = str_replace(DHFR.DHPS, "dhfr-IRNI dhps-SGKAA", "IRNI.SGKAA")) %>%
  mutate(DHFR.DHPS = str_replace(DHFR.DHPS, "dhfr-NRNI dhps-AGKAA", "NRNI.AGKAA")) %>%
  mutate(DHFR.DHPS = str_replace(DHFR.DHPS, "dhfr-NRNI dhps-AAKAA", "NRNI.AAKAA"))

data_GRC_plot <- geno_qc_nocontrols %>%
  select(all_of(cols_GRC_plot))


##### Genotypes

### prepare data for plotting
data_geno_plot2 <- data_geno_plot %>%
  gather(Genotype, Status, -c(sample_id)) %>% # go from wide to long format
  group_by(Genotype) %>%
  count(Status) %>%
  rename(Count = n) %>%
  filter(Genotype != "DHFR.DHPS")
data_geno_plot2 <- drop_na(data_geno_plot2)
head(data_geno_plot2,10)

### Control plot variables
unique(data_geno_plot2$Status)

legend_labels <- c("IRNI", # dhfr triple mutant
                   "NRNI", # dhfr double mutant
                   "NCSI",
                   "AGKAA", # dhps double mutant
                   "SGKAA", # dhps single mutant
                   "AGKAS", # dhps triple mutant
                   "AAKAA")

break_ord <- c("IRNI",
               "NRNI",
               "NCSI",
               "AGKAA",
               "SGKAA",
               "AGKAS",
               "AAKAA")

order <- c("DHFR",
           "DHPS")

### Generate the plots

unique(data_geno_plot2$Status)
n_distinct(data_geno_plot2$Status)

p1_1 <- data_geno_plot2 %>%
  ggplot(aes(x=Genotype, y=Count, fill=Status)) +
  scale_x_discrete(limits = order) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") +
  ylab("Count") + 
  ylim(0,22) +
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=30, hjust=1),
        axis.title=element_text(size=11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = 'white', colour = NA),
        plot.background = element_rect(fill="white", colour = NA)) +
  scale_fill_brewer(palette = "Set2", name = "Genotype", labels = legend_labels, breaks = break_ord)
p1_1


##### Drug resistance predictions
data_GRC_plot2 <- data_GRC_plot %>%
  gather(Drug, Status, -c(sample_id)) %>% # go from wide to long format
  group_by(Drug) %>%
  count(Status) %>%
  rename(Count = n)

data_GRC_plot2 <- drop_na(data_GRC_plot2)

head(data_GRC_plot2, 10)

### Control plot variables
legend_labels_dr <- c("Sensitive",
                      "Resistant")

break_ord_dr <- c("S",
                  "R")

order_dr <- c("CQ",
              "SX",
              "PYR",
              "SP.Rx",
              "SP.IPTp",
              "ART")

cols_dr <- c("#0072B2", # Blue; Sensitive
             "#F0E442") # Yellow; Resistant
# If need an 'Unknown' then use Green, #009E73

p2_1 <- data_GRC_plot2 %>%
  ggplot(aes(x=Drug, y=Count, fill=Status)) +
  scale_x_discrete(limits = order_dr) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") +
  ylab("Count") + 
  ylim(0,22) +
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=30, hjust=1),
        axis.title=element_text(size=11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = 'white', colour = NA),
        plot.background = element_rect(fill="white", colour = NA)) +
  scale_fill_manual(values = cols_dr, name = "Status", labels = legend_labels_dr, breaks = break_ord_dr)
p2_1

p2_1.2 <- data_GRC_plot2 %>%
  ggplot(aes(x=Drug, y=Count, fill=Status)) +
  scale_x_discrete(limits = order_dr) +
  geom_bar(position="fill", stat="identity") +
  theme_minimal() +
  xlab("") +
  ylab("Ratio") + 
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=30, hjust=1),
        axis.title=element_text(size=11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10) +
  scale_fill_manual(values = cols_dr))
p2_1.2


######### Save outputs
ggsave(plot=p1_1, file=plot_dhfr_dhps_haplo_fn, dpi=200, height=7, width=4, units="in")
ggsave(plot=p2_1, file=plot_GRC_fn, dpi=200, height=4.5, width=7, units="in")

