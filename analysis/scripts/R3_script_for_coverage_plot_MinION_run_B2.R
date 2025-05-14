############### R script for coverage plot

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
library("readxl")
library("tidyverse")


# Output path
analysis_dr <- paste0(analysis_dr,'/',MinION_run_name,'/',genotyping_method)

### Output filename
plot_cov_by_amplicon_all_fn <- paste0(analysis_dr, "/plots", "/",MinION_run_name,"_",genotyping_method,"_coverage_plot_",analysis_date,".png")

### Import data
## Coverage for each sample within each MinION run
# define paths
run_cov_fn <- paste0(analysis_dr,"/coverage/",MinION_run_name,"_",genotyping_method,"_coverage_by_run_sample",analysis_date,".csv")

# import data
cov_data <- read.csv(run_cov_fn)

### Annotations for sample type
# Highlight controls
cov_data$is_control <-
  ifelse( (cov_data$sample_id=="Control_Dd2" |
             cov_data$sample_id=="Control_KH2" |
             cov_data$sample_id=="Control_HB3"), "TRUE", "FALSE")

head(cov_data)
cov_data %>% count(is_control)

colnames(cov_data)
cols_to_use <- c("sample_id",
                 "is_control",
                 "crt_coverage_median",
                 "dhfr_coverage_median",
                 "dhps_coverage_median",
                 "mdr1_coverage_median",
                 "k13_coverage_median",
                 "csp_coverage_median")

cov_data_v2 <- cov_data %>%
  select(all_of(cols_to_use))

head(cov_data_v2)


###### Define sample sets
cov_data_all_samps <- cov_data_v2

cov_data_no_controls <- cov_data_v2 %>%
  filter(is_control=="FALSE")

dim(cov_data_all_samps)
dim(cov_data_no_controls)


##### Prepare for plotting
## Prepare the plot features
p_order <- c("crt",
              "dhfr",
              "dhps",
              "mdr1",
              "kelch13",
              "csp")

p_labels <-  c("crt",
               "dhfr",
               "dhps",
               "mdr1",
               "kelch13",
               "csp")


## Prepare data for plotting
cov_data_no_controls_plt <- cov_data_no_controls %>%
  pivot_longer(
    cols = ends_with("coverage_median"),
    values_to = "coverage",
    names_to = "gene_target") %>%
  mutate(gene_target = str_replace(gene_target, "crt_coverage_median", "crt")) %>%
  mutate(gene_target = str_replace(gene_target, "dhfr_coverage_median", "dhfr")) %>%
  mutate(gene_target = str_replace(gene_target, "dhps_coverage_median", "dhps")) %>%
  mutate(gene_target = str_replace(gene_target, "mdr1_coverage_median", "mdr1")) %>%
  mutate(gene_target = str_replace(gene_target, "k13_coverage_median", "kelch13")) %>%
  mutate(gene_target = str_replace(gene_target, "csp_coverage_median", "csp"))

##### Stats

cov_summ_crt <- cov_data_no_controls_plt %>%
  filter(gene_target=="crt")

head(cov_summ_crt)

crt_med_cov <- median(cov_summ_crt$coverage, na.rm=T)
dhfr_med_cov <- median(cov_data_no_controls$dhfr_coverage_median)
dhps_med_cov <- median(cov_data_no_controls$dhps_coverage_median)
mdr1_med_cov <- median(cov_data_no_controls$mdr1_coverage_median)
k13_med_cov <- median(cov_data_no_controls$k13_coverage_median)
csp_med_cov <- median(cov_data_no_controls$csp_coverage_median)

print(paste0("Median coverage for crt: ", crt_med_cov))
print(paste0("Median coverage for dhfr: ", dhfr_med_cov))
print(paste0("Median coverage for dhps: ", dhps_med_cov))
print(paste0("Median coverage for mdr1: ", mdr1_med_cov))
print(paste0("Median coverage for k13: ", k13_med_cov))
print(paste0("Median coverage for csp: ", csp_med_cov))

##### Produce plots
##Excluding controls
p1 <- cov_data_no_controls_plt %>%
  ggplot( aes(x=gene_target, y=coverage)) +
  geom_boxplot(outlier.shape = NA, fill="#e0e0e0") +
  theme_classic() +
  xlab("") +
  ylab("Coverage depth") +
  scale_y_continuous(limits = c(0,20000), expand = c(0, 0)) +
  theme(axis.text.x=element_text(size=11),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=11),
        axis.text=element_text(colour="black")) +
  geom_jitter(width = 0.1, height = 0, size=1.3, alpha=0.6) +
  scale_x_discrete(limits=p_order,labels=p_labels) +
  theme(plot.margin=unit(c(0.6,0.4,0,0.4),"cm")) # t, r, b, l
p1


#### Save plots
# Ensure the directory exists
dir.create(dirname(plot_cov_by_amplicon_all_fn), recursive = TRUE, showWarnings = FALSE)
ggsave(plot=p1, file=plot_cov_by_amplicon_all_fn, dpi=300, height=4, width=7, units="in")

  