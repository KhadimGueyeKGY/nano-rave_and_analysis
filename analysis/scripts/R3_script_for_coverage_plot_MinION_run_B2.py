import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

# Argument parser setup
parser = argparse.ArgumentParser(
    description="Process amplicon coverage from NanoRAVE Clair3 diploid output"
)

parser.add_argument("-m", "--minion_run_name", required=True, help="MinION run name (e.g., csp_setA_18aprl2025)")
parser.add_argument("-g", "--genotyping_method", required=True, help="Genotyping method (e.g., clair3_diploid)")
parser.add_argument("-d", "--analysis_date", required=True, help="Analysis date (format: YYYYMMDD)")
parser.add_argument("-a", "--analysis_dr", required=True, help="Path to the output directory")

args = parser.parse_args()

# Assign variables
MinION_run_name = args.minion_run_name
genotyping_method = args.genotyping_method
analysis_date = args.analysis_date
analysis_dr = args.analysis_dr

print("📦 Running coverage script with:")
print(f"   MinION run name     : {MinION_run_name}")
print(f"   Genotyping method   : {genotyping_method}")
print(f"   Analysis date       : {analysis_date}")
print(f"   Output directory    : {analysis_dr}")

# Define paths
analysis_path = f"{analysis_dr}/{MinION_run_name}/{genotyping_method}"
coverage_file = f"{analysis_path}/coverage/{MinION_run_name}_{genotyping_method}_coverage_by_run_sample{analysis_date}.csv"
output_plot = f"{analysis_path}/plots/{MinION_run_name}_{genotyping_method}_coverage_plot_{analysis_date}.png"

# Ensure plot directory exists
os.makedirs(os.path.dirname(output_plot), exist_ok=True)

# Load coverage data
cov_data = pd.read_csv(coverage_file)

# Annotate controls
control_samples = ["Control_Dd2", "Control_KH2", "Control_HB3"]
cov_data["is_control"] = cov_data["sample_id"].apply(lambda x: "TRUE" if x in control_samples else "FALSE")

# Identify all gene coverage columns dynamically
coverage_cols = [col for col in cov_data.columns if col.startswith("coverage_median_")]

if not coverage_cols:
    raise ValueError("No coverage_median_* columns found in the input file.")

# Melt the dataframe for seaborn
cov_long = cov_data.melt(id_vars=["sample_id", "is_control"], value_vars=coverage_cols,
                         var_name="gene_target", value_name="coverage")

# Clean gene target names
cov_long["gene_target"] = cov_long["gene_target"].str.replace("coverage_median_", "", regex=False)

# Exclude controls
cov_plot_data = cov_long[cov_long["is_control"] == "FALSE"]

# Plot configuration
sns.set_style("ticks")  # removes grid
plt.figure(figsize=(10, 6))

# Boxplot
sns.boxplot(x="gene_target", y="coverage", data=cov_plot_data, palette="Blues", showfliers=False)

# Jitter plot
sns.stripplot(x="gene_target", y="coverage", data=cov_plot_data, color="black", alpha=0.6, jitter=True)

# Labels & custom styling
plt.xlabel("\nGene Target", fontsize=16)
plt.ylabel("Coverage Depth\n", fontsize=16)
#plt.title("Coverage Depth per Gene") #, fontsize=14)
plt.ylim(0)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
sns.despine(top=True, right=True)
plt.tight_layout()

# Save the plot
plt.savefig(output_plot, dpi=300)
print(f"✅ Plot saved at: {output_plot}")

# Show the plot
plt.show()
