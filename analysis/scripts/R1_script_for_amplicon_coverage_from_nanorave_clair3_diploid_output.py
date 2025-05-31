import os
import pandas as pd
import numpy as np
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

bedfile_dir = f"{analysis_dr}{MinION_run_name}/{genotyping_method}/genome_coverage"
amplicon_positions_file = "../data/analysis_dependencies/amplicon_gene_positions.csv"
metadata_file = "../data/metadata/metadata_ont_multiplex_barcodes.xlsx"
min_cov_threshold = 50

# Load amplicon positions and metadata
amp_pos_raw = pd.read_csv(amplicon_positions_file)
meta = pd.read_excel(metadata_file, sheet_name="multiplex_B2")
samp_names = meta[["ont_barcode", "sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date"]]

# Function to extract barcode from filename
def extract_barcode(filename):
    barcode = next((item for item in filename.split("_") if "barcode" in item), None)
    return barcode

# Function to process gene coverage
def process_gene_coverage(gene):
    filenames = [f for f in os.listdir(bedfile_dir) if f.endswith(f"_{gene}.bedGraph")]
    if not filenames:
        print(f"⚠️ No BED files found for {gene}, skipping...")
        return None

    gene_coverages = []

    for filename in filenames:
        barcode = extract_barcode(filename)
        file_path = os.path.join(bedfile_dir, filename)
        num_barcode = filename.split("_")
        num_barcode = next((item for item in num_barcode if "barcode" in item), None)
        print(f"🔄 {num_barcode} - Processing coverage for {gene}...")
        
        df = pd.read_csv(file_path, sep="\t", header=None, names=["gene_id", "start", "end", "depth"])
        df = df.astype({"start": int, "end": int, "depth": float})

        amp_pos_gene = amp_pos_raw[amp_pos_raw["gene"] == gene]
        if amp_pos_gene.empty:
            print(f"⚠️ Amplicon positions not found for {gene}, skipping...")
            continue
        amp_pos_gene = amp_pos_gene.iloc[0]

        df_filtered = df[(df["start"] > amp_pos_gene["amplicon_start"]) & 
                         (df["end"] < amp_pos_gene["amplicon_end"])]
        # Check if filtered data is empty before processing
        if df_filtered.empty:
            print(f"⚠️ No data found for barcode {barcode} - gene {gene}, skipping...")
            continue
        
        df_filtered.loc[:, "depth"] = df_filtered["depth"].dropna()

        # Compute summary statistics
        stats = df_filtered["depth"].describe()

        sample_info = samp_names[samp_names["ont_barcode"] == barcode]
        if sample_info.empty:
            print(f"⚠️ Metadata not found for camgen {barcode} - gene {gene}, skipping...")
            continue
        row = sample_info.iloc[0].to_dict()

        coverage_data = {
            "sample_id": row["sample_id"],
            "patient_id": row["patient_id"],
            "ont_multiplex_group": row["ont_multiplex_group"],
            "ont_seq_date": row["ont_seq_date"],
            "ont_barcode": barcode,
            "gene_target": gene,
            "coverage_min": stats["min"],
            "coverage_Q1": stats["25%"],
            "coverage_median": stats["50%"],
            "coverage_mean": stats["mean"],
            "coverage_Q3": stats["75%"],
            "coverage_max": stats["max"],
            "coverage_above_threshold": stats["50%"] > min_cov_threshold,
            "coverage_below_5pcnt": stats["50%"] < np.percentile(df_filtered["depth"], 5),
            "coverage_below_10pcnt": stats["50%"] < np.percentile(df_filtered["depth"], 10)
        }
        gene_coverages.append(coverage_data)

    return pd.DataFrame(gene_coverages)

# Process all genes
target_genes = ["crt", "dhfr", "dhps", "mdr1", "k13", "csp"]

coverage_dfs = [process_gene_coverage(gene) for gene in target_genes if process_gene_coverage(gene) is not None]
coverage_df = pd.concat(coverage_dfs, ignore_index=True)

# Generate Amplicon Run Summary
amplicon_summary = coverage_df.groupby("gene_target")[["coverage_median", "coverage_Q1", "coverage_Q3"]].median()
amplicon_summary["coverage_IQR"] = amplicon_summary["coverage_Q3"] - amplicon_summary["coverage_Q1"]
amplicon_summary["coverage_lower_5pcnt"] = coverage_df.groupby("gene_target")["coverage_median"].quantile(0.05)
amplicon_summary["coverage_lower_10pcnt"] = coverage_df.groupby("gene_target")["coverage_median"].quantile(0.10)

# Generate Per Sample Coverage Report
sample_coverage = coverage_df.pivot(index="sample_id", columns="gene_target", values=["coverage_median", "coverage_below_5pcnt", "coverage_below_10pcnt", "coverage_above_threshold"])
sample_coverage.columns = ["_".join(col) for col in sample_coverage.columns]

# Save both files
output_dir = f"{analysis_dr}{MinION_run_name}/{genotyping_method}/coverage"
os.makedirs(output_dir, exist_ok=True)

amplicon_summary.to_csv(f"{output_dir}/{MinION_run_name}_{genotyping_method}_coverage_amplicon_run_summary_{analysis_date}.csv")
sample_coverage.to_csv(f"{output_dir}/{MinION_run_name}_{genotyping_method}_coverage_by_run_sample{analysis_date}.csv")

print("✅ Coverage summary generated successfully!\n\n")

print(f"✅ Saving files to:  {output_dir}/{MinION_run_name}_{genotyping_method}_coverage_amplicon_run_summary_{analysis_date}.csv")
print(f"✅ Saving files to:  {output_dir}/{MinION_run_name}_{genotyping_method}_coverage_by_run_sample_{analysis_date}.csv")

