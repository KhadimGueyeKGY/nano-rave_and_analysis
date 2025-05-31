import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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
analysis_dir = os.path.join(analysis_dr, MinION_run_name, genotyping_method)
plot_dir = os.path.join(analysis_dir, "plots")
os.makedirs(plot_dir, exist_ok=True)

plot_dhfr_dhps_haplo_fn = os.path.join(plot_dir, f"{MinION_run_name}_{genotyping_method}_dhfr_dhps_haplotype_counts_summary_qc_{analysis_date}.png")
plot_GRC_fn = os.path.join(plot_dir, f"{MinION_run_name}_{genotyping_method}_drug_resistance_counts_summary_qc_{analysis_date}.png")

run_cov_fn = os.path.join(analysis_dir, "coverage", f"{MinION_run_name}_{genotyping_method}_coverage_by_run_sample{analysis_date}.csv")
run_geno_fn = os.path.join(analysis_dir, "genotypes", f"{MinION_run_name}_genotype_calls_haplotypes_samplerows_v2.csv")

cov_data = pd.read_csv(run_cov_fn)
geno_raw = pd.read_csv(run_geno_fn)

geno_cov = pd.merge(geno_raw, cov_data, on="sample_id", how="left")

replace_conditions = {
    "crt_227_AC": "crt_227_coverage_above_threshold",
    "dhfr_152_AT": "dhfr_coverage_above_threshold",
    "dhfr_175_TC": "dhfr_coverage_above_threshold",
    "dhfr_323_GA": "dhfr_coverage_above_threshold",
    "dhfr_323_GC": "dhfr_coverage_above_threshold",
    "dhfr_490_AT": "dhfr_coverage_above_threshold",
    "dhfr_492_AG": "dhfr_coverage_above_threshold",
    "dhps_1306_TG": "dhps_coverage_above_threshold",
    "dhps_1307_CT": "dhps_coverage_above_threshold",
    "dhps_1618_AG": "dhps_coverage_above_threshold",
    "dhps_1620_AT": "dhps_coverage_above_threshold",
    "dhps_1742_CG": "dhps_coverage_above_threshold",
    "dhps_1837_GT": "dhps_coverage_above_threshold",
    "dhps_1837_GA": "dhps_coverage_above_threshold",
    "mdr1_256_AT": "mdr1_coverage_above_threshold",
    "mdr1_551_AT": "mdr1_coverage_above_threshold",
    "dhfr_haplotype": "dhfr_coverage_above_threshold",
    "dhps_haplotype": "dhps_coverage_above_threshold",
    "dhfr_dhps_haplotype": ["dhfr_coverage_above_threshold", "dhps_coverage_above_threshold"],
    "CQ": "crt_227_coverage_above_threshold",
    "PYR": "dhfr_coverage_above_threshold",
    "SX": "dhps_coverage_above_threshold",
    "SP.Rx": ["dhfr_coverage_above_threshold", "dhps_coverage_above_threshold"],
    "SP.IPTp": ["dhfr_coverage_above_threshold", "dhps_coverage_above_threshold"],
    "ART": "k13_coverage_above_threshold"
}

for col, cond in replace_conditions.items():
    if isinstance(cond, list):
        if all(c in geno_cov.columns for c in cond):
            mask = ~(geno_cov[cond[0]] == True) | ~(geno_cov[cond[1]] == True)
        else:
            continue  # Skip if any required column is missing
    else:
        if cond in geno_cov.columns:
            mask = ~(geno_cov[cond] == True)
        else:
            continue  # Skip if required column is missing
    geno_cov[col] = geno_cov[col].mask(mask)

cols = list(replace_conditions.keys()) + ["sample_id"]
geno_cov_v3 = geno_cov[cols].copy()


geno_cov_v3["control"] = geno_cov_v3["sample_id"].isin(["Control_HB3", "Control_Dd2", "Control_KH2"])
geno_qc_nocontrols = geno_cov_v3[geno_cov_v3["control"] == False].copy()

geno_qc_nocontrols["crt_227_AC"] = geno_qc_nocontrols["crt_227_AC"].astype(str)
geno_qc_nocontrols["mdr1_256_AT"] = geno_qc_nocontrols["mdr1_256_AT"].astype(str)
geno_qc_nocontrols["mdr1_551_AT"] = geno_qc_nocontrols["mdr1_551_AT"].astype(str)

data_geno_plot = geno_qc_nocontrols[["sample_id", "dhfr_haplotype", "dhps_haplotype", "dhfr_dhps_haplotype"]].rename(columns={
    "dhfr_haplotype": "DHFR",
    "dhps_haplotype": "DHPS",
    "dhfr_dhps_haplotype": "DHFR.DHPS"
})

replacements = {
    "dhfr-IRNI dhps-AGKAA": "IRNI.AGKAA",
    "dhfr-IRNI dhps-SGEAA": "IRNI.SGEAA",
    "dhfr-IRNI dhps-AAKAA": "IRNI.AAKAA",
    "dhfr-IRNI dhps-SGKAA": "IRNI.SGKAA",
    "dhfr-NRNI dhps-AGKAA": "NRNI.AGKAA",
    "dhfr-NRNI dhps-AAKAA": "NRNI.AAKAA"
}
data_geno_plot["DHFR.DHPS"] = data_geno_plot["DHFR.DHPS"].replace(replacements)

data_GRC_plot = geno_qc_nocontrols[["sample_id", "CQ", "SX", "PYR", "SP.Rx", "SP.IPTp", "ART"]]

data_geno_plot2 = data_geno_plot.melt(id_vars="sample_id", var_name="Genotype", value_name="Status")
data_geno_plot2 = data_geno_plot2[data_geno_plot2["Genotype"] != "DHFR.DHPS"]
data_geno_plot2 = data_geno_plot2.dropna().groupby(["Genotype", "Status"]).size().reset_index(name="Count")

order = ["DHFR", "DHPS"]

plt.figure(figsize=(12, 8))  # Large, wide, readable figure
sns.barplot(data=data_geno_plot2, x="Genotype", y="Count", hue="Status", order=order)
plt.ylim(0)
plt.xlabel("\nGenotype", fontsize=16)
plt.ylabel("Sample Count\n", fontsize=16)
plt.xticks(rotation=30, ha="right", fontsize=14)
plt.yticks(fontsize=14)
plt.legend(title="Genotype", loc="center right", bbox_to_anchor=(1.18, 0.5), ncol=1, fontsize=14, title_fontproperties={"weight": "bold"})
sns.despine(top=True, right=True)
plt.tight_layout()
plt.savefig(plot_dhfr_dhps_haplo_fn, dpi=300, bbox_inches="tight")  # High DPI, tight bounding box
print(f'✅ : {plot_dhfr_dhps_haplo_fn}')
plt.close()

data_GRC_plot2 = data_GRC_plot.melt(id_vars="sample_id", var_name="Drug", value_name="Status")
data_GRC_plot2 = data_GRC_plot2.dropna().groupby(["Drug", "Status"]).size().reset_index(name="Count")

order_dr = ["CQ", "SX", "PYR", "SP.Rx", "SP.IPTp", "ART"]
cols_dr = {"S": "#0072B2", "R": "#F0E442"}
legend_labels_dr = ["Sensitive", "Resistant"]
break_ord_dr = ["S", "R"]

plt.figure(figsize=(16, 8))  # Even wider for drug resistance plot
ax = sns.barplot(
    data=data_GRC_plot2,
    x="Drug",
    y="Count",
    hue="Status",
    order=order_dr,
    hue_order=break_ord_dr,
    palette=cols_dr
)
plt.ylim(0)
plt.xlabel("\nDrug", fontsize=16)
plt.ylabel("Sample Count\n", fontsize=16)
plt.xticks(rotation=30, ha="right", fontsize=14)
plt.yticks(fontsize=14)

# Set legend with custom labels, place at center right, slightly outside
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles,
    legend_labels_dr,
    title="Resistance Status",
    loc="center right",
    bbox_to_anchor=(1.18, 0.5),
    fontsize=14,
    title_fontproperties={"weight": "bold"}
)
sns.despine(top=True, right=True)
plt.tight_layout()
plt.savefig(plot_GRC_fn, dpi=300, bbox_inches="tight")  # High DPI, tight bounding box
print(f'✅ : {plot_GRC_fn}')
plt.close()


print(f"✅ Plots saved")

