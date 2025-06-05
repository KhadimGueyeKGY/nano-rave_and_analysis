#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

def detect_snp_positions(input_dir, gene):
    """Automatically detect all SNP positions from VCF files for a given gene."""
    vcf_files = [f for f in os.listdir(input_dir) if gene in f and f.endswith(".vcf")]
    snp_positions = set()

    for vcf in vcf_files:
        with open(os.path.join(input_dir, vcf)) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                pos, ref, alt = parts[1], parts[3], parts[4]
                if alt == "." or alt == "" or ref == "." or ref == "":
                    continue
                if len(ref) == 1 and len(alt) == 1:  # valid SNP
                    snp_positions.add(pos)
    return sorted(snp_positions)

def extract_haplotypes(input_dir, gene, snp_positions):
    """Extract haplotypes from VCFs based on detected SNP positions."""
    vcf_files = [f for f in os.listdir(input_dir) if gene in f and f.endswith(".vcf")]
    haplotypes = []

    for vcf in vcf_files:
        sample_hap = {}
        with open(os.path.join(input_dir, vcf)) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue
                pos = parts[1]
                ref = parts[3]
                alt = parts[4]
                fmt = parts[8].split(":")
                values = parts[9].split(":")

                if pos not in snp_positions:
                    continue

                if len(ref) != 1 or len(alt) != 1:
                    continue  # skip indels

                fmt_dict = dict(zip(fmt, values))
                gt = fmt_dict.get("GT", "./.")
                if gt == "0/0":
                    allele = ref
                elif gt in ["0/1", "1/0"]:
                    allele = ref + "/" + alt
                elif gt == "1/1":
                    allele = alt
                else:
                    allele = "NA"

                sample_hap[pos] = allele

        hap_string = "-".join([sample_hap.get(pos, "NA") for pos in snp_positions])
        haplotypes.append(hap_string)

    haplo_counter = Counter(haplotypes)
    return haplo_counter, len(vcf_files)

def plot_haplotype_frequencies(haplo_counts, output_plot_file, gene, sample_count):
    if not haplo_counts:
        print(f"No haplotypes found for gene {gene}.")
        return

    df = pd.DataFrame(haplo_counts.items(), columns=["Haplotype", "Count"])
    df = df.sort_values(by="Count", ascending=False).reset_index(drop=True)
    df["Frequency"] = df["Count"] / sample_count
    df["Rank"] = df.index + 1

    plt.figure(figsize=(max(10, len(df) * 0.5), 6))
    sns.barplot(data=df, x="Rank", y="Frequency", hue="Haplotype", dodge=False, palette="Set2")
    plt.xlabel("Haplotype Rank", fontsize=14)
    plt.ylabel("Haplotype Frequency", fontsize=14)
    plt.title(f"Frequency Distribution for {gene} Haplotypes - ONT", fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(title="Haplotype", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_plot_file, dpi=300)
    plt.show()
    plt.close()
    print(f"✅ Plot saved to {output_plot_file}")

def main():
    parser = argparse.ArgumentParser(description="Automatically extract and plot haplotype frequencies from VCFs")
    parser.add_argument("-e", "--gene", required=True, help="Gene name (e.g., mdr1)")
    parser.add_argument("-m", "--minion_run_name", required=True, help="MinION run name")
    parser.add_argument("-g", "--genotyping_method", required=True, help="Genotyping method (e.g., clair3_diploid)")
    parser.add_argument("-d", "--analysis_date", required=True, help="Analysis date (e.g., 20250601)")
    parser.add_argument("-a", "--analysis_dr", required=True, help="Base directory for analysis")

    args = parser.parse_args()

    input_dir = os.path.join(args.analysis_dr, args.minion_run_name, args.genotyping_method, "variant_calling_unzip")
    output_dir = os.path.join(args.analysis_dr, args.minion_run_name, args.genotyping_method, "plots")
    os.makedirs(output_dir, exist_ok=True)

    print(f"🔍 Detecting SNP positions in {args.gene}...")
    snp_positions = detect_snp_positions(input_dir, args.gene)
    print(f"📌 {len(snp_positions)} SNP positions found: {', '.join(snp_positions)}")

    haplo_counts, sample_count = extract_haplotypes(input_dir, args.gene, snp_positions)
    print(f"🧬 Extracted {len(haplo_counts)} haplotypes for gene {args.gene} from {sample_count} VCF(s).")

    output_plot_file = os.path.join(
        output_dir,
        f"{args.minion_run_name}_{args.genotyping_method}_{args.gene}_haplotype_frequencies_{args.analysis_date}.png"
    )

    plot_haplotype_frequencies(haplo_counts, output_plot_file, args.gene, sample_count)

if __name__ == "__main__":
    main()
