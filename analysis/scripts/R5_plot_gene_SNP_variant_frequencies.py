#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter


def extract_gene_variants(input_dir, gene):
    vcf_files = [f for f in os.listdir(input_dir) if gene in f and f.endswith(".vcf")]
    snp_list = []
    for vcf in vcf_files:
        with open(os.path.join(input_dir, vcf)) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    pos = parts[1]
                    ref = parts[3]
                    alt = parts[4]
                    if alt == "." or alt == "" or len(ref) != 1 or len(alt) != 1:
                        continue  # Only SNPs (1 base)
                    try:
                        pos_int = int(pos)
                        snp = (pos_int, ref, alt)
                        snp_list.append(snp)
                    except ValueError:
                        continue
    return Counter(snp_list), len(vcf_files)


def plot_variant_frequencies(variant_counts, output_plot_file, gene, sample_count):
    if not variant_counts:
        print(f"No SNPs found for gene {gene}.")
        return

    df = pd.DataFrame(variant_counts.items(), columns=["SNP_info", "Count"])
    df[["Position", "Ref", "Alt"]] = pd.DataFrame(df["SNP_info"].tolist(), index=df.index)
    df = df.sort_values(by="Position")
    df["Frequency"] = df["Count"] / sample_count
    df["SNP"] = df.apply(lambda row: f"{row['Position']} {row['Ref']}:{row['Alt']}", axis=1)

    x_positions = df["Position"]
    x_labels = df["SNP"]
    y_freqs = df["Frequency"]
    min_x = max(min(x_positions) - 10, 0)

    fig_height = 6
    fig_width = max(10, len(df) * 0.4)
    plt.figure(figsize=(fig_width, fig_height))
    plt.bar(x_positions, y_freqs, color='skyblue', width=3)
    plt.xticks(ticks=x_positions, labels=x_labels, rotation=90, fontsize=10)
    plt.yticks(fontsize=12)

    plt.xlabel("Genomic Position (SNP)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.title(f"SNP Frequencies for {gene}", fontsize=16)
    plt.xlim(left=min_x)
    sns.despine(top=True, right=True)
    plt.tight_layout()
    plt.savefig(output_plot_file, dpi=300)
    #plt.show()
    print(f"✅ Plot saved to {output_plot_file}")
    plt.close()

    return df[["Gene", "SNP", "Frequency"]] if "Gene" in df.columns else df[["Position", "Ref", "Alt", "SNP", "Frequency"]]


def main():
    parser = argparse.ArgumentParser(description="Extract and plot SNP frequencies for a gene")
    parser.add_argument("-e", "--gene", required=True)
    parser.add_argument("-m", "--minion_run_name", required=True)
    parser.add_argument("-g", "--genotyping_method", required=True)
    parser.add_argument("-d", "--analysis_date", required=True)
    parser.add_argument("-a", "--analysis_dr", required=True)

    args = parser.parse_args()
    input_dir = os.path.join(args.analysis_dr, args.minion_run_name, args.genotyping_method, "variant_calling_unzip")
    output_dir = os.path.join(args.analysis_dr, args.minion_run_name, args.genotyping_method, "plots")
    os.makedirs(output_dir, exist_ok=True)

    print(f"🔍 Extracting SNPs from {args.gene}...")
    variant_counts, sample_count = extract_gene_variants(input_dir, args.gene)
    print(f"🧬 Extracted {len(variant_counts)} SNPs from {sample_count} VCF files.")

    output_plot_file = os.path.join(
        output_dir,
        f"{args.minion_run_name}_{args.genotyping_method}_{args.gene}_variant_frequencies_{args.analysis_date}.png"
    )

    snp_df = plot_variant_frequencies(variant_counts, output_plot_file, args.gene, sample_count)

    if snp_df is not None:
        snp_df.insert(0, "Gene", args.gene)
        output_table_file = os.path.join(
            output_dir,
            f"{args.minion_run_name}_{args.genotyping_method}_{args.gene}_snp_variant_table_{args.analysis_date}.csv"
        )
        snp_df[["Gene", "SNP", "Frequency"]].to_csv(output_table_file, index=False)
        print(f"📄 SNP table saved to {output_table_file}")


if __name__ == "__main__":
    main()
