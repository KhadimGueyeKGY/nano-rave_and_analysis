#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

def extract_gene_indels(input_dir, gene):
    vcf_files = [f for f in os.listdir(input_dir) if gene in f and f.endswith(".vcf")]
    indel_list = []
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
                    if alt == "." or alt == "" or ref == "." or ref == "":
                        continue
                    if len(ref) == 1 and len(alt) == 1:
                        continue  # Skip SNPs
                    try:
                        pos_int = int(pos)
                        indel = (pos_int, f"{ref}-{pos}-{alt}")
                        indel_list.append(indel)
                    except ValueError:
                        continue
    return Counter(indel_list), len(vcf_files)

def plot_indel_frequencies(indel_counts, output_plot_file, gene, sample_count):
    if not indel_counts:
        print(f"No INDELs found for gene {gene}.")
        return None

    df = pd.DataFrame(indel_counts.items(), columns=["Position_INDEL", "Count"])
    df[['Position', 'INDEL']] = pd.DataFrame(df['Position_INDEL'].tolist(), index=df.index)
    df = df.sort_values(by="Position")
    df["Frequency"] = df["Count"] / sample_count

    x_positions = df["Position"]
    x_labels = df["INDEL"]
    y_freqs = df["Frequency"]

    min_x = max(min(x_positions) - 10, 0)
    fig_height = 6
    fig_width = max(10, len(df) * 0.4)

    plt.figure(figsize=(fig_width, fig_height))
    plt.bar(x_positions, y_freqs, color='lightcoral', width=3)
    plt.xticks(ticks=x_positions, labels=x_labels, rotation=90, fontsize=6)
    plt.yticks(fontsize=12)

    plt.xlabel("Genomic Position (INDEL)", fontsize=16)
    plt.ylabel("Frequency\n", fontsize=16)
    plt.title(f"INDEL Frequencies for {gene}", fontsize=16)
    plt.xlim(left=min_x)
    sns.despine(top=True, right=True)

    plt.tight_layout(rect=[0.02, 0.01, 0.99, 0.99])
    plt.savefig(output_plot_file, dpi=300)
    #plt.show()
    print(f"✅ Plot saved to {output_plot_file}")
    plt.close()

    return df[["Position", "INDEL", "Frequency"]]

def main():
    parser = argparse.ArgumentParser(description="Extract and plot INDEL frequencies for a gene")
    parser.add_argument("-e", "--gene", required=True)
    parser.add_argument("-m", "--minion_run_name", required=True)
    parser.add_argument("-g", "--genotyping_method", required=True)
    parser.add_argument("-d", "--analysis_date", required=True)
    parser.add_argument("-a", "--analysis_dr", required=True)

    args = parser.parse_args()
    input_dir = os.path.join(args.analysis_dr, args.minion_run_name, args.genotyping_method, "variant_calling_unzip")
    output_dir = os.path.join(args.analysis_dr, args.minion_run_name, args.genotyping_method, "plots")
    os.makedirs(output_dir, exist_ok=True)

    print(f"🔍 Extracting INDELs from {args.gene}...")
    indel_counts, sample_count = extract_gene_indels(input_dir, args.gene)
    print(f"🧬 Extracted {len(indel_counts)} INDELs from {sample_count} VCF files.")

    output_plot_file = os.path.join(
        output_dir,
        f"{args.minion_run_name}_{args.genotyping_method}_{args.gene}_indel_frequencies_{args.analysis_date}.png"
    )

    indel_df = plot_indel_frequencies(indel_counts, output_plot_file, args.gene, sample_count)

    if indel_df is not None:
        indel_df.insert(0, "Gene", args.gene)
        output_csv_file = os.path.join(
            output_dir,
            f"{args.minion_run_name}_{args.genotyping_method}_{args.gene}_indel_table_{args.analysis_date}.csv"
        )
        indel_df.to_csv(output_csv_file, index=False)
        print(f"📄 INDEL table saved to {output_csv_file}")

if __name__ == "__main__":
    main()
