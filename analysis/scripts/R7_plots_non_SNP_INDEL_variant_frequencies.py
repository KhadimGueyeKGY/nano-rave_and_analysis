import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter


def extract_other_variants(input_dir, gene):
    vcf_files = [f for f in os.listdir(input_dir) if gene in f and f.endswith(".vcf")]
    variant_list = []
    for vcf in vcf_files:
        with open(os.path.join(input_dir, vcf)) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue

                pos = parts[1]
                ref = parts[3]
                alt = parts[4]

                if alt in [".", ""] or ref in [".", ""]:
                    continue

                if "," in alt:
                    alt_alleles = alt.split(",")
                else:
                    alt_alleles = [alt]

                for alt_allele in alt_alleles:
                    if not (len(ref) == 1 and len(alt_allele) == 1):  # not SNP
                        if not (len(ref) != len(alt_allele)):  # not INDEL
                            try:
                                pos_int = int(pos)
                                variant = (pos_int, f"{ref}-{pos}-{alt_allele}")
                                variant_list.append(variant)
                            except ValueError:
                                continue
    return Counter(variant_list), len(vcf_files)


def plot_variant_frequencies(variant_counts, output_plot_file, gene, sample_count):
    if not variant_counts:
        print(f"No non-SNP/INDEL variants found for gene {gene}.")
        return

    df = pd.DataFrame(variant_counts.items(), columns=["Position_Variant", "Count"])
    df[['Position', 'Variant']] = pd.DataFrame(df['Position_Variant'].tolist(), index=df.index)
    df = df.sort_values(by="Position")
    df["Frequency"] = df["Count"] / sample_count

    x_positions = df["Position"]
    x_labels = df["Variant"]
    y_freqs = df["Frequency"]
    min_x = max(min(x_positions) - 10, 0)

    fig_height = 6
    fig_width = max(12, len(df) * 0.7)  # Increase width per label for better spacing

    plt.figure(figsize=(fig_width, fig_height))
    plt.bar(x_positions, y_freqs, color='salmon', width=3)
    plt.xticks(ticks=x_positions, labels=x_labels, rotation=90, fontsize=12, ha='center')
    plt.yticks(fontsize=12)

    plt.xlabel("Genomic Position (Other Variant)", fontsize=16)
    plt.ylabel("Frequency\n", fontsize=16)
    plt.title(f"Other Variant Frequencies for {gene}", fontsize=16)
    plt.xlim(left=min_x)
    sns.despine(top=True, right=True)
    plt.tight_layout()  # Let matplotlib auto-adjust the layout
    plt.savefig(output_plot_file, dpi=300, bbox_inches='tight')  # Ensure nothing is cut off
    plt.show()
    print(f"Plot saved to {output_plot_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Extract and plot non-SNP/INDEL variant frequencies")
    parser.add_argument("-e", "--gene", required=True)
    parser.add_argument("-m", "--minion_run_name", required=True)
    parser.add_argument("-g", "--genotyping_method", required=True)
    parser.add_argument("-d", "--analysis_date", required=True)
    parser.add_argument("-a", "--analysis_dr", required=True)

    args = parser.parse_args()
    input_dir = os.path.join(args.analysis_dr, args.minion_run_name, args.genotyping_method, "variant_calling_unzip")
    output_dir = os.path.join(args.analysis_dr, args.minion_run_name, args.genotyping_method, "plots")
    os.makedirs(output_dir, exist_ok=True)

    variant_counts, sample_count = extract_other_variants(input_dir, args.gene)
    print(f"Extracted {len(variant_counts)} non-SNP/INDEL variants for gene {args.gene} from {sample_count} VCFs.")

    output_plot_file = os.path.join(
        output_dir,
        f"{args.minion_run_name}_{args.genotyping_method}_{args.gene}_other_variant_frequencies_{args.analysis_date}.png"
    )

    plot_variant_frequencies(variant_counts, output_plot_file, args.gene, sample_count)


if __name__ == "__main__":
    main()
