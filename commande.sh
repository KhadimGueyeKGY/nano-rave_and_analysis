#!/bin/bash
# reference prep

cd /mnt/c/Users/gueye/Documents/works/Francis-job/nano-rave/data/inputs/reference 

cp ../../../../nanopore/nanopore/analysis/nanorave_files/*.csv .
cp ../../../../nanopore/nanopore/analysis/nanorave_files/*.txt .
cp ../../../../nanopore/nanopore/analysis/nanorave_files/*.fasta .

cd target_gene_cds_seqs/

cat ref_target_gene_cds_seq_crt.fasta ref_target_gene_cds_seq_dhfr.fasta > temp.fasta


============================================================================================================================================================


# Nano-rave

## 1. Medaka_haploid

nextflow run main.nf --sequencing_manifest data/inputs/test_manifest.csv --reference_manifest data/inputs/reference_manifest.csv --variant_caller medaka_haploid --min_barcode_dir_size 5 --results_dir output/medaka_haploid -with-trace

### Clean up files after variant calling

mv medaka_haploid/variant_calling/ medaka_haploid/variant_calling_unzip
rm -r medaka_haploid/variant_calling_unzip/*.tbi
gunzip medaka_haploid/variant_calling_unzip/*.vcf.gz


------------------------------------------------------------------------------------------------------------------------------------------------------------
## 2. Clair3 haploid_precise

nextflow run main.nf --sequencing_manifest data/inputs/test_manifest.csv --reference_manifest data/inputs/reference_manifest.csv --variant_caller clair3 --clair3_args "--model_path /opt/models/r941_prom_sup_g5014 --no_phasing_for_fa --include_all_ctgs --haploid_precise" --results_dir output/clair3_haploid_precise -with-trace

### Clean up files after variant calling

cd output/clair3_haploid_precise
cp -r variant_calling/vcf variant_calling_unzip
rm -r variant_calling_unzip/*.tbi
gunzip variant_calling_unzip/*.vcf.gz
cd ../../

------------------------------------------------------------------------------------------------------------------------------------------------------------
## 3. Clair3 haploid_sensitive

nextflow run main.nf --sequencing_manifest data/inputs/test_manifest.csv --reference_manifest data/inputs/reference_manifest.csv --variant_caller clair3 --clair3_args "--model_path /opt/models/r941_prom_sup_g5014 --no_phasing_for_fa --include_all_ctgs --haploid_sensitive" --results_dir output/clair3_haploid_sens_v2 -with-trace

### Clean up files after variant calling

cd output/clair3_haploid_sens_v2
cp -r variant_calling/vcf variant_calling_unzip
rm -r variant_calling_unzip/*.tbi
gunzip variant_calling_unzip/*.vcf.gz
cd ../../

------------------------------------------------------------------------------------------------------------------------------------------------------------
## 4. Clair3 - diploids

nextflow run main.nf --sequencing_manifest data/inputs/test_manifest.csv --reference_manifest data/inputs/reference_manifest.csv --variant_caller clair3 --clair3_args "--model_path /opt/models/r941_prom_sup_g5014 --no_phasing_for_fa --include_all_ctgs" --results_dir output/clair3_diploid -with-trace

### Clean up files after variant calling

cd output/clair3_diploid
cp -r variant_calling/vcf variant_calling_unzip
rm -r variant_calling_unzip/*.tbi
gunzip variant_calling_unzip/*.vcf.gz
cd ../../


============================================================================================================================================================

# analysis

cd analysis/scripts

Rscript 'GAMP_22906B2 - R script for amplicon coverage, from nanorave Clair3 output (example).R'

Rscript 'R script for coverage plot (MinION run B2).R'

