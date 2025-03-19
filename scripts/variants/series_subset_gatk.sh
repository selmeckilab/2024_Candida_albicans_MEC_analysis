#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msilarge,msismall
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=2-20

set -ue
set -o pipefail

export BCFTOOLS_PLUGINS=/home/selmecki/shared/software/software.install/bcftools/1.17/plugins

# Input variables
line=${SLURM_ARRAY_TASK_ID}
file_list=series_list.txt
full_vcf=../Calbicans_bwa_filtered_annotated.vcf.gz

# Get file and series variables
file=$(awk -v val="$line" 'NR == val { print $0}' "$file_list")
series=$(basename "$file" | cut -d "_" -f 1)
sample_count=$(wc -l < "${file}")
max_alt=$(( sample_count - 1 ))
subset_vcf="${series}/${series}_unfixed.vcf.gz"
genotype_table="${series}/${series}_genotypes.txt"

mkdir -p "${series}"

# Load modules
module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load htslib

# Subset to only unfixed variants within each series
# and to biallelic SNPs only
bcftools view -S "${file}" -c 1 -C "${max_alt}" "${full_vcf}" \
    | bcftools view -m2 -M2 -v snps \
    | bcftools view -i "(GT='het' & FORMAT/VAF > 0.14 & FORMAT/VAF < 0.85) | (GT='AA' & FORMAT/VAF > 0.98)" \
  -Oz -o "${subset_vcf}"

tabix "${subset_vcf}"

# Output a genotype table of these unfixed vars for use in upset plot
bcftools view -e 'GT="mis"' "${subset_vcf}" \
    | bcftools query -H -f '%CHROM\t%POS[\t%GT]\n' >> "${genotype_table}"

# Clean up genotype table header
sed -i '1s/\[[0-9]\+\]//g' "${genotype_table}"
sed -i '1s/\:GT//g' "${genotype_table}"
sed -i '1s/^\# //' "${genotype_table}"
