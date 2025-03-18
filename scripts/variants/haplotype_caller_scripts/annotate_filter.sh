#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msilarge,msismall
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

export BCFTOOLS_PLUGINS=/home/selmecki/shared/software/software.install/bcftools/1.17/libexec/bcftools/

regions=SC5314_A21_regions.bed
species=Calbicans
ref=SC5314-A21
gene_file=Calbicans_SC5314_A21_sorted_genes.bed.gz
gene_vcf=Calbicans_genes.vcf.gz
fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta
bcftools_out=Calbicans_filtered.vcf.gz
snpeff=/home/selmecki/shared/software/snpEff/snpEff.jar
snpeff_config=/home/selmecki/shared/software/snpEff/snpEff.config
snpeff_db=SC5314_s02m09r08
annotate_vcf=Calbicans_bwa_filtered_annotated.vcf
genotype_table=Calbicans_bwa_genotypes.txt

# Load modules
module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load htslib

# Remove intermediate files.
function finish {
  rm "${gene_vcf}"*
  rm "${species}_${ref}"_regions.vcf*
  rm "${bcftools_out}"
 }

 trap finish EXIT

# Filter out repetitive regions, compress, index
# because un-indexed VCFs do not always work
bcftools view \
    -R "${regions}" \
    "${species}"_"${ref}"_merged.vcf.gz \
    -Oz \
    -o "${species}"_"${ref}"_regions.vcf.gz

bcftools index "${species}"_"${ref}"_regions.vcf.gz

# Add gene name annotation as info column.
bcftools annotate \
    -a "${gene_file}" \
    -c CHROM,FROM,TO,GENE \
    -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \
    -Oz \
    -o "${gene_vcf}" \
    "${species}"_"${ref}"_regions.vcf.gz

bcftools index "${gene_vcf}"

# Quality and frequency hard filtering.
bcftools norm -f "${fasta}" -m -indels \
    "${gene_vcf}" \
    | bcftools +fill-tags - -- -t "FORMAT/VAF" \
    | bcftools view -i "(GT='het' & FORMAT/VAF > 0.14 & FORMAT/VAF < 0.85) | (GT='AA' & FORMAT/VAF > 0.98)" \
    | bcftools view -e 'GT[0]="alt"' -Ou \
    | bcftools view -c 1 \
    -Oz -o "${bcftools_out}"

# Annotate using snpeff with manually built database.
java -Xmx9g -jar "${snpeff}" -c "${snpeff_config}" "${snpeff_db}" \
"${bcftools_out}" >  "${annotate_vcf}"

# Zip and index output.
bgzip "${annotate_vcf}"
tabix "${annotate_vcf}".gz

# Subset to biallelic SNPs, output tab-delimited genotype file,
# use for multiple correspondence analysis

bcftools view -m2 -M2 -v snps "${annotate_vcf}" \
    | bcftools view -e 'GT="mis"' \
    | bcftools query -H -f '%CHROM\t%POS[\t%GT]\n' >> "${genotype_table}"

sed -i '1s/\[[0-9]\+\]//g' "${genotype_table}"
sed -i '1s/\:GT//g' "${genotype_table}"
sed -i '1s/^\# //' "${genotype_table}"

