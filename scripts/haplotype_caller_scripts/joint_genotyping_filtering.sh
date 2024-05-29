#!/bin/bash
#SBATCH --mem=32gb
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=48:00:00
#SBATCH -p msismall,msibigmem

set -ue
set -o pipefail

# Sample and reference variables
ploidy_number=2
reference_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta
species=Calbicans
ref=SC5314-A21

mkdir -p genotyped_vcf
mkdir -p filtered_vcf

# Load modules
module load gatk/4.4.0
module load htslib/1.9
module load bcftools/1.10.2

# Joint genotyping
gatk --java-options "-Xmx29g" GenotypeGVCFs \
 -R "${reference_fasta}" \
 -V gendb://db/"${species}_${ref}" \
 -ploidy $((ploidy_number)) \
 -O genotyped_vcf/"${species}"_"${ref}"_haplotypecaller_genotype.g.vcf

# Subset to SNPs-only callset with SelectVariants
gatk --java-options "-Xmx29g"  SelectVariants \
    -V genotyped_vcf/"${species}"_"${ref}"_haplotypecaller_genotype.g.vcf \
    -select-type SNP \
    -O genotyped_vcf/"${species}"_"${ref}"_snps.vcf

# Subset to indels-only callset with SelectVariants
gatk --java-options "-Xmx29g" SelectVariants \
    -V genotyped_vcf/"${species}"_"${ref}"_haplotypecaller_genotype.g.vcf \
    -select-type INDEL \
    -O genotyped_vcf/"${species}"_"${ref}"_indels.vcf

# Filter SNPS
gatk --java-options "-Xmx29g" VariantFiltration \
-R "${reference_fasta}" \
-V genotyped_vcf/"${species}"_"${ref}"_snps.vcf \
-O filtered_vcf/"${species}"_"${ref}"_snps_filtered.vcf \
--filter-name "QD2" \
--filter-expression "QD < 2.00" \
--filter-name "QUAL30" \
--filter-expression "QUAL < 30.0" \
--filter-name "SOR3" \
--filter-expression "SOR > 3.0" \
--filter-name "FS60" \
--filter-expression "FS > 60.0" \
--filter-name "MQ40" \
--filter-expression "MQ < 40.0" \
--filter-name "MQRankSum-12.5" \
--filter-expression "MQRankSum < -12.5" \
--filter-name "ReadPosRankSum-8" \
--filter-expression "ReadPosRankSum < -8.0" \

# Output only passing calls
gatk SelectVariants \
-R "${reference_fasta}" \
-V filtered_vcf/"${species}"_"${ref}"_snps_filtered.vcf \
-O filtered_vcf/"${species}"_"${ref}"_snps_select.vcf \
-select "vc.isNotFiltered()"

# Filter indels
gatk VariantFiltration \
-R "${reference_fasta}" \
-V genotyped_vcf/"${species}"_"${ref}"_indels.vcf \
-O filtered_vcf/"${species}"_"${ref}"_indels_filtered.vcf \
--filter-name "QD2" \
--filter-expression "QD < 2.0" \
--filter-name "QUAL30" \
--filter-expression "QUAL < 30.0" \
--filter-name "FS200" \
--filter-expression "FS > 200.0" \
--filter-name "ReadPosRankSum-20" \
--filter-expression "ReadPosRankSum < -20.0" \

# Output only passing calls
gatk SelectVariants \
-R "${reference_fasta}" \
-V filtered_vcf/"${species}"_"${ref}"_indels_filtered.vcf \
-O filtered_vcf/"${species}"_"${ref}"_indels_select.vcf \
-select "vc.isNotFiltered()"

# Merge final vars
gatk MergeVcfs \
    --INPUT filtered_vcf/"${species}"_"${ref}"_snps_select.vcf \
    --INPUT filtered_vcf/"${species}"_"${ref}"_indels_select.vcf \
    --OUTPUT "${species}"_"${ref}"_filtered_merged.vcf

bgzip "${species}"_"${ref}"_filtered_merged.vcf
tabix "${species}"_"${ref}"_filtered_merged.vcf.gz

