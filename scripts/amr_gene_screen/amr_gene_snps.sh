#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 20
#SBATCH -p msilarge,msismall
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

regions=Calbicans_amr_genes.bed.gz
input=Calbicans_MEC_bwa_filtered_annotated.vcf.gz
output=Calbicans_MEC_bwa_filtered_annotated_amr.vcf

module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load htslib

bcftools view -R "${regions}" \
    "${input}" > "${output}"

bgzip "${output}"
tabix "${output}".gz
