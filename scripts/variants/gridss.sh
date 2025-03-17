#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=36gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=48:00:00
#SBATCH -p msismall,msibigmem
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

ref=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta
vcf=Calbicans_MEC_SV.vcf
bam=$(find ../combined_data_sets/subsample_bam/ -name "*.bam" | sort)

# Load modules
module use /home/selmecki/shared/software/modulefiles.local

module load gridss

gridss \
    --reference "${ref}" \
    --output "${vcf}" \
    "${bam}"
