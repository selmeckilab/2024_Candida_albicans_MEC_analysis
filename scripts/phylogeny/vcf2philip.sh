#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msismall,msibigmem,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

vcf="../../haplotype_caller/Calbicans_SC5314-A21_filtered_merged.vcf.gz"
min=310

module use /home/selmecki/shared/software/modulefiles.local
module load vcf2phylip

vcf2phylip.py \
    -m "${min}" \
    -i "${vcf}"
