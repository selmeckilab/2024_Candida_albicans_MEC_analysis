#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=2gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=2-100

set -ue
set -o pipefail

module use /home/selmecki/shared/software/modulefiles.local
module load dellyi

# Variables
line=${SLURM_ARRAY_TASK_ID}
sample_file=bam.files
ref_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes//SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta

# Get file and sample IDs
bam=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)
strain=$(basename "$bam" | cut -d "_" -f 1)

# Call SVs
delly call -q 20 -g "${ref_fasta}" -o "${strain}".vcf "${bam}"
