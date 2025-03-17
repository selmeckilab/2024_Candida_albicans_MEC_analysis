#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=12gb
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=6:00:00
#SBATCH -p msismall,msilarge
#SBATCH --array=1

set -ue
set -o pipefail

# Sample and reference variables
sample_file=bam.files
line=${SLURM_ARRAY_TASK_ID}
ploidy_number=2
reference_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta

# Get strain ID froms sample file line equal to array task ID
bam=$(awk -v val="${line}" 'NR == val { print $1}' "${sample_file}")
strain=$(basename $bam | cut -d "_" -f 1)

mkdir -p gvcf

# Load modules
module load python
module load gatk

# HaplotypeCaller and Genotype
gatk --java-options "-Xmx8g -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R "${reference_fasta}" \
 -I "${bam}" \
 -O gvcf/"${strain}"_haplotypecaller.g.vcf \
 -ploidy $((ploidy_number)) \
 -ERC GVCF
