#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=5gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=2-303

set -ue
set -o pipefail

line=${SLURM_ARRAY_TASK_ID}
sample_file=bams_to_subsample.txt

# Read sample file line corresponding to array task ID and get variables
read1=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)
strain=$(basename "$read1" | cut -d "_" -f 1)

# Load modules for trimming and aligning
module use /home/selmecki/shared/software/modulefiles.local

module load bbmap
module load samtools/1.10
module load qualimap/20231012

reformat.sh in="${read1}" \
out=combined_data_sets/subsample_bam/"$strain"_subset.bam \
sampleseed=2 samplereadstarget=6000000

# Reindex
samtools index combined_data_sets/subsample_bam/"${strain}"_subset.bam

# Basic qc
samtools flagstat combined_data_sets/subsample_bam/"${strain}"_subset.bam \
   >  combined_data_sets/logs/"${strain}".log

unset DISPLAY # Qualimap won't work on cluster without this
qualimap bamqc -bam combined_data_sets/subsample_bam/"${strain}"_subset.bam \
    -outdir combined_data_sets/subsample_bam/"${strain}" --java-mem-size=4G
