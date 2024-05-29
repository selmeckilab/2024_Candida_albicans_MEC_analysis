#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=12gb
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=12:00:00
#SBATCH -p msismall,msilarge

set -ue
set -o pipefail

species=Calbicans
ref=SC5314-A21
sample_list=Calbicans_envsample_map

# Load modules
module load gatk/4.4.0

gatk --java-options "-Xmx9g" GenomicsDBImport \
  --genomicsdb-update-workspace-path db/"${species}"_${ref} \
  --sample-name-map "${sample_list}"
