#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 15
#SBATCH -p msismall,msilarge

set -ue
set -o pipefail

species=Calbicans

# Create sample map
find gvcf -name '*.g.vcf' -type f -printf '%p\n' \
 | awk -F'[/_]' 'BEGIN{OFS="\t"} {print $2, $0}'> "${species}_envsample_map"
