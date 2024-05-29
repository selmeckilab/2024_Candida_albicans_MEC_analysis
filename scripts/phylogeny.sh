#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=48:00:00
#SBATCH -p msilarge,msismall
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

module load raxml/8.2.11_pthread

raxml -f a -m GTRGAMMA -p 12345 -x 12345 -# 1000 -s Calbicans_MEC_bwa_filtered_annotated.min4.phy -T 16 -n T2
