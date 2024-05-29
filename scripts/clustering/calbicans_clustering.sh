#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=36:00:00
#SBATCH -p msismall
#SBATCH -o %j.out
#SBATCH -e %j.err
set -ue
set -o pipefail

module load R/4.1.0

Rscript --vanilla calbicans_mca.R
