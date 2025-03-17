#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=64gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=24:00:00
#SBATCH -p msismall,msibigmem,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

seed=310
bootstrap=100
in_file=Calbicans_SC5314-A21_filtered_merged.min310.phy
threads=126
out_file=Calbicans_310

module load raxml/8.2.11_pthread

raxml \
    -f a \
    -m GTRGAMMA \
    -p "${seed}" \
    -x "${seed}" \
    -# "${bootstrap}" \
    -s "${in_file}" \
    -T "${threads}" \
    -n "${out_file}"
