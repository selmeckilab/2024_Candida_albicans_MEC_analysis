#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=20
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err

job1=$(sbatch --parsable haplotypecaller_1.sh)
job2=$(sbatch --parsable --dependency=afterok:${job1} samplemap.sh)
job3=$(sbatch --parsable --dependency=afterok:${job1}:${job2} genomicsdb_new.sh)
job4=$(sbatch --parsable --dependency=afterok:${job1}:${job2}:${job3} joint_genotyping_filtering.sh)

