# Variant calling with haplotype caller

One of the benefits of using GATK haplotype caller is that it allows joint
genotyping of samples added at different times, through the use of a database.
Joint genotyping for a population of samples is best practice - it's not good
to separately call variants on samples individually and then combine them.

## Workflow:

0. If samples have widely varying sequencind depth, consider subsampling bams for
more consistency (see "array_subsample.sh").

1. Variants are called on individual samples. See "haplotypecaller_1.sh" which is
a slurm array.

2. A sample map is generated, see "samplemap.sh". You can also just manually make
this, it's a tabbed file with sample name and path to individual vcf

3. Samples are added to the database. To create a new database, see
"genomicsdb_new.sh", and to add samples to an existing database, see
"genomicsdb_update.sh". These scripts use the "intervals.list" file which in our
case is just the chromosome names and full length of chromosomes from 1-NNN.

4. Joint genotyping of all samples is performed. See
"joint_genotyping_filtering.sh", which uses hard filtering parameters as
recommended by GATK for non-model organisms.

All of the scripts can be run sequentially, see "depend_haplotype.sh".

For further filtering, annotation with gene names (if available for your
species) and further filtering using bcftools, see "annotate_filter.sh".
