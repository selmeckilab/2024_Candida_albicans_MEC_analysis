# MEC Candidemia Study *C. albicans* analysis

### 2024-06-28-24
#### Global alignment of MEC fastq files to A21 ref
Finalizing CNV data for dissertation chapter and reviewing potential gene
deletions. Aligned reads to ref using bbmap to see how much soft clipping by
bwa is affecting results for some poorly-aligning regions. See
array_bbduk_bbmap.sh. 

### 2024-06-28
#### MTL locus checks
Aligned isolates with potential MTL homozygosity based on coverage depth (not
equal to 1) to A22 ref genome to see alignments at a and alpha. See bam files in 
align_variants_SC5314_A21_ref/A22_ref_bam/.

### 2024-06-27
#### CNV breakpoint analysis
Searched for small-scale breakpoints based on 500-bp sliding windows using
samtools mpileup files from genome plotting and subsetting to intersect genes
found in A21 gff. Visually inspected potential deletions and amplifications in
IGV. Calculated CNV distance to nearest repeat region using the A21 repeat GFF
created by Travis Stratton (made a copy of the GFF and subset to remove
duplicates). See cnv_repeat_proximit.sh.

### 2024-05-19
#### Expanded data set - Ropars, Hirakawa, MEC
Ran RAxML on 299 samples with 1000 bootstraps and no missing genotypes using
filtered freebayes VCF (to keep topology most similar to the MEC tree built
using freebayes).

### 2024-05-07
#### Expanded data set
Ran RAxML on 307 samples with 100 bootstraps and default vcf2phylip params (min
4 samples per site) using filtered haplotype_caller VCF.

### 2024-04-30
#### Expanded data set - Davis lab env samples
Haplotype caller workflow (initial var calls and genomicsdb update). See
scripts/haplotype_caller_scripts.

### 2024-04-26
#### Expanded data set - Davis lab env samples
Preliminary trimming, alignment and QC of DAD32,35,36.

### 2024-03-25
#### SV annotation
Annotated all MEC SVs for simple events with StructuralVariantAnnotation (see
sv_annotation.R). Split multisample VCF into individual sample VCFs and subset
regions to exclude repetitive regions (subtelomeric, MRS, LTR) and
retrotransposons (as denoted in A21 gff) and then visualized results with ggbio
and karyoploter (see array_split_vcf.sh and ind_SVs.R).

### 2024-03-21
#### SV annotation and filtering
Ran gridss for joint SV calling. Will filter for all passing quality, min 0.15
AF for at least one sample. To annotate and filter further will use
StructuralVariantAnnotation in Bioconductor and maybe visualize in karyoploteR.
#### Samplot for quick SV visualization
Remembering why I did not pursue this cl tool further - it does not accept
non-standard contig names. Can only use after running reheader and changing chr
names to things like "chr1" etc. With only BND annotation in gridss vcf also
cannot annotate image created anyways.

### 2024-03-05
#### Revisit VCF filtering
Running a "parameter sweep" to understand impacts of filtering on INFO fields
SAR, SAF, RPR, RPL, variant type (complex) and VAF for het vs homozygous
positions.

### 2024-02-16
#### LOH visualization
Using saved genome plot excel sheets to generate a combined LOH heatmap for all
MEC isolates. R package pheatmap plots the heatmap well, but converting to
matrix and transposing samples vs sites is not working well. Ggplot geom_tile
is plotting the x-axis incorrectly whenever more than one chrom is plotted,
though it can work as facets.

### 2024-02-15
#### SV calling
Finished perSVade for all MEC isolates. Will use sv_functions.py as an
interactive module called by the singularity container to try integrating
samples as described in the github repo FAQs.

### 2024-02-10
#### Continue SV calling
perSVade memory usage is reasonable but it creates thousands of intermediate
files. Multiple samples failed (run as an array) with write errors. Re-running
in smaller batches to decrease file count.
Singularity container setup seems to be having some issues with file
permissions in the shared directory, will run in home dir.

### 2024-02-09
#### MEC isolates - SV calling
Running perSVade v1.02.7 as singularity container to call CNVs
(coverage-based), calculate gene coverage, call SVs, integrate CNVs and SVs,
and annotate. Also running optimize_parameters module for each isolate using
homologous regions.

### 2024-01-30
#### Combined data set - GATK variant calling
Divide the 304 vcfs into groups of 20-30, use GATK v4.4.0, and incrementally
build the db. Each group takes ~6 hours to add.

### 2024-01-24
#### SRA upload
All 100 MEC isolates uploaded to NCBI SRA, BioProject PRJNA1068683, via
aspera.

### 2024-01-23
#### Additional external data files
Added 5 isolates sequenced by Davis lab: 2 animal fecal isolates (DDS005,
DDS006) and 3 commensal isolates (DDS001, DDS002, DDS003). Ran standard
trimming, alignment and QC. Sample with highest coverage has 6.5 million reads
so will not do any downsampling but will simply add to all other isolates and
run variant calling with Haplotypecaller.

### 2024-01-15
#### Combined data set variant calling
Chromosomes 1,2,3 and R timed out (48 hours). Restarting with 96 hours, will
further divide chrs if this job fails.

### 2024-01-13
#### External data set processing
Performed standard trimming, alignment to SC5314 A21 and qc for Ropars strains.
Downsampled all isolates to max of 6million reads (~50x coverage) (see
"array_subsample.sh"). Re-ran Qualimap and MultiQC for review. Removed 3
strains for coverage less than 20x (SRR6669970, SRR6669880, SRR6669899). Lowest
coverage is an MEC isolate with ~27x coverage.
Started Freebayes for all isolates, split into array by chromosomes.

### 2024-01-12
#### Add external *C. albicans* data for phylogenetic tree building
Hirakawa data set:
Paired-end WGS data is present in External_data directory, fastq trimming,
alignment to SC5314 A21  and QC has been performed. BAM files copied back from 
second-tier storage.

Ropars data set:
Only some of the WGS data for 2018 Ropars et al is present in the External_data
directory and no indication of how or when it was downloaded. The Bioproject ID
is PRJNA432884. Downloaded list of SRR IDs and metadata for all available
samples (NCBI website says 182 samples selected but only 181 accession numbers
and metadata lines print). Downloaded all paired-end reads, split and
compressed with SRA toolkit.

### 2024-01-09
#### RAxML tree repeat
1000x bootstrapping timed out (24h).

### 2024-01-08
#### RAxML tree repeat - 1000 bootstrap support
Rerun with 1000x boostrapping - same number used by Ropars 2018 paper. They
also took a no-outgroup, midpoint-rooting approach.

### 2024-01-04
#### RAxML tree generation
Run RAxML with GTRGAMMA model and fast bootstrapping (100x), no outgroup
included. See phylogeny.sh. Visualize output in R with ggtree. Note for install:
version conflict (R version vs ggplot2 version), downgraded ggplot package to 
work with versions available on MSI. See 'tree_viz.R' in 'pop_structure' repo.

### 2024-01-02
#### VCF filtering and phylip file generation
Concatentate new freebayes chr vcfs, remove SNPs called in AMS2401, remove
AMS2401, filter on quality, subset to exclude known repetitive regions,
annotate for final MEC population list. Use vcf2phylip.py to create phylip file
(default parameters) for RAxML phylogenetic tree building.

### 2023-12-04
#### AMR SNP upload to redcap
See 'redcap_variant_upload.R' in 'clinical_isolate_redcap' repo.

### 2023-11-30
#### AMR gene annotation update
Re-filtered again to same AMR genes. Subset and exported csv for redcap upload.

### 2023-11-29
#### Freebayes parameter updates
Rerun variant calling with no min allele frequency and drop the minimum
supporting reads to 5. Filtering for depth and AF will be done post-calling
with BCFtools.
#### Genome plot saving updates
Resize, set device to png with dpi of 300 for now. Need to figure out best
saving settings for publication-quality/Illustrator read figures. 

### 2023-09-25
#### AMR gene annotation and filtering
Re-filtered variant file with SC5314 included, and subset to regions of
interest (Calbicans_amr_genes.bed) before removing all variants found in
SC5314. Will subset out synonymous variants before uploading to redcap.

### 2023-09-21
#### AMR gene annotation and filtering
Subset filtered and annotated VCF to genes of interest
(Calbicans_amr_genes.bed.gz). Removed indels, split multiallelic SNPs,
calculated per-sample VAF `bcftools +fill-tags -- -t FORMAT/VAF`. Reviewed in
IGV, and added SC5314 BAM. Probably (?) moderate and high-impact SNPs are not
shared by SC5314 but will re-call the population with it as a variant, and then
just use bcftools to remove all variants shared by SC5314 before
continuing with annotation and upload to redcap.

### 2023-04-25
#### Variant-based clustering
Re-did SNP-based MCA without factoextra wrapper for plotting - plotted 1st 2
dimensions of eigen coords with ggplot. Script "calbicans_mca.R".

### 2023-04-20 - 2023-05-01
#### Genome-wide copy number, SNP heterozygosity and allelic frequency plotting
Ran new workflow ("calbicans_gc_correct.sh", "calbicans_ymap.sh"):

1. gc correction with deeptools
2. read depth per base with samtools depth
3. SNPs per base with samtools mpileup and berman lab python script for counts per site
4. horizontal genome plot with "allele_ratio_genome_vis.R"

### 2022-12-12
#### Variant-based clustering
More troubleshooting of clustering and dendrograms.
Errors in patient list - reviewed redcap series and rewrote sample sheet.
MCA (default params) first 2 dimensions only explain ~10 and ~8% of variation.
K-means clustering doesn't show an obvious number of clustering when plotting
BIC values.

### 2022-12-11
#### Variant-based clustering
Troubleshooting clustering (R packages poppr and FactoMineR).

### 2022-12-09
#### Variant-based clustering
Generated tab-delimited genotype table (columns: chr, pos, per sample GT; rows:
variants), removing fixed SNPs or missing genotypes (bcftools).
Created csv of sample IDs and patient IDs (for stratification by patient).
R script refactoring for clustering (MCA and k-means) and dendrogram creation
(no outgroups).

### 2022-12-05
#### Variant calling, filtering annotating to SC5314 A21 assembly
Concatenated chr vcf files, filtered (see script "albicans_filter_annotate.sh"
for params). Annotated with snpeff (custom database, SC5314_s02m09r08, yeast
alternative codon table).

### 2022-12-02
#### Variant calling, filtering annotating to SC5314 A21 assembly
Repeated variant calling with MEC005 added.
Same parameters, previous vcf files deleted.

### 2022-11-17
#### Variant calling, filtering annotating to SC5314 A21 assembly
Population level variant calling with freebayes (script
"array_alb_freebayes.sh").
Array jobs split on chromosomes, calling all 100 samples per.
Freebayes parameters: -C 10, -F 0.4, -p 2

### 2022-11-21
#### Data pre-processing, alignment to SC5314 A21 assembly, QC
MEC005 is *C. albicans* (labeled as *C. dubliniensis*). Updated Calbicans
sample sheet and sequence paths, trimmed, aligned to reference.

### 2022-11-13
#### Data pre-processing, alignment to SC5314 A21 assembly, QC
Basic QC (trimmed reads, alignment, classifier) as batch job ("basic_qc.sh").

### 2022-11-11
#### Data pre-processing, alignment to SC5314 A21 assembly, QC
Re-did adapter and quality trimming with BBDuk, aligned newly trimmed reads to
SC5314 assembly 21 reference genome (script"array_bbduk_align_sort.sh").

### 2022-07-15
#### Data pre-processing, alignment to SC5314 A21 assembly, QC
MEC008 (labeled as *C. albicans*) had ~75% mapped reads to SC5314 A21 reference
genome and was called as *C. dubliensis* by centrifuge (~90% of reads). Raw
sequencing data was moved to proper subdirectory, intermediate files removed
from this project directory, sample list and redcap information updated.

### 2022-07-12
#### Data pre-processing, alignment to SC5314 A21 assembly, QC
Reviewed multiqc for logs and centrifuge output. Re-running centrifuge for MEC
131,132,133,134 on the trimmed (1P and 2P) output from trimmomatic. Qualimap
shows short insert sizes and FastQC shows high nextera transposase adapters -
so the Centrifuge output with Ralstonia reads are probably
read-through/contamination from adapters.

### 2022-07-11
All MEC isolates have been sequenced by MiGS/SeqCenter - gDNA prep and
sequencing info in REDcap. Current fastq file paths are in albicans_samples.txt
Performed QC, trimming with trimmomatic, alignment, sorting/marking duplicates per standard
Selmecki workflows (script "array_trim_align_sort.sh"). MEC103_1 is purified
single colony subculture (MEC103-1 in REDcap) from original MEC103 mixed
culture (received as *C. albicans*) that was sequenced in November 2021. Ran
centrifuge (interactive, conda env) with ncbi-nt db.
