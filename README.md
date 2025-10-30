# *Candida albicans* MEC isolate analysis

## Study isolates
There are 101 *C. albicans* MEC isolates analyzed for this project (see
"Calbicans_sequencing_paths.txt" for list of isolate IDs). All MEC isolates are
single colony subcultures from blood cultures collected from UMMC and affiliated 
hospitals. These were collected under IRB IDs STUDY00006473 and STUDY00021428.

## Metadata
Clinical metadata is available for study members in a Box spreadsheet and UMN 
data shelter. Isolate data is available for study members in a REDCap 
database.

## DNA extraction and Sequencing
DNA extraction was performed with the standard Selmecki lab phenol protocol
(bead beating setting may have varied) in multiple batches between 2021 and
2024.
Illumina sequencing for all isolates was performed at SeqCenter (formerly MiGS)
on a NextSeq 2000 (the 400 Mb sequencing package, paired end 151 bp reads).
Insert sizes range from ~100 bp to > 400bp and median coverage ranges from 27x
to 111x.

## Sequencing reads
All files have been submitted to SRA (NCBI BioProject PRJNA1068683).

## Trimming, alignment, QC
Samples were trimmed with BBDuk (BBMap v38.94) and aligned to the SC5314 A21
reference genome with bwa mem v0.7.17, default parameters, followed by sorting
and duplicate marking with samtools v1.10 (see "scripts/classify_aliqn_qc/array_bbduk_align_sort.sh").
Basic stats were generated with samtools flagstats. Trimmed reads and
alignments were QC'd with FASTQC v0.11.9  and Qualimap v2.2.2-dev, and results
compiled with MultiQC v1.13 (see "scripts/classify_align_qc/basic_qc.sh").
Isolates with low percentage of aligned reads were investigated by
classification of fastq files using Centrifuge v1.0.4 (see
"scripts/classify_align_qc/array_centrifuge_classifier.sh"), removed from this 
analysis and updated in the study database if found to be a different species 
(MEC008 removed, MEC005 added).

## Multilocus sequence typing
"Consensus sequences" of the 7 MLST loci were generated for all 100 MEC isolates
(see "scripts/mlst/"). The PubMLST database was queried via API to
assign allele IDs and sequence types and to upload to the study database. New 
alleles were submitted to PubMLST for assignment of
an ID.

## Genome visualization
A local version of YMAP's genome-wide LOH and CNV plotting was performed  using
Deeptools for GC bias correction and samtools for depth and SNP calling. See
"scripts/genome_vis/".

## Variant calling, annotation, filtering
Variants were called with freebayes v1.1.0, with parameters of -p 2 and -C 5,
using each chromosome as a separate region (see "scripts/variants/array_alb_freebayes.sh").
Chromosome concatenation, quality filtering and annotation were performed with
BCFtools v1.17 and SnpEff v5.0e (see "scripts/variants/concat_filter_annotate.sh").

## Clustering
Variant sites with biallelic SNPs and no missing genotypes were used to
generate a matrix for clustering ("concat_filter_annotate.sh"). Multiple
correspondence analysis was performed using the FactoMiner package v2.9 (see
"scripts/clustering/"). 

## Publicly available data
### Sequencing data
Sequencing data from BioProjects PRJNA193498 (WGS SRA runs with library name 
“Pond-NNN”) and PRJNA432884 (all available SRA runs) were downloaded from NCBI 
Sequence Read Archive using the SRA toolkit (v3.0.0). Trimming, alignment and QC
were performed as described above. Samples with < 20x mean coverage were 
discarded. 

### HaplotypeCaller variant calling, annotation, filtering
Joint genotyping is recommended for variant calling of multiple samples (i.e.,
not combining VCF files of samples that were called separately). To facilitate 
variant calling of publicly available sequencing data in addition to our study 
isolates, we used GATK's HaplotypeCaller workflow (see 
"scripts/haplotype_caller_scripts/"). Variants were called for all isolates with 
GATK v4.1.2. A GenomicsDB object was created and VCF data from all isolates was 
imported. Joint genotyping and filtering was performed following recommendations
from GATK. Subsequent isolates were added to the GenomicsDB object and joint
genotyping/filtering was repeated as needed.

### Phylogenetic tree construction, cluster determination and visualization
The filtered VCF file was subset to include only SNP sites with no missing data 
and converted to PHYLIP format using vcf2phylip v2.8 (see 
"scripts/phylogeny/vcf2philip.sh"). A maximum-likelihood tree was built using 
RAxML v8.2.11 with the GTR+GAMMA model of nucleotide substitution and rapid 
bootstrapping (100x) (see "scripts/phylogeny/phylogeny.sh"). No outgroup was
specified. To define clusters, we adapted the approach described by 
Schikora-Tamarit and Gabaldón (2024). We defined cluster nodes as internal nodes
with bootstrap support >95%, having subtending branch lengths above a minimum 
threshold, and having no child nodes that also met the definition for cluster 
nodes. We compared a range of minimum threshold values to identify a value 
(branch length > 0.03) that had the most agreement with existing MLST and 
WGS-based clusters (Odds et al., 2007; Ropars et al., 2018) and minimized the 
number of singleton isolates. Midpoint rooting of the RAxML bipartition file was
performed with the R package phangorn v2.11.1. Isolate metadata was added to 
create a treedata object which was visualized using the R packages ggtree v0.4.6
and tidytree v3.10.1 (see "scripts/figures/Calbicans_gatk_expanded_data_tree.R").



