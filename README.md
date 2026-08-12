# Don-Hare Analysis Pipelines

The repository contains two independent analysis lines for hare genomic data.
The population-genetic line produces filtered VCF files, PCA, and ADMIXTURE
results. The BUSCO line extracts BUSCO gene sequences, filters and aligns them,
and builds gene and species trees.

The samples `1k`, `3k`, `5k`, and `5kS8` are treated as a separate comparison
group. Both lines produce results first without these samples and then with all
samples.

## Line 1: VCF and Population Analyses

Run the scripts in this order:

### `scr_11_bam_to_vcf.sh`

Calls variants from Paleomix BAM files, restricts variant calling to regions
listed in `chr.txt`, merges per-sample VCF files, and creates two filtered VCF
datasets: one without `1k`, `3k`, `5k`, `5kS8`, and one with all samples.

Filters: `QUAL >= 20` and `MIN(FMT/DP) >= 3`.

### `scr_12_pca_admixture.sh`

Runs PLINK PCA and ADMIXTURE for both VCF datasets. PLINK applies `--geno 0.2`
and calculates four principal components. ADMIXTURE tests `K=2..15` and saves
the K with the lowest cross-validation error.

### `scr_13_population_plots.R`

Creates PCA scatter plots and ADMIXTURE ancestry-proportion plots for both
datasets using the outputs from `scr_12_pca_admixture.sh`.

### ABBA-BABA and fbranch

These analyses are not included yet. They require defined groups `P1`, `P2`,
`P3`, `P4`, an outgroup, and a selected implementation such as Dsuite.

## Line 2: BUSCO Gene Analysis and Trees

Run the scripts in this order:

### `scr_21_bam_to_busco.sh`

Extracts BUSCO gene loci from BAM files and prepares two merged FASTA datasets:
one without `1k`, `3k`, `5k`, `5kS8`, and one with all samples.

Filters and settings: `MIN_COVERAGE=2`, discard sequences shorter than
`MIN_SEQ_LENGTH=10`, and use BUSCO coordinates from `BED_FILE`.

### `scr_22_mafft.sh`

Aligns both merged FASTA datasets with MAFFT using nine parallel threads.

### `scr_23_filter_by_gaps.sh`

Filters aligned BUSCO genes independently for both datasets and writes lists of
genes that pass the filters.

Filters: maximum sequence length `5000`, minimum ATGC count per sample `30`,
and maximum gap percentage `30%`. The dataset without the four excluded
samples has no required-sample check. In the all-sample dataset, `1k`, `3k`,
`5k`, and `5kS8` are all required. If any required sample is absent from a
gene FASTA file, that gene is excluded from the gene list and is not analysed
further.

### `scr_24_iqtree_astral.sh`

Builds IQ-TREE gene trees with the `GTR+G` model for both filtered datasets,
then runs ASTRAL separately to produce one species tree without the excluded
samples and one species tree with all samples.

## Other Files

### `scr_98_sex.sh`

An independent sex-analysis script based on X-chromosome and autosomal
coverage. BAM filtering and indexing commands are currently commented out and
must be enabled or revised before use.

### `scr_99_old_ngs_analysis`

Archived older analysis material; it is not part of the two main workflows.
