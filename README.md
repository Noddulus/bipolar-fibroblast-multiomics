# Fibroblasts as an in vitro model of circadian genetic and genomic studies
## Project Overview
Bipolar disorder (BD) is a highly heritable psychiatric disorder characterized by shifts in mood and systemic circadian abnormalities. This repository contains the complete computational pipeline for investigating the genetic architecture of circadian rhythms in an in vitro fibroblast model, assessing its contribution to the polygenic risk of BD.

Using primary cell lines from 6 healthy individuals, we collected temporal genomic features over a 48-hour period, integrating transcriptomic (RNA-seq) and open chromatin (ATAC-seq) data. The analytical workflow characterizes the biological pathways activated in this model, evaluating their relevance within the known genetic architecture of BD and related psychiatric traits through stratified linkage disequilibrium score regression (sLDSR).

## Publication:
Francia, M. et al. (2024). Fibroblasts as an in vitro model of circadian genetic and genomic studies. Mammalian Genome. https://doi.org/10.1007/s00335-024-10050-7

## Software & Dependencies
<ins>Command Line Tools (Linux/HPC Environment)</ins>
* Quality Control & Trimming: FastQC, TrimGalore, Cutadapt

* Alignment & Processing: STAR (GRCh38), Bowtie2, Samtools, featureCounts

* Peak Calling & Manipulation: MACS2, bedtools

* Motif Enrichment: HOMER

* Genetic Architecture: LDSC (ldsc.py)

<ins>R Packages (Version 4.1.1)</ins>
* Differential Expression & Normalization: limma, Glimma, edgeR

* Network Analysis: WGCNA

* Circadian Detection: MetaCycle (ARSER, JTK, LS), RAIN

* Temporal Clustering & Modeling: splines, dtwclust

* Annotation: Homo.sapiens

## Pipeline Execution
<ins>Step 1: RNA-seq Processing and Analysis</ins>
* Quality Control & Alignment: Raw reads are evaluated with FastQC and trimmed using TrimGalore/Cutadapt. Reads are aligned to the human reference genome (GRCh38) using STAR.

* Quantification: Samtools indexes the aligned files, and featureCounts associates read counts with NCBI GRCh38 gene annotations.

* Normalization & Filtering: Genes with low counts are filtered out and normalized (CPM) using limma, Glimma, and edgeR.

<ins>Network & Circadian Analysis:</ins>

* WGCNA classifies genes with similar temporal patterns (power = 12, signed network).

* MetaCycle and RAIN identify circadian expression patterns.

* Cubic splines (5 degrees of freedom) are fitted using the splines package to model temporal dynamics, followed by FDR correction and time-series clustering with dtwclust.

<ins>Step 2: ATAC-seq Processing and Analysis</ins>
* Alignment & Processing: Following the ENCODE pipeline, reads are trimmed and aligned to GRCh38 using Bowtie2 (2 kb insert size, up to 4 alignments). Samtools removes blacklisted regions and PCR duplicates.

* Peak Calling: Open chromatin regions (OCRs) are identified using MACS2 `-g hs -q 0.01 --nomodel --shift -100 --extsize 200`. Overlapping peaks across timepoints/subjects are merged into a consensus bed file using bedtools.

* Quantification & Network Analysis: featureCounts assigns reads to consensus regions (normalized by RPM). WGCNA clusters peaks with similar temporal accessibility patterns.

* Motif Enrichment: HOMER identifies enriched transcription factor motifs (e.g., bHLH, glucocorticoid receptor target motifs) within dynamically accessible peak modules.

<ins>Stratified Linkage Disequilibrium Score Regression (sLDSR)</ins>
To assess if the dynamically accessible chromatin regions contribute to BD heritability:

* ATAC-seq modules (increasing and decreasing accessibility) are extended by 1 kb and 10 kb genomic windows.

* ldsc.py computes partitioned SNP-based heritability against GWAS summary statistics.

* Models use the 1000 Genomes Project phase 3 reference panel, accounting for the full baseline model.

## Data Availability
Due to patient privacy and ethical restrictions, raw .fastq sequencing files and identifying patient information are not publicly hosted in this repository.
Summary statistics, processed read count matrices, and dynamic ATAC-seq peak files utilized for clustering and sLDSR are available in the supplementary materials of the published manuscript, and in the following repositories: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263711 & https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263713
