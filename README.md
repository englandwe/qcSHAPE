# qcSHAPE
Preliminary QC analysis of SHAPE data

This script streamlines the process of running QC checks on icSHAPE, icLASER, or similar RT-stop-based RNA structure data.  Works on outputs from the icSHAPE pipeline (https://github.com/qczhang/icSHAPE); run this script in its output folder.  Parallelizable.

Outputs include:
- Per sample
  - Count & percentage of reads mapped
  - Count & percentage of RT stops per base
  - Count of reads mapped per sequence biotype
- Per experiment
  -   Correlations between samples (RPKM & RT stops)
  -   Reactivity by region (5'UTR, 3'UTR, ORF, start & stop codons)
  -   Enrichment scores by base
  -   Per-position enrichment scores for 18S & 28S rRNA
 
## Usage: 

qcSHAPE.py THREADS IDFILE GTF FASTA STOPS,RPKM,SCORE 

THREADS: number of threads to use

IDFILE: two-column tab-delimited file.  First column is sample ID; second column is group (i.e. treatment or control).

GTF: .gtf file of reference transcriptome

FASTA: .fasta file of reference genome

STOPS,RPKM,SCORE: minimum cutoffs for number of RT stops, RPKM, and enrichment score, as a comma-separated list


  
