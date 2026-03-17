# MAVE Intra-Residue Correlation Analysis

This repository contains scripts used for the analyses in:

**“Determining the intra-residue correlation of missense variant impact using MAVE scores: implications for the ACMG/AMP PM5 criterion for DNA variant classification” (Dai et al., 2026)**

## Contents

- `inputs/`: Example input files containing missense variants derived from MAVE datasets.  
- `annotations/`: Additional files used by scripts. **Large files are not included** in the repository and must be downloaded separately:  
  - [`AlphaMissense_hg38.tsv`](https://github.com/google-deepmind/alphamissense): AlphaMissense scores  
  - [`HUMAN_9606_idmapping.dat`](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz): UniProt ID mapping  
  - [`variant_summary.txt`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz): ClinVar variant summary  

- `MAVE_primary_analysis.R`: Performs residue-level analysis of missense variants using MAVE scores, calculating LOF proportions, likelihood ratios, evidence strength, and generating summary plots.  

- `MAVE_secondary_analysis.R`: Performs variant-level annotation and integration of genomic coordinates, ClinVar, in silico predictors, and AlphaMissense scores for MAVE missense variants.  

- `MAVE_secondary_analysis_results.R`: Analyzes MAVE variant data to evaluate Grantham distances, compare functional classifications to in silico predictors and ClinVar annotations, and generate concordance and ROC visualizations.  

- `new_ACMG_rule.R`: Calculates a new ACMG-like score for MAVE variants by integrating ClinVar annotations, nucleotide/amino acid differences, Grantham distances, and in silico predictor scores, and outputs annotated MAVE and ClinVar lookup tables.  

- `modified_ACMG_rule.R`: Generates a modified ACMG-like score for MAVE variants by integrating LOF/FUNC variant concordance, amino acid differences with Grantham distances, and REVEL scores, producing both a filtered MAVE lookup table and a scored MAVE dataset.  

- `ACMG_rule_results.R`: Generates summary visualizations for ACMG-like scoring of missense variants, including stacked bar plots of total ACMG and PM5 scores, boxplots of predictor-specific VEP scores, and faceted heatmaps comparing grouped predictor scores with PM5 scoring, for both ClinVar-filtered and MAVE-derived datasets.  

- `VEP_discordant_analysis.R`: Computes per-gene and domain-stratified AUCs for VEP tools on MAVE missense variants, generating heatmaps and bar plots to visualize predictive performance and discordance.

## Usage

### 1. Set up your R environment
Install required packages (if not already installed):

```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "stringr", "httr", "jsonlite", "purrr"))
```

### 2. Prepare input files

Input file containing all missense variants from MAVE datasets are located in `inputs/`.  

Download the large annotation files into `annotations/`:

- [`AlphaMissense_hg38.tsv`](https://github.com/google-deepmind/alphamissense)  
- [`HUMAN_9606_idmapping.dat`](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz)  
- [`variant_summary.txt`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)

### 3. Run the scripts

Execute the analysis scripts in order, depending on your goal:

- `MAVE_primary_analysis.R` — Residue-level analysis
- `MAVE_secondary_analysis.R` — Variant-level annotation
- `MAVE_secondary_analysis_results.R` — Evaluate Grantham distances, ROC, concordance
- `new_ACMG_rule.R` — Compute ACMG-like scores
- `modified_ACMG_rule.R` — Compute modified ACMG scores
- `ACMG_rule_results.R` — Generate summary plots
- `VEP_discordant_analysis.R` — Compute per-gene/domain AUCs and heatmaps