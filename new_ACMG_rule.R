# Load libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(httr)
library(jsonlite)
library(data.table)
library(Rsamtools)

# Read in MAVE data
maveData = fread("secondary_results/data.csv")
geneList = unique(maveData$gene_only)

# Make sure maveData is a data.table
maveData = as.data.table(maveData)

# Create new column
maveData$new_ACMG_score = 0

# Read in ClinVar annotation files
varSummary = fread("annotations/variant_summary.txt")

# Filter by genome build and genes in geneList
varSummary_filtered = varSummary[varSummary$Assembly == "GRCh38" & varSummary$GeneSymbol %in% geneList & varSummary$OriginSimple == "germline",] 

# Split Name column 
varSummary_filtered = varSummary_filtered %>% 
  separate(col = Name,
           into = c("Name", "protein_variant"),
           sep = " ", remove = TRUE,
           convert = TRUE) # split Name column

# Format protein_variant column
varSummary_filtered = varSummary_filtered %>%
  mutate(protein_variant = gsub("[()]", "", protein_variant))

# Create p_format column 
varSummary_filtered$p_format = NA
varSummary_filtered$p_format = ifelse(
  is.na(varSummary_filtered$protein_variant),
  NA,  # or "" if you prefer empty string
  paste0(varSummary_filtered$GeneSymbol, ":", varSummary_filtered$protein_variant)
)

# Define terms to exclude
exclude_terms = c(
  "Conflicting classifications of pathogenicity",
  "Uncertain significance",
  "Conflicting classifications of pathogenicity; risk factor",
  "Conflicting classifications of pathogenicity; other; risk factor",
  "not provided",
  "Uncertain significance/Uncertain risk allele",
  "Uncertain risk allele",
  "Uncertain significance; risk factor"
)

# Select only needed columns and exclude unwanted rows
clinvarTable = varSummary_filtered %>%
  dplyr::select(
    Type,
    Name,
    protein_variant,
    GeneSymbol,
    ClinicalSignificance,
    Chromosome,
    Start,
    Stop,
    ReferenceAlleleVCF,
    AlternateAlleleVCF,
    p_format
  ) %>%
  dplyr::filter(
    !is.na(protein_variant),
    !ClinicalSignificance %in% exclude_terms
  )

# Extract amino acid position from protein_variant
clinvarTable = clinvarTable %>%
  mutate(
    aa_pos = as.integer(str_extract(protein_variant, "\\d+")),
    ref_AA = str_extract(protein_variant, "(?<=p\\.)[A-Z][a-z]{2}(?=\\d)"),
    alt_AA = str_extract(protein_variant, "(?<=\\d)[A-Z][a-z]{2}$")
  )

# Create a lookup table for 3-letter to 1-letter amino acids
aa_lookup = c(
  Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
  Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
  Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
  Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
  Ter = "*", Stop = "*" # for stop codons
)

# Convert ref_AA and alt_AA to 1-letter codes
clinvarTable = clinvarTable %>%
  mutate(
    ref_AA_1 = aa_lookup[ref_AA],
    alt_AA_1 = aa_lookup[alt_AA]
  )

# Exclude rows where alt_AA_1 is nonsense or other (* or NA)
clinvarTable = clinvarTable %>%
  filter(!is.na(alt_AA_1), alt_AA_1 != "*", !is.na(ref_AA_1), ref_AA_1 != "*")

# Annotate with Grantham distance
# Load in Grantham matrix
grantham_data = read.delim("annotations/grantham.matrix", header = T, row.names = 1)
grantham_matrix = as.matrix(grantham_data)

# Annotate
clinvarTable = clinvarTable %>%
  rowwise() %>%
  mutate(
    grantham_distance = ifelse(
      !is.na(ref_AA_1) & !is.na(alt_AA_1) &
        ref_AA_1 %in% rownames(grantham_matrix) &
        alt_AA_1 %in% colnames(grantham_matrix),
      grantham_matrix[ref_AA_1, alt_AA_1],
      NA
    )
  ) %>%
  ungroup()

# Get ClinVar P/LP table
clinvar_PLP = clinvarTable %>%
  filter(ClinicalSignificance %in% c("Pathogenic", 
                                     "Pathogenic/Likely pathogenic", 
                                     "Likely pathogenic", 
                                     "Likely pathogenic/Likely risk allele", 
                                     "Pathogenic/Likely pathogenic/Likely risk allele")) %>%
  mutate(ClinicalSignificance = ifelse(grepl("Pathogenic", ClinicalSignificance), "P", "LP"))

# Get ClinVar B/LB table
clinvar_BLB = clinvarTable %>%
  filter(ClinicalSignificance %in% c("Benign", 
                                     "Benign/Likely benign",
                                     "Likely benign", 
                                     "Benign; association")) %>%
  mutate(ClinicalSignificance = ifelse(grepl("Benign", ClinicalSignificance), "B", "LB"))

### CONDITION : DIFF NT, SAME AA CHANGE P/LP AND B/LB ###
# Convert to data.table
DT  = as.data.table(maveData)
PLP = as.data.table(clinvar_PLP)
BLB = as.data.table(clinvar_BLB)

# Keep row index for scoring
DT[, i_row := .I]

# -------- P/LP --------
# Join PLP to DT on p_format (amino acid match)
PLP_join = PLP[DT, on = .(p_format), nomatch = 0, allow.cartesian = TRUE]

PLP_join[, dt_pformat := p_format]  # save the p_format from DT

# Calculate nucleotide differences with explicit p_format check
PLP_join[, diff_nt := (p_format == dt_pformat) & (
  (Chromosome != chr) |
    (Start      != start) |
    (Stop       != stop) |
    (ReferenceAlleleVCF != ref_plus) |
    (AlternateAlleleVCF != alt_plus)
)]


# Keep only rows with same AA but different NT
PLP_diff = PLP_join[diff_nt == TRUE]

# Compute score per MAVE row
# Base: P = 4, LP = 2
# Bonus: +2 if >1 matching P/LP entry
PLP_scores = PLP_diff[, .(
  score = fcase(
    any(ClinicalSignificance == "P") & .N > 1, 6L,  # 4 + 2
    any(ClinicalSignificance == "P"),               4L,
    any(ClinicalSignificance == "LP") & .N > 1, 4L, # 2 + 2
    any(ClinicalSignificance == "LP"),              2L,
    default = 0L
  )
), by = i_row]

# -------- B/LB --------
# Join BLB to DT on p_format
BLB_join = BLB[DT, on = .(p_format), nomatch = 0, allow.cartesian = TRUE]

# Immediately save DT's p_format in a new column
BLB_join[, dt_pformat := p_format]  # this keeps the MAVE row's p_format

# Calculate nucleotide differences only for rows with matching amino acid change
BLB_join[, diff_nt := (p_format == dt_pformat) & (
  (Chromosome != chr) |
    (Start      != start) |
    (Stop       != stop) |
    (ReferenceAlleleVCF != ref_plus) |
    (AlternateAlleleVCF != alt_plus)
)]

# Keep only rows where same AA but different NT
BLB_diff = BLB_join[diff_nt == TRUE]

# Compute score per MAVE row
BLB_scores = BLB_diff[, {
  hasB  = any(ClinicalSignificance == "B")
  hasLB = any(ClinicalSignificance == "LB")
  n_blb = sum(ClinicalSignificance %in% c("B","LB"))
  
  # Base penalty
  base  = if (hasB) -4L else if (hasLB) -2L else 0L
  
  # Flat extra penalty if >1 B/LB entry
  bonus = if (base < 0L && n_blb > 1L) -2L else 0L
  
  .(score = base + bonus)
}, by = i_row]

# -------- Combine scores --------
DT[, diff_nt_same_aa := NA_real_]

# Add P/LP and B/LB scores (NA treated as 0)
DT[PLP_scores$i_row, diff_nt_same_aa := fcoalesce(diff_nt_same_aa, 0) + PLP_scores$score]
DT[BLB_scores$i_row, diff_nt_same_aa := fcoalesce(diff_nt_same_aa, 0) + BLB_scores$score]

# Clean up
DT[, i_row := NULL]

# Replace maveData
maveData = as.data.frame(DT)

### CONDITION: DIFF AA CHANGE & GRANTHAM DISTANCE ###
### --- Ensure data.tables ---
setDT(maveData)
setDT(clinvar_PLP)
setDT(clinvar_BLB)

# Initialize column
maveData[, diff_aa_grantham_dis := NA_real_]

### --- P/LP scoring ---
# Join maveData with clinvar_PLP
PLP_matches = merge(
  maveData,
  clinvar_PLP,
  by.x = c("gene", "aa_pos", "ref_aa1"),
  by.y = c("GeneSymbol", "aa_pos", "ref_AA_1"),
  all.x = FALSE,
  all.y = FALSE,
  suffixes = c(".mave", ".clin")
)

# Keep rows where alt AA differs AND ClinVar has smaller Grantham distance
PLP_matches = PLP_matches[
  alt_AA_1 != alt_aa1 & grantham_distance.clin < grantham_distance.mave
]

# Score per MAVE variant
PLP_score = PLP_matches[, {
  hasP   = any(ClinicalSignificance == "P")
  hasLP  = any(ClinicalSignificance == "LP")
  n_plp  = sum(ClinicalSignificance %in% c("P","LP"))
  
  base   = if (hasP) 2L else if (hasLP) 1L else 0L
  bonus  = if (base > 0L && n_plp > 1L) 1L else 0L
  
  .(plp_score = base + bonus)
}, by = .(gene, aa_pos, ref_aa1, alt_aa1)]

# Add back to maveData (treat NA as 0)
maveData[PLP_score, diff_aa_grantham_dis := fcoalesce(diff_aa_grantham_dis, 0) + plp_score,
         on = .(gene, aa_pos, ref_aa1, alt_aa1)]

### --- B/LB scoring ---
# Join with clinvar_BLB
BLB_matches = merge(
  maveData,
  clinvar_BLB,
  by.x = c("gene", "aa_pos", "ref_aa1"),
  by.y = c("GeneSymbol", "aa_pos", "ref_AA_1"),
  all.x = FALSE,
  all.y = FALSE,
  suffixes = c(".mave", ".clin")
)

# Keep rows where alt AAs differ and ClinVar has larger Grantham distance
BLB_matches = BLB_matches[
  alt_AA_1 != alt_aa1 & grantham_distance.clin > grantham_distance.mave
]

# Score per MAVE variant
BLB_score = BLB_matches[, {
  hasB   = any(ClinicalSignificance == "B")
  hasLB  = any(ClinicalSignificance == "LB")
  n_blb  = sum(ClinicalSignificance %in% c("B","LB"))
  
  base   = if (hasB) -2L else if (hasLB) -1L else 0L
  bonus  = if (base < 0L && n_blb > 1L) -1L else 0L
  
  .(blb_score = base + bonus)
}, by = .(gene, aa_pos, ref_aa1, alt_aa1)]

# Add back to maveData (treat NA as 0)
maveData[BLB_score, diff_aa_grantham_dis := fcoalesce(diff_aa_grantham_dis, 0) + blb_score,
         on = .(gene, aa_pos, ref_aa1, alt_aa1)]

## FOR EACH PREDICTOR: REVEL ##
maveData$revel_eval <- dplyr::case_when(
  maveData$REVEL_score <= 0.016 ~ -4,
  maveData$REVEL_score > 0.016 & maveData$REVEL_score <= 0.052 ~ -3,
  maveData$REVEL_score > 0.052 & maveData$REVEL_score <= 0.184 ~ -2,
  maveData$REVEL_score > 0.184 & maveData$REVEL_score <= 0.290 ~ -1,
  maveData$REVEL_score > 0.290 & maveData$REVEL_score <= 0.643 ~  0,
  maveData$REVEL_score > 0.643 & maveData$REVEL_score <= 0.772 ~  1,
  maveData$REVEL_score > 0.772 & maveData$REVEL_score <= 0.878 ~  2,
  maveData$REVEL_score > 0.878 & maveData$REVEL_score <= 0.931 ~  3,
  maveData$REVEL_score > 0.931 ~  4,
  TRUE ~ NA_real_
)

## FOR EACH PREDICTOR: BAYESDEL ##
maveData$bayesdel_eval <- dplyr::case_when(
  maveData$BayesDel_noAF_score <= -0.52 ~ -3,
  maveData$BayesDel_noAF_score > -0.52 & maveData$BayesDel_noAF_score <= -0.34 ~ -2,
  maveData$BayesDel_noAF_score > -0.34 & maveData$BayesDel_noAF_score <= -0.16 ~ -1,
  maveData$BayesDel_noAF_score > -0.16 & maveData$BayesDel_noAF_score <= 0.12 ~  0,
  maveData$BayesDel_noAF_score > 0.12 & maveData$BayesDel_noAF_score <= 0.26 ~  1,
  maveData$BayesDel_noAF_score > 0.26 & maveData$BayesDel_noAF_score <= 0.40 ~  2,
  maveData$BayesDel_noAF_score > 0.40 & maveData$BayesDel_noAF_score <= 0.49 ~  3,
  maveData$BayesDel_noAF_score > 0.49 ~  4,
  TRUE ~ NA_real_
)

## FOR EACH PREDICTOR: MUTPRED2 ##
maveData$mutpred_eval <- dplyr::case_when(
  maveData$MutPred_score <= 0.010 ~ -4,
  maveData$MutPred_score > 0.010 & maveData$MutPred_score <= 0.0318 ~ -3,
  maveData$MutPred_score > 0.0318 & maveData$MutPred_score <= 0.197 ~ -2,
  maveData$MutPred_score > 0.197 & maveData$MutPred_score <= 0.391 ~ -1,
  maveData$MutPred_score > 0.391 & maveData$MutPred_score <= 0.736 ~  0,
  maveData$MutPred_score > 0.736 & maveData$MutPred_score <= 0.828 ~  1,
  maveData$MutPred_score > 0.828 & maveData$MutPred_score <= 0.894 ~  2,
  maveData$MutPred_score > 0.894 & maveData$MutPred_score <= 0.931 ~  3,
  maveData$MutPred_score > 0.931 ~  4,
  TRUE ~ NA_real_
)

## FOR EACH PREDICTOR: VEST ##
maveData$vest_eval <- dplyr::case_when(
  maveData$VEST4_score <= 0.077 ~ -3,
  maveData$VEST4_score > 0.077 & maveData$VEST4_score <= 0.302 ~ -2,
  maveData$VEST4_score > 0.302 & maveData$VEST4_score <= 0.449 ~ -1,
  maveData$VEST4_score > 0.449 & maveData$VEST4_score <= 0.763 ~  0,
  maveData$VEST4_score > 0.763 & maveData$VEST4_score <= 0.860 ~  1,
  maveData$VEST4_score > 0.860 & maveData$VEST4_score <= 0.908 ~  2,
  maveData$VEST4_score > 0.908 & maveData$VEST4_score <= 0.964 ~  3,
  maveData$VEST4_score > 0.964 ~  4,
  TRUE ~ NA_real_
)

## CALCULATE NEW_ACMG_SCORE
maveData$new_ACMG_score = rowSums(
  maveData[, c("diff_aa_grantham_dis", "revel_eval")],
  na.rm = TRUE
)

# Replace rows where all three were NA
# Compute rowSums ignoring NA
maveData$new_ACMG_score = rowSums(
  as.matrix(maveData[, c("diff_aa_grantham_dis", "revel_eval")]),
  na.rm = TRUE
)

# Identify rows where all three are NA
all_na = rowSums(!is.na(as.matrix(maveData[, c("diff_nt_same_aa", "diff_aa_grantham_dis", "revel_eval")]))) == 0

# Replace those with NA
maveData$new_ACMG_score[all_na] = NA_real_ 

# Ensure data.tables
setDT(maveData)
setDT(clinvarTable)

# Create in_clinvar column by exact match
maveData[, in_clinvar := (chr %in% clinvarTable$Chromosome) &
           (start %in% clinvarTable$Start) &
           (stop %in% clinvarTable$Stop) &
           (ref_plus %in% clinvarTable$ReferenceAlleleVCF) &
           (alt_plus %in% clinvarTable$AlternateAlleleVCF)]

write.csv(clinvarTable, "new_acmg/clinvarTable_lookup.csv", row.names = F) 
write.csv(maveData, "new_acmg/new_ACMG_scoring_results.csv", row.names = F)
