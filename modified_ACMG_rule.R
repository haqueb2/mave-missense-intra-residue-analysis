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
geneList = unique(maveData$gene)

# Make sure maveData is a data.table
maveData = as.data.table(maveData)

# Create new column
maveData$new_ACMG_score_mod = 0

# Select only needed columns and exclude unwanted rows
maveTable = maveData %>%
  dplyr::select(
    gene,
    protein_variant,
    func_classification,
    variant_category,
    aa_pos,
    variant_count_per_residue,
    p_format,
    ref_aa3,
    alt_aa3,
    ref_aa1,
    alt_aa1,
    grantham_distance,
    chr,
    start,
    stop,
    ref,
    alt,
    REVEL_score,
    REVEL_pred,
    ref_plus,
    alt_plus,
  ) %>%
  dplyr::filter(
    !is.na(protein_variant)
  )

# Get MAVE LOF table
maveLOF = maveTable %>%
  filter(func_classification %in% c("LOF"))

# Get MAVE FUNC table
maveFUNC = maveTable %>%
  filter(func_classification %in% c("FUNC"))

### CONDITION : DIFF NT, SAME AA CHANGE LOF AND FUNC ###
# Convert to data.table
DT = as.data.table(maveData)
LOF = as.data.table(maveLOF)
FUNC = as.data.table(maveFUNC)

# Add row index to join back scores later
DT[, i_row := .I]

# -------- LOF --------
setkey(DT, p_format)
setkey(LOF, p_format)

# Join on amino acid change
LOF_join = LOF[DT, nomatch=0, allow.cartesian=TRUE]

# Flag rows that differ at the nucleotide level
LOF_join[maveData, 
         diff_nt := (chr != i.chr) | 
           (start != i.start) | 
           (stop != i.stop) | 
           (ref_plus != i.ref_plus) | 
           (alt_plus != i.alt_plus),
         on = .(p_format)]

LOF_diff = LOF_join[diff_nt == TRUE]

# Compute score per maveData row
# Base: LOF = 4, else 0
# Bonus: +2 if there is >1 matching LOF entry (flat, not per extra)
LOF_scores <- LOF_diff[, .(
  score = fcase(
    any(func_classification == "LOF") & .N > 1, 6L, # 4 + 2
    any(func_classification == "LOF"),              4L,
    default = 0L
  )
), by = i_row]

# -------- FUNC --------
setkey(FUNC, p_format)
FUNC_join = FUNC[DT, nomatch=0, allow.cartesian=TRUE]

# Flag nucleotide differences
FUNC_join[maveData,
         diff_nt := (chr != i.chr) |
           (start != i.start) |
           (stop != i.stop) |
           (ref_plus != i.ref_plus) |
           (alt_plus != i.alt_plus),
         on = .(p_format)]

FUNC_diff = FUNC_join[diff_nt == TRUE]

FUNC_scores = FUNC_diff[, .(
  score = fcase(
    sum(func_classification == "FUNC") > 1, -6L, # -4 + -2
    any(func_classification == "FUNC"),    -4L,
    default = 0L
  )
), by = i_row]

# -------- Combine scores --------
# initialize the column
DT[, diff_nt_same_aa := NA_real_]

# Add LOF and FUNC scores (treat NA as 0)
DT[LOF_scores, on = .(i_row), diff_nt_same_aa := fcoalesce(diff_nt_same_aa, 0) + i.score]
DT[FUNC_scores, on = .(i_row), diff_nt_same_aa := fcoalesce(diff_nt_same_aa, 0) + i.score]

# Clean up
DT[, i_row := NULL]

# Replace maveData
maveData = as.data.frame(DT)

### CONDITION: DIFF AA CHANGE & GRANTHAM DISTANCE ###
# Ensure data.tables
setDT(maveData)
setDT(maveLOF)

# Initialize score column
maveData[, diff_aa_grantham_dis := 0L]

# --- Merge on gene, aa_pos, ref_aa1 to compare same position ---
LOF_matches <- merge(
  maveData,
  maveLOF,
  by = c("gene", "aa_pos", "ref_aa1"),
  allow.cartesian = TRUE,
  suffixes = c(".mave", ".LOF")
)

# --- Keep rows where alt AA differs and MAVE Grantham < LOF Grantham ---
LOF_matches <- LOF_matches[
  alt_aa1.mave != alt_aa1.LOF & 
    grantham_distance.mave < grantham_distance.LOF
]

# --- Compute LOF score per MAVE variant ---
LOF_scores <- LOF_matches[, .(
  LOF_score = fcase(
    sum(func_classification.LOF == "LOF") > 1, 3L,  # 2 + 1 bonus
    any(func_classification.LOF == "LOF"),           2L,
    default = 0L
  )
), by = .(gene, aa_pos, ref_aa1, alt_aa1.mave)]

# --- Join back to main MAVE data ---
maveData[LOF_scores, 
         on = .(gene, aa_pos, ref_aa1, alt_aa1 = alt_aa1.mave),
         diff_aa_grantham_dis := fcoalesce(diff_aa_grantham_dis, 0L) + i.LOF_score]

### --- FUNC scoring ---
# Merge MAVE data with FUNC table
FUNC_matches <- merge(
  maveData,
  maveFUNC,
  by = c("gene", "aa_pos", "ref_aa1"),
  allow.cartesian = TRUE,
  suffixes = c(".mave", ".FUNC")
)

# Keep rows where alt AA differs and MAVE Grantham < FUNC Grantham
FUNC_matches <- FUNC_matches[
  alt_aa1.mave != alt_aa1.FUNC & grantham_distance.FUNC > grantham_distance.mave
]

# Compute FUNC score per MAVE variant
FUNC_score <- FUNC_matches[, .(
  FUNC_score = fcase(
    sum(func_classification.FUNC == "FUNC") > 1, -3L,  # -2 base + -1 bonus
    any(func_classification.FUNC == "FUNC"),           -2L,
    default = 0L
  )
), by = .(gene, aa_pos, ref_aa1, alt_aa1.mave)]

# Add back to maveData
maveData[FUNC_score, 
         on = .(gene, aa_pos, ref_aa1, alt_aa1 = alt_aa1.mave),
         diff_aa_grantham_dis := fcoalesce(diff_aa_grantham_dis, 0L) + i.FUNC_score]

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
maveData$new_ACMG_score_mod = rowSums(
  maveData[, c("diff_nt_same_aa", "diff_aa_grantham_dis", "revel_eval")],
  na.rm = TRUE
)

# Replace rows where all three were NA
# Compute rowSums ignoring NA
maveData$new_ACMG_score_mod = rowSums(
  as.matrix(maveData[, c("diff_nt_same_aa", "diff_aa_grantham_dis", "revel_eval")]),
  na.rm = TRUE
)

# Identify rows where all three are NA
all_na = rowSums(!is.na(as.matrix(maveData[, c("diff_nt_same_aa", "diff_aa_grantham_dis", "revel_eval")]))) == 0

# Replace those with NA
maveData$new_ACMG_score_mod[all_na] = NA_real_

write.csv(maveTable,"new_acmg/maveTable_lookup.csv", row.names = F)
write.csv(maveData, "new_acmg/new_ACMG_scoring_MAVE_modified.csv", row.names = F)