# Set working directory
setwd("/Volumes/T9/MAVE Data Project")

# Load libraries 
library(data.table)  
library(dplyr)     
library(pROC)        
library(tidyr)
library(viridis)
library(ggplot2)

# Read in secondary analysis data 
maveData = fread("secondary_results/final.csv")

# Only FUNC vs LOF
maveData_binary = maveData %>%
  filter(func_classification %in% c("FUNC", "LOF")) %>%
  mutate(
    lof_label = ifelse(func_classification == "LOF", 1, 0)
  )

# Convert VEP score columns to numeric and replace "." with NA
vep_cols = c("BayesDel_noAF_score", "VEST4_score", "REVEL_score",
              "MutPred_score", "am_pathogenicity", "VARITY_R_score")

for (col in vep_cols) {
  maveData_binary[[col]] = as.numeric(ifelse(maveData_binary[[col]] == ".", NA, maveData_binary[[col]]))
}

# Compute AUC per gene for each VEP
auc_per_gene = maveData_binary %>%
  filter(!is.na(lof_label)) %>%  # ensure label exists
  group_by(gene) %>%
  summarise(across(all_of(vep_cols),
                   ~if(sum(!is.na(.x)) > 0) {  # require at least one non-NA
                     as.numeric(pROC::auc(lof_label, .x))
                   } else {
                     NA_real_
                   },
                   .names = "{.col}_AUC"),
            .groups = "drop")

# Reshape AUC table to long format
auc_long = auc_per_gene %>%
  pivot_longer(
    cols = -gene,
    names_to = "tool",
    values_to = "AUC"
  )

# Order genes by mean AUC (makes patterns easier to see)
gene_order = auc_long %>%
  group_by(gene) %>%
  summarise(mean_auc = mean(AUC, na.rm = TRUE)) %>%
  arrange(mean_auc) %>%
  pull(gene)

auc_long$gene = factor(auc_long$gene, levels = gene_order)

# Clean tool labels for plotting
auc_long$tool = recode(
  auc_long$tool,
  BayesDel_noAF_score_AUC = "BayesDel",
  VEST4_score_AUC    = "VEST4",
  REVEL_score_AUC    = "REVEL",
  MutPred_score_AUC  = "MutPred2",
  am_pathogenicity_AUC       = "AlphaMissense",
  VARITY_R_score_AUC   = "VARITY"
)

# Plot: heatmap
ggplot(auc_long, aes(x = tool, y = gene, fill = AUC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  
  # Add AUC labels ONLY for weak performance (<0.7)
  # geom_text(
  #   data = subset(auc_long, !is.na(AUC) & AUC < 0.7),
  #   aes(label = sprintf("%.2f", AUC)),
  #   size = 2.3,
  #   color = "white"
  # ) +
  # 
  # Color gradient with viridis, midpoint at 0.7
  scale_fill_viridis(
    option = "mako",
    direction = 1,
    na.value = "grey85",
    limits = c(0.58, 1),
    oob = scales::squish,
    name = "AUC"
  ) +
  
  labs(
    x = "Variant Effect Predictor (VEP)",
    y = "Gene"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Read in secondary analysis data 
maveData = fread("Revisions/Secondary/by_variant_secondary_results(06Jan2026).csv")

# Only FUNC vs LOF
maveData_clinvar = maveData %>%
  filter(func_classification %in% c("FUNC", "LOF")) %>%
  mutate(
    lof_label = ifelse(func_classification == "LOF", 1, 0)
  )

# Filter for data where MAVE and ClinVar are concordant 
# LOF = LP/P, FUNC = LB/B
maveData_clinvar = maveData_clinvar %>%
  filter(
    (!is.na(clinical_significance_sim)) & 
      ((lof_label == 1 & clinical_significance_sim == "P/LP") |
         (lof_label == 0 & clinical_significance_sim == "B/LB"))
  )

# Convert VEP score columns to numeric and replace "." with NA
vep_cols = c("BayesDel_noAF_score", "VEST4_score", "REVEL_score",
             "MutPred_score", "am_pathogenicity", "VARITY_R_score")

for (col in vep_cols) {
  maveData_clinvar[[col]] = as.numeric(
    ifelse(maveData_clinvar[[col]] == ".", NA, maveData_clinvar[[col]])
  )
}

# Compute AUC per gene for each VEP
auc_per_gene = maveData_clinvar %>%
  filter(!is.na(lof_label)) %>%  # ensure label exists
  group_by(gene) %>%
  summarise(across(all_of(vep_cols),
                   ~if(length(unique(lof_label[!is.na(.x)])) == 2) {  # require both classes
                     as.numeric(pROC::auc(lof_label, .x))
                   } else {
                     NA_real_
                   },
                   .names = "{.col}_AUC"),
            .groups = "drop")

# Reshape AUC table to long format
auc_long = auc_per_gene %>%
  pivot_longer(
    cols = -gene,
    names_to = "tool",
    values_to = "AUC"
  )

# Order genes by mean AUC (makes patterns easier to see)
gene_order = auc_long %>%
  group_by(gene) %>%
  summarise(mean_auc = mean(AUC, na.rm = TRUE)) %>%
  arrange(mean_auc) %>%
  pull(gene)

auc_long$gene = factor(auc_long$gene, levels = gene_order)

# Clean tool labels for plotting
auc_long$tool = recode(
  auc_long$tool,
  BayesDel_noAF_score_AUC = "BayesDel",
  VEST4_score_AUC    = "VEST4",
  REVEL_score_AUC    = "REVEL",
  MutPred_score_AUC  = "MutPred2",
  am_pathogenicity_AUC       = "AlphaMissense",
  VARITY_R_score_AUC   = "VARITY"
)

# Number of variants in ClinVar-concordant subset
n_variants_concordant = nrow(maveData_clinvar)
cat("Number of ClinVar-concordant variants:", n_variants_concordant, "\n")

# Number of genes in ClinVar-concordant subset
genes_concordant = maveData_clinvar %>%
  group_by(gene_only) %>%
  summarise(n = n()) %>%
  nrow()
cat("Number of genes with ClinVar-concordant variants:", genes_concordant, "\n")

# Alternatively, get summary by gene
variant_summary = maveData_clinvar %>%
  group_by(gene) %>%
  summarise(n_variants = n()) %>%
  arrange(desc(n_variants))

print(variant_summary)

# Plot: heatmap
ggplot(auc_long, aes(x = tool, y = gene, fill = AUC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  
  # Add AUC labels ONLY for weak performance (<0.7)
  # geom_text(
  #   data = subset(auc_long, !is.na(AUC) & AUC < 0.7),
  #   aes(label = sprintf("%.2f", AUC)),
  #   size = 2.3,
  #   color = "white"
  # ) +
  # 
  # Color gradient with viridis, midpoint at 0.7
  scale_fill_viridis(
    option = "mako",
    direction = 1,
    na.value = "grey85",
    limits = c(0.50, 1),
    oob = scales::squish,
    name = "AUC"
  ) +
  
  labs(
    x = "Variant Effect Predictor (VEP)",
    y = "Gene"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Read in secondary analysis data 
maveData = fread("Revisions/Secondary/by_variant_secondary_results(06Jan2026).csv")

# Only FUNC vs LOF and create lof_label
maveData_protein = maveData %>%
  filter(func_classification %in% c("FUNC", "LOF")) %>%
  mutate(
    lof_label = ifelse(func_classification == "LOF", 1, 0)
  )

# Make a new column protein_domain = 1 if interpro_domain is NOT NA, else 0
maveData_protein = maveData_protein %>%
  mutate(
    protein_domain = ifelse(!is.na(interpro_domain), 1, 0)
  )

# List of VEP score columns
vep_cols = c("BayesDel_noAF_score", "VEST4_score", "REVEL_score",
             "MutPred_score", "am_pathogenicity", "VARITY_R_score")

# Convert VEP scores to numeric and replace "." with NA
for (col in vep_cols) {
  maveData_protein[[col]] = as.numeric(ifelse(maveData_protein[[col]] == ".", NA, maveData_protein[[col]]))
}

# Function to compute AUC safely (requires both classes present)
compute_auc = function(response, score) {
  if (length(unique(response[!is.na(score)])) == 2) {
    as.numeric(pROC::auc(response, score))
  } else {
    NA_real_
  }
}

# Compute AUC per VEP, separately for protein_domain = 0 vs 1
auc_by_domain = maveData_protein %>%
  group_by(protein_domain) %>%
  summarise(across(all_of(vep_cols),
                   ~compute_auc(lof_label, .x),
                   .names = "{.col}_AUC"),
            .groups = "drop") %>%
  pivot_longer(
    cols = -protein_domain,
    names_to = "tool",
    values_to = "AUC"
  ) %>%
  mutate(
    # Clean tool labels for plotting
    tool = recode(
      tool,
      BayesDel_noAF_score_AUC = "BayesDel",
      VEST4_score_AUC = "VEST4",
      REVEL_score_AUC = "REVEL",
      MutPred_score_AUC = "MutPred2",
      am_pathogenicity_AUC = "AlphaMissense",
      VARITY_R_score_AUC = "VARITY"
    ),
    # Label protein domain groups
    protein_domain = factor(protein_domain, levels = c(0,1),
                            labels = c("Outside Domain", "Inside Domain"))
  )

# Double bar graph: protein_domain vs VEP tool
ggplot(auc_by_domain, aes(x = tool, y = AUC, fill = protein_domain)) +
  
  # Bar geometry, dodge to put two bars side by side per VEP
  geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7) +

  
  # Viridis color scale for domain groups
  scale_fill_viridis(
    option = "mako",
    discrete = TRUE,
    direction = 1,
    name = "Protein Domain"
  ) +
  
  # Axis labels and title
  labs(
    x = "Variant Effect Predictor (VEP)",
    y = "AUC",
    title = "AUC by VEP tool and protein domain"
  ) +
  
  # Theme adjustments
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

