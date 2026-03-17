# Set working direcory
setwd("/Volumes/T9/MAVE Data Project")

# Load libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(data.table)
library(scales)
library(pROC)

maveData = fread("secondary_results/final.csv")

# Read in MANE_Select transcripts 
maneSelect = read_excel("annotations/MANE_Select_NM.xlsx")

# ------------------------------------------------------------
# Add canonical column
# ------------------------------------------------------------

# Step 0: Initialize canonical as logical FALSE for all rows
maveData = maveData %>%
  mutate(canonical = FALSE)

# Step 1: Set canonical = TRUE if MANE_SELECT transcript
# matches an NM number in the MANE Select reference
maveData = maveData %>%
  mutate(
    canonical = MANE_SELECT %in% maneSelect$NM_number
  )

# Step 2: For rows not already canonical,
# set canonical = TRUE if MANE_SELECT starts with "E"
# (e.g. Ensembl-style transcripts)
maveData = maveData %>%
  mutate(
    canonical = if_else(
      canonical == FALSE &
        !is.na(MANE_SELECT) &
        str_starts(MANE_SELECT, "E"),
      TRUE,
      canonical,
      missing = canonical
    )
  )

# Step 3: For rows still not canonical AND with missing MANE_SELECT,
# set canonical = TRUE if the transcript matches
# a MANE Select transcript format
maveData = maveData %>%
  mutate(
    canonical = if_else(
      canonical == FALSE &
        is.na(MANE_SELECT) &
        transcript %in% maneSelect$format,
      TRUE,
      canonical,
      missing = canonical
    )
  )

## GRANTHAM DISTANCE ANALYSIS ##
# Replace "." with NA and convert grantham_distance column to numeric
maveData$grantham_distance = as.numeric(ifelse(maveData$grantham_distance == ".", NA, maveData$grantham_distance))

# Summarize Grantham distance per gene and aa_pos
grantham_summary = maveData %>%
  filter(!is.na(grantham_distance)) %>%
  group_by(gene, aa_pos) %>%
  summarize(
    min_LOF = ifelse(any(func_classification == "LOF"),
                     min(grantham_distance[func_classification == "LOF"], na.rm = TRUE),
                     NA_real_),
    max_FUNC = ifelse(any(func_classification == "FUNC"),
                      max(grantham_distance[func_classification == "FUNC"], na.rm = TRUE),
                      NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    LOF_higher = min_LOF > max_FUNC
  )

# Summarize per gene
gene_summary = grantham_summary %>%
  group_by(gene) %>%
  summarize(
    total_positions = n(),
    LOF_higher_count = sum(LOF_higher, na.rm = TRUE),
    LOF_higher_pct = LOF_higher_count / total_positions * 100,
    .groups = "drop"
  )

# Print table
print(gene_summary)

# Bar plot: % of positions where LOF > FUNC per gene
ggplot(gene_summary, aes(x = reorder(gene, -LOF_higher_pct), y = LOF_higher_pct)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(
    x = "Gene",
    y = "% Positions (min Grantham distance LOF > max Grantham distance FUNC)"
  ) +
  theme(
    panel.grid = element_blank(),                # remove grid lines
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,
      face = "italic"                            # italicize gene names
    )
  )

## IN SILICO SCORES COMPARED TO MAVE CLASS ##
plot_data = maveData %>%
  mutate(
    VEST4_score = as.numeric(ifelse(VEST4_score == ".", NA, VEST4_score)),
    REVEL_score = as.numeric(ifelse(REVEL_score == ".", NA, REVEL_score)),
    MutPred_score = as.numeric(ifelse(MutPred_score == ".", NA, MutPred_score)),
    BayesDel_noAF_score = as.numeric(ifelse(BayesDel_noAF_score == ".", NA, BayesDel_noAF_score))
  ) %>%
  select(func_classification, VEST4_score, REVEL_score, MutPred_score, BayesDel_noAF_score) %>%
  pivot_longer(
    cols = -func_classification,
    names_to = "Predictor",
    values_to = "Score"
  ) %>%
  filter(!is.na(Score))

# Thresholds table
thresholds = data.frame(
  Predictor = c("VEST4_score", "REVEL_score", "MutPred_score", "BayesDel_noAF_score"),
  Threshold = c(0.764, 0.644, 0.737, 0.13)
)

ggplot(plot_data, aes(x = func_classification, y = Score)) +
  geom_violin(trim = FALSE, fill = "lightblue", color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  facet_wrap(~ Predictor, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "MAVE Functional Classification",
    y = "VEP Score"
  ) +
  theme(
    panel.grid = element_blank(),                   # removes major & minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)  
  ) +
  geom_hline(
    data = thresholds,
    aes(yintercept = Threshold),
    linetype = "dashed",
    color = "grey"
  )

## PLOTTING CONCORDANCE OF IN SILICO TOOLS AND MAVE CLASS ##
# Filter: only FUNC or LOF, exclude INT
maveData_binary = maveData %>%
  filter(func_classification %in% c("FUNC", "LOF")) %>%
  mutate(
    func_classification = factor(func_classification, levels = c("FUNC", "LOF")),
    BayesDel_noAF_score = factor(BayesDel_noAF_pred, levels = c("T", "D")),
    VEST4_pred = factor(VEST4_pred, levels = c("T", "D")),
    REVEL_pred = factor(REVEL_pred, levels = c("T", "D")),
    MutPred_pred = factor(MutPred_pred, levels = c("T", "D"))
  )

# Reshape data for plotting
confusion_counts_all = maveData_binary %>%
  pivot_longer(
    cols = c(BayesDel_noAF_score, VEST4_pred, REVEL_pred, MutPred_pred),
    names_to = "tool",
    values_to = "prediction"
  ) %>%
  filter(!is.na(prediction)) %>%
  group_by(tool, func_classification) %>%
  count(prediction, name = "count") %>%
  group_by(tool) %>%                       # <<< key change
  mutate(
    total = sum(count),                    # total per tool
    percent = round((count / total) * 100, 1)
  ) %>%
  ungroup()

# Heatmap with counts
ggplot(confusion_counts_all, aes(x = func_classification, y = prediction, fill = percent)) +
  geom_tile(color = "white") +
  
  # Counts (bigger, on top)
  geom_text(
    aes(label = comma(count)),
    color = "black",
    size = 4,
    vjust = -0.2
  ) +
  
  # Percentages (slightly smaller, below counts)
  geom_text(
    aes(label = paste0("(", percent, "%)")),
    color = "black",
    size = 3.5,
    vjust = 1.5
  ) +
  
  scale_fill_gradient(low = "lightyellow", high = "red") +
  labs(
    x = "MAVE Classification",
    y = "Prediction",
    fill = "Percent"
  ) +
  facet_wrap(~tool) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10)
  )

## ADDITIONAL IN SILICO SCORES COMPARED TO MAVE CLASS ##
plot_data = maveData %>%
  # Convert "." to NA and coerce to numeric
  mutate(
    AM_score = as.numeric(ifelse(am_pathogenicity == ".", NA, am_pathogenicity)),
    VARITY_R_score = as.numeric(ifelse(VARITY_R_score == ".", NA, VARITY_R_score))
  ) %>%
  select(func_classification, AM_score, VARITY_R_score) %>%
  pivot_longer(
    cols = -func_classification,
    names_to = "Predictor",
    values_to = "Score"
  ) %>%
  filter(!is.na(Score))

# Thresholds for plotting (T and D for AM, D for VARITY)
thresholds = data.frame(
  Predictor = c("AM_score", "AM_score", "VARITY_R_score"),
  Threshold = c(0.34, 0.564, 0.75)
)

# Violin plot for AM and VARITY scores
ggplot(plot_data, aes(x = func_classification, y = Score)) +
  geom_violin(trim = FALSE, fill = "lightgreen", color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  facet_wrap(~ Predictor, scales = "free_y") +
  geom_hline(data = thresholds, aes(yintercept = Threshold), linetype = "dashed", color = "grey") +
  labs(x = "MAVE Functional Classification", y = "Score") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

# Prediction labels
maveData_binary = maveData %>%
  filter(func_classification %in% c("FUNC", "LOF")) %>%
  mutate(
    func_classification = factor(func_classification, levels = c("FUNC", "LOF")),
    # AM_pred with 3 classes: T / A / D
    AM_pred = case_when(
      is.na(am_pathogenicity) ~ NA_character_,  # keep NA as NA
      am_pathogenicity > 0.564 ~ "D",
      am_pathogenicity < 0.34  ~ "T",
      TRUE                     ~ "A"
    ),
    VARITY_pred = ifelse(VARITY_R_score > 0.75, "D", "T"),
    AM_pred = factor(AM_pred, levels = c("T", "A", "D")),
    VARITY_pred = factor(VARITY_pred, levels = c("T", "D"))
  )

# Reshape data and create full grids for each tool separately
confusion_counts_all <- maveData_binary %>%
  pivot_longer(
    cols = c(AM_pred, VARITY_pred),
    names_to = "tool",
    values_to = "prediction"
  ) %>%
  filter(!is.na(prediction)) %>%
  group_by(tool, func_classification, prediction) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tool) %>%
  group_modify(~{
    if (.y$tool == "AM_pred") {
      complete(.x,
               func_classification = c("FUNC", "LOF"),
               prediction = c("T", "A", "D"),
               fill = list(count = 0)) %>%
        mutate(prediction = factor(prediction, levels = c("T", "A", "D")))
    } else {
      complete(.x,
               func_classification = c("FUNC", "LOF"),
               prediction = c("T", "D"),
               fill = list(count = 0)) %>%
        mutate(prediction = factor(prediction, levels = c("T", "D")))
    }
  }) %>%
  ungroup() %>%
  group_by(tool) %>%                         # <<< percent of total per tool
  mutate(
    total = sum(count),
    percent = round(count / total * 100, 1)
  ) %>%
  ungroup()


# Heatmap
ggplot(confusion_counts_all, aes(x = func_classification, y = prediction, fill = percent)) +
  geom_tile(color = "white") +
  # Two geom_text layers: one for count, one for percent (slightly smaller)
  geom_text(
    aes(label = comma(count)),
    color = "black",
    size = 4,
    vjust = -0.2
  ) +
  geom_text(
    aes(label = paste0("(", percent, "%)")),
    color = "black",
    size = 3.5,
    vjust = 1.5
  ) +
  scale_fill_gradient(low = "lightyellow", high = "red") +
  facet_wrap(~tool, scales = "free_y") +
  labs(
    x = "MAVE Classification",
    y = "Prediction",
    fill = "Percent"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # facet labels bigger
    axis.text.x = element_text(size = 10),                # bigger x tick labels
    axis.text.y = element_text(size = 10)                 # bigger y tick labels
  )

## MAKING ROC CURVE WITH ALL 6 PREDICTORS ##
# Keep numeric scores for all predictors 
roc_data = maveData %>%
  filter(func_classification %in% c("FUNC", "LOF")) %>%
  mutate(
    lof_label = ifelse(func_classification == "LOF", 1, 0),
    VEST4_score       = as.numeric(ifelse(VEST4_score == ".", NA, VEST4_score)),
    REVEL_score       = as.numeric(ifelse(REVEL_score == ".", NA, REVEL_score)),
    MutPred_score     = as.numeric(ifelse(MutPred_score == ".", NA, MutPred_score)),
    BayesDel_noAF_score = as.numeric(ifelse(BayesDel_noAF_score == ".", NA, BayesDel_noAF_score)),
    AM_score          = as.numeric(ifelse(am_pathogenicity == ".", NA, am_pathogenicity)),
    VARITY_R_score    = as.numeric(ifelse(VARITY_R_score == ".", NA, VARITY_R_score))
  )

# Build ROC curves for all predictors 
roc_list = list(
  BayesDel = roc(roc_data$lof_label, roc_data$BayesDel_noAF_score, quiet = TRUE),
  VEST4    = roc(roc_data$lof_label, roc_data$VEST4_score, quiet = TRUE),
  REVEL    = roc(roc_data$lof_label, roc_data$REVEL_score, quiet = TRUE),
  MutPred  = roc(roc_data$lof_label, roc_data$MutPred_score, quiet = TRUE),
  AM       = roc(roc_data$lof_label, roc_data$AM_score, quiet = TRUE),
  VARITY   = roc(roc_data$lof_label, roc_data$VARITY_R_score, quiet = TRUE)
)

# Cover into a tidy dataframe for plotting
roc_df = bind_rows(
  lapply(names(roc_list), function(name) {
    data.frame(
      tool = name,
      fpr = 1 - roc_list[[name]]$specificities,
      tpr = roc_list[[name]]$sensitivities,
      auc = round(pROC::auc(roc_list[[name]]), 3)
    )
  })
)

# Plot all ROC curves into one ggplot and include AUC in legend 
roc_df = roc_df %>%
  group_by(tool) %>%
  mutate(tool_auc = paste0(tool, " (AUC=", unique(auc), ")"))

# Create a named color palette matching the tool_auc values
color_map = setNames(
  c("orange", "deepskyblue", "springgreen3", "hotpink", "purple", "red"),
  unique(roc_df$tool_auc)
)

# Plot with manual colors
ggplot(roc_df, aes(x = fpr, y = tpr, color = tool_auc)) +
  geom_line(size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray") +
  scale_color_manual(values = color_map) +   # <-- THIS LINE CONTROLS COLORS
  labs(
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Predictor"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2),
    axis.text = element_text(size = 10),
    axis.title = element_text(margin = margin(t = 12))
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

## PLOTTING CONCORDANCE OF CLINVAR CLASS AND MAVE CLASS ##
# Create new column 
maveData$clinical_significance_sim = NA

# Assign P/LP if "pathogenic" is present but NOT "pathogenicity"
maveData$clinical_significance_sim[grepl("pathogenic", maveData$clinical_significance, ignore.case = TRUE) & 
                                     !grepl("pathogenicity", maveData$clinical_significance, ignore.case = TRUE)] = "P/LP"

# Assign B/LB if "benign" is present
maveData$clinical_significance_sim[grepl("benign", maveData$clinical_significance, ignore.case = TRUE)] = "B/LB"

# Prepare data for plotting
maveData_binary = maveData %>%
  filter(func_classification %in% c("FUNC", "LOF", "INT")) %>%
  mutate(
    func_classification = factor(func_classification, levels = c("LOF", "INT", "FUNC")),
    clinical_significance_sim = factor(clinical_significance_sim, levels = c("P/LP", "B/LB"))
  )

# Reshape data for plotting
confusion_counts_all = maveData_binary %>%
  filter(!is.na(clinical_significance_sim)) %>%
  group_by(clinical_significance_sim, func_classification) %>%
  summarise(count = n(), .groups = "drop") %>%   # <<< drop all grouping
  mutate(
    total = sum(count),                          # <<< global total
    percent = round(count / total * 100, 1)
  )

# Heatmap with percentages
ggplot(confusion_counts_all, aes(x = func_classification, y = clinical_significance_sim, fill = percent)) +
  geom_tile(color = "white") +
  # Count on top
  geom_text(aes(label = comma(count)), color = "black", size = 4, vjust = -0.2) +
  # Percent below in smaller text
  geom_text(aes(label = paste0("(", percent, "%)")), color = "black", size = 3.5, vjust = 1.5) +
  scale_fill_gradient(low = "lightyellow", high = "red") +
  labs(
    x = "MAVE Classification",
    y = "ClinVar Classification",
    fill = "Percentage"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),      # bigger tick labels
    axis.title.x = element_text(margin = margin(t = 20)),  # space above x-axis title
    axis.title.y = element_text(margin = margin(r = 20))   # space right of y-axis title
  )

## PLOTTING CONCORDANCE OF CLINVAR CLASS (≥2 STARS) AND MAVE CLASS ##
# Create simplified ClinVar classification
maveData$clinical_significance_sim = NA

# Assign P/LP if "pathogenic" is present but NOT "pathogenicity"
maveData$clinical_significance_sim[
  grepl("pathogenic", maveData$clinical_significance, ignore.case = TRUE) &
    !grepl("pathogenicity", maveData$clinical_significance, ignore.case = TRUE)
] = "P/LP"

# Assign B/LB if "benign" is present
maveData$clinical_significance_sim[
  grepl("benign", maveData$clinical_significance, ignore.case = TRUE)
] = "B/LB"


# Prepare data for plotting
maveData_binary = maveData %>%
  filter(
    func_classification %in% c("FUNC", "LOF", "INT"),
    !is.na(clinical_significance_sim),
    star_rating >= 2                     # <<< KEY reviewer-driven filter
  ) %>%
  mutate(
    func_classification = factor(func_classification,
                                 levels = c("LOF", "INT", "FUNC")),
    clinical_significance_sim = factor(clinical_significance_sim,
                                       levels = c("P/LP", "B/LB"))
  )


# Summarise counts and compute percent of total
confusion_counts_all = maveData_binary %>%
  group_by(clinical_significance_sim, func_classification) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    total = sum(count),                  # global total (≥2-star variants)
    percent = round(count / total * 100, 1)
  )


# Heatmap
ggplot(confusion_counts_all, aes(x = func_classification, y = clinical_significance_sim, fill = percent)) +
  geom_tile(color = "white") +
  # Count on top
  geom_text(aes(label = comma(count)), color = "black", size = 4, vjust = -0.2) +
  # Percent below in smaller text
  geom_text(aes(label = paste0("(", percent, "%)")), color = "black", size = 3.5, vjust = 1.5) +
  scale_fill_gradient(low = "lightyellow", high = "red") +
  labs(
    x = "MAVE Classification",
    y = "ClinVar Classification",
    fill = "Percentage"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),      # bigger tick labels
    axis.title.x = element_text(margin = margin(t = 20)),  # space above x-axis title
    axis.title.y = element_text(margin = margin(r = 20))   # space right of y-axis title
  )

