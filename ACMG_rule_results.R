# Set working directory
setwd("/Volumes/T9/MAVE Data Project")

# Load libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Read in data that uses CLINVAR
maveData = read.csv("new_acmg/new_ACMG_scoring_results.csv")

# Re-format data 
filteredData = maveData[maveData$canonical == "TRUE",] # canonical transcripts
filteredData = filteredData[filteredData$in_clinvar == "FALSE", ] # only variants that are not reported in clinvar

# Cap new_ACMG_score to range -4 to 8
filteredData$new_ACMG_score = pmax(pmin(filteredData$new_ACMG_score, 8), -4)

# Remove rows where new_ACMG_score is NA
filteredData = filteredData[!is.na(filteredData$new_ACMG_score), ]

# Plot new_ACMG_score and func_classification
ggplot(filteredData, aes(x = factor(new_ACMG_score), fill = func_classification)) +
  geom_bar(position = "fill") +  # Stacked
  labs(
    x = "New total ACMG score (IMP_MSS)",
    y = "Proportion of MAVE missense variants",
    fill = "MAVE functional classification",
  ) +
  theme_minimal() +
  theme(
   panel.grid = element_blank(),
  axis.text = element_text(size = 10)
  )

# Remove rows where diff_aa_grantham_dis is NA
diffAA_plot_data = filteredData[!is.na(filteredData$diff_aa_grantham_dis), ]
  
# Plot diff_aa_grantham_dis and func_classification
ggplot(diffAA_plot_data, aes(x = factor(diff_aa_grantham_dis), fill = func_classification)) +
  geom_bar(position = "fill") +  # Stacked
  labs(
    x = "New PM5 ACMG score",
    y = "Proportion of MAVE missense variants",
    fill = "MAVE functional classification",
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10)
  )

# PLOT NEW ACMG SCORES FOR EACH PREDICTOR
# Select the columns you need
plot_data = filteredData %>%
  select(REVEL_score, revel_eval, REVEL_pred,
         BayesDel_noAF_score, bayesdel_eval, BayesDel_noAF_pred,
         MutPred_score, mutpred_eval, MutPred_pred,
         VEST4_score, vest_eval, VEST4_pred)

# Reshape to long format
plot_long = plot_data %>%
  pivot_longer(
    cols = c(revel_eval, bayesdel_eval, mutpred_eval, vest_eval),
    names_to = "Predictor",
    values_to = "Integer_Score"
  ) %>%
  mutate(
    Predictor = factor(Predictor, 
                       levels = c("revel_eval", "bayesdel_eval", "mutpred_eval", "vest_eval"),
                       labels = c("REVEL", "BayesDel", "MutPred2", "VEST4")),
    # Add continuous scores
    Continuous_Score = case_when(
      Predictor == "REVEL" ~ REVEL_score,
      Predictor == "BayesDel" ~ BayesDel_noAF_score,
      Predictor == "MutPred2" ~ MutPred_score,
      Predictor == "VEST4" ~ VEST4_score
    ),
    # Add prediction label
    Pred_Label = case_when(
      Predictor == "REVEL" ~ REVEL_pred,
      Predictor == "BayesDel" ~ BayesDel_noAF_pred,
      Predictor == "MutPred2" ~ MutPred_pred,
      Predictor == "VEST4" ~ VEST4_pred
    )
  ) %>%
  filter(!is.na(Integer_Score))

# Plot boxplot, color-coded by Pred_Label
ggplot(plot_long, aes(x = factor(Integer_Score), y = Continuous_Score, fill = Pred_Label)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.6, color = "black") +
  facet_wrap(~ Predictor, scales = "free_y") +
  scale_fill_manual(values = c("T" = "#4daf4a", "D" = "#e41a1c")) +  # green = T, red = D
  labs(
    x = "New VEP ACMG score",
    y = "In silico VEP score",
    fill = "VEP Label"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

# Reshape data, calculate counts + percentages, and fill missing combinations with 0
heatmap_data = filteredData %>%
  select(revel_eval, bayesdel_eval, mutpred_eval, vest_eval, diff_aa_grantham_dis) %>%
  pivot_longer(
    cols = c(revel_eval, bayesdel_eval, mutpred_eval, vest_eval),
    names_to = "Predictor",
    values_to = "Eval_Score"
  ) %>%
  mutate(
    # Group predictor-specific eval scores
    eval_group = case_when(
      Eval_Score < 0 ~ "Negative",
      Eval_Score == 0 ~ "Zero",
      Eval_Score > 0 ~ "Positive",
      TRUE ~ NA_character_
    ),
    # Group grantham difference
    grantham_group = case_when(
      diff_aa_grantham_dis < 0 ~ "Negative",
      diff_aa_grantham_dis == 0 ~ "Zero",
      diff_aa_grantham_dis > 0 ~ "Positive",
      TRUE ~ NA_character_
    ),
    # Nice predictor labels
    Predictor = factor(Predictor,
                       levels = c("revel_eval", "bayesdel_eval", "mutpred_eval", "vest_eval"),
                       labels = c("REVEL", "BayesDel", "MutPred2", "VEST4")),
    # Factor groups so they appear Negative > Zero > Positive
    eval_group = factor(eval_group, levels = c("Negative", "Zero", "Positive")),
    grantham_group = factor(grantham_group, levels = c("Negative", "Zero", "Positive"))
  ) %>%
  filter(!is.na(eval_group), !is.na(grantham_group)) %>%
  count(Predictor, grantham_group, eval_group) %>%
  # Fill in missing combinations with 0
  complete(Predictor, grantham_group, eval_group, fill = list(n = 0)) %>%
  group_by(Predictor) %>%
  mutate(
    pct = n / sum(n) * 100,
    label_count = n,
    label_pct = paste0("(", round(pct, 1), "%)")
  ) %>%
  ungroup()

# Plot faceted heatmap
ggplot(heatmap_data, aes(x = eval_group, y = grantham_group, fill = pct)) +
  geom_tile(color = "white") +
  # First line: counts
  geom_text(aes(label = label_count), color = "black", size = 4, vjust = 0) +
  # Second line: percentages, pushed down a bit
  geom_text(aes(label = label_pct), color = "black", size = 3.5, vjust = 1.5) +
  scale_fill_gradient(low = "#f0f0f0", high = "#08519c") +
  facet_wrap(~ Predictor) +
  labs(
    x = "New VEP ACMG score (grouped)",
    y = "New PM5 ACMG score (grouped)",
    fill = "Percentage"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

# Read in data that uses MAVE
maveData = read.csv("new_acmg/new_ACMG_scoring_MAVE_modified.csv") 

# Re-format data 
filteredData = maveData[maveData$canonical == "TRUE",] # canonical transcripts

# Cap new_ACMG_score_mod to range -4 to 8
filteredData$new_ACMG_score_mod = pmax(pmin(filteredData$new_ACMG_score_mod, 8), -4)

# Remove rows where new_ACMG_score_mod is NA
filteredData = filteredData[!is.na(filteredData$new_ACMG_score_mod), ]

# Plot new_ACMG_score_mod and func_classification
ggplot(filteredData, aes(x = factor(new_ACMG_score_mod), fill = func_classification)) +
  geom_bar(position = "fill") +  # Stacked
  labs(
    x = "New total ACMG score (IMP_MSS)",
    y = "Proportion of MAVE missense variants",
    fill = "MAVE functional classification",
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10)
  )

# Remove rows where diff_aa_grantham_dis is NA
diffAA_plot_data = filteredData[!is.na(filteredData$diff_aa_grantham_dis), ]

# Plot diff_aa_grantham_dis and func_classification
ggplot(diffAA_plot_data, aes(x = factor(diff_aa_grantham_dis), fill = func_classification)) +
  geom_bar(position = "fill") +  # Stacked
  labs(
    x = "New PM5 ACMG score",
    y = "Proportion of MAVE missense variants",
    fill = "MAVE functional classification",
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10)
  )

# PLOT NEW ACMG SCORES FOR EACH PREDICTOR
# Select the columns you need
plot_data = filteredData %>%
  select(REVEL_score, revel_eval, REVEL_pred,
         BayesDel_noAF_score, bayesdel_eval, BayesDel_noAF_pred,
         MutPred_score, mutpred_eval, MutPred_pred,
         VEST4_score, vest_eval, VEST4_pred)

# Reshape to long format
plot_long = plot_data %>%
  pivot_longer(
    cols = c(revel_eval, bayesdel_eval, mutpred_eval, vest_eval),
    names_to = "Predictor",
    values_to = "Integer_Score"
  ) %>%
  mutate(
    Predictor = factor(Predictor, 
                       levels = c("revel_eval", "bayesdel_eval", "mutpred_eval", "vest_eval"),
                       labels = c("REVEL", "BayesDel", "MutPred2", "VEST4")),
    # Add continuous scores
    Continuous_Score = case_when(
      Predictor == "REVEL" ~ REVEL_score,
      Predictor == "BayesDel" ~ BayesDel_noAF_score,
      Predictor == "MutPred2" ~ MutPred_score,
      Predictor == "VEST4" ~ VEST4_score
    ),
    # Add prediction label
    Pred_Label = case_when(
      Predictor == "REVEL" ~ REVEL_pred,
      Predictor == "BayesDel" ~ BayesDel_noAF_pred,
      Predictor == "MutPred2" ~ MutPred_pred,
      Predictor == "VEST4" ~ VEST4_pred
    )
  ) %>%
  filter(!is.na(Integer_Score))

# Plot boxplot, color-coded by Pred_Label
ggplot(plot_long, aes(x = factor(Integer_Score), y = Continuous_Score, fill = Pred_Label)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.6, color = "black") +
  facet_wrap(~ Predictor, scales = "free_y") +
  scale_fill_manual(values = c("T" = "#4daf4a", "D" = "#e41a1c")) +  # green = T, red = D
  labs(
    x = "New VEP ACMG score",
    y = "In silico VEP score",
    fill = "VEP Label"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

# Reshape data, calculate counts + percentages, and fill missing combinations with 0
heatmap_data = filteredData %>%
  select(revel_eval, bayesdel_eval, mutpred_eval, vest_eval, diff_aa_grantham_dis) %>%
  pivot_longer(
    cols = c(revel_eval, bayesdel_eval, mutpred_eval, vest_eval),
    names_to = "Predictor",
    values_to = "Eval_Score"
  ) %>%
  mutate(
    # Group predictor-specific eval scores
    eval_group = case_when(
      Eval_Score < 0 ~ "Negative",
      Eval_Score == 0 ~ "Zero",
      Eval_Score > 0 ~ "Positive",
      TRUE ~ NA_character_
    ),
    # Group grantham difference
    grantham_group = case_when(
      diff_aa_grantham_dis < 0 ~ "Negative",
      diff_aa_grantham_dis == 0 ~ "Zero",
      diff_aa_grantham_dis > 0 ~ "Positive",
      TRUE ~ NA_character_
    ),
    # Nice predictor labels
    Predictor = factor(Predictor,
                       levels = c("revel_eval", "bayesdel_eval", "mutpred_eval", "vest_eval"),
                       labels = c("REVEL", "BayesDel", "MutPred2", "VEST4")),
    # Factor groups so they appear Negative > Zero > Positive
    eval_group = factor(eval_group, levels = c("Negative", "Zero", "Positive")),
    grantham_group = factor(grantham_group, levels = c("Negative", "Zero", "Positive"))
  ) %>%
  filter(!is.na(eval_group), !is.na(grantham_group)) %>%
  count(Predictor, grantham_group, eval_group) %>%
  # Fill in missing combinations with 0
  complete(Predictor, grantham_group, eval_group, fill = list(n = 0)) %>%
  group_by(Predictor) %>%
  mutate(
    pct = n / sum(n) * 100,
    label_count = n,
    label_pct = paste0("(", round(pct, 1), "%)")
  ) %>%
  ungroup()

# Plot faceted heatmap
ggplot(heatmap_data, aes(x = eval_group, y = grantham_group, fill = pct)) +
  geom_tile(color = "white") +
  # First line: counts
  geom_text(aes(label = label_count), color = "black", size = 4, vjust = 0) +
  # Second line: percentages, pushed down a bit
  geom_text(aes(label = label_pct), color = "black", size = 3.5, vjust = 1.5) +
  scale_fill_gradient(low = "#f0f0f0", high = "#08519c") +
  facet_wrap(~ Predictor) +
  labs(
    x = "New VEP ACMG score (grouped)",
    y = "New PM5 ACMG score (grouped)",
    fill = "Percentage"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

# Read in data that uses MAVE (V2)
maveData = read.csv("") 

# Re-format data 
filteredData = maveData[maveData$canonical == "TRUE",] # canonical transcripts

# Cap new_ACMG_score_mod to range -4 to 8
filteredData$new_ACMG_score_mod = pmax(pmin(filteredData$new_ACMG_score_mod, 8), -4)

# Remove rows where new_ACMG_score_mod is NA
filteredData = filteredData[!is.na(filteredData$new_ACMG_score_mod), ]

# Plot new_ACMG_score_mod and func_classification
ggplot(filteredData, aes(x = factor(new_ACMG_score_mod), fill = func_classification)) +
  geom_bar(position = "fill") +  # Stacked
  labs(
    x = "New total ACMG score (IMP_MSS)",
    y = "Proportion of MAVE missense variants",
    fill = "MAVE functional classification",
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10)
  )

# Remove rows where diff_aa_grantham_dis is NA
diffAA_plot_data = filteredData[!is.na(filteredData$diff_aa_grantham_dis), ]

# Plot diff_aa_grantham_dis and func_classification
ggplot(diffAA_plot_data, aes(x = factor(diff_aa_grantham_dis), fill = func_classification)) +
  geom_bar(position = "fill") +  # Stacked
  labs(
    x = "New PM5 ACMG score",
    y = "Proportion of MAVE missense variants",
    fill = "MAVE functional classification",
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10)
  )

# PLOT NEW ACMG SCORES FOR EACH PREDICTOR
# Select the columns you need
plot_data = filteredData %>%
  select(REVEL_score, revel_eval, REVEL_pred,
         BayesDel_noAF_score, bayesdel_eval, BayesDel_noAF_pred,
         MutPred_score, mutpred_eval, MutPred_pred,
         VEST4_score, vest_eval, VEST4_pred)

# Reshape to long format
plot_long = plot_data %>%
  pivot_longer(
    cols = c(revel_eval, bayesdel_eval, mutpred_eval, vest_eval),
    names_to = "Predictor",
    values_to = "Integer_Score"
  ) %>%
  mutate(
    Predictor = factor(Predictor, 
                       levels = c("revel_eval", "bayesdel_eval", "mutpred_eval", "vest_eval"),
                       labels = c("REVEL", "BayesDel", "MutPred2", "VEST4")),
    # Add continuous scores
    Continuous_Score = case_when(
      Predictor == "REVEL" ~ REVEL_score,
      Predictor == "BayesDel" ~ BayesDel_noAF_score,
      Predictor == "MutPred2" ~ MutPred_score,
      Predictor == "VEST4" ~ VEST4_score
    ),
    # Add prediction label
    Pred_Label = case_when(
      Predictor == "REVEL" ~ REVEL_pred,
      Predictor == "BayesDel" ~ BayesDel_noAF_pred,
      Predictor == "MutPred2" ~ MutPred_pred,
      Predictor == "VEST4" ~ VEST4_pred
    )
  ) %>%
  filter(!is.na(Integer_Score))

# Plot boxplot, color-coded by Pred_Label
ggplot(plot_long, aes(x = factor(Integer_Score), y = Continuous_Score, fill = Pred_Label)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.6, color = "black") +
  facet_wrap(~ Predictor, scales = "free_y") +
  scale_fill_manual(values = c("T" = "#4daf4a", "D" = "#e41a1c")) +  # green = T, red = D
  labs(
    x = "New VEP ACMG score",
    y = "In silico VEP score",
    fill = "VEP Label"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

# Reshape data, calculate counts + percentages, and fill missing combinations with 0
heatmap_data = filteredData %>%
  select(revel_eval, bayesdel_eval, mutpred_eval, vest_eval, diff_aa_grantham_dis) %>%
  pivot_longer(
    cols = c(revel_eval, bayesdel_eval, mutpred_eval, vest_eval),
    names_to = "Predictor",
    values_to = "Eval_Score"
  ) %>%
  mutate(
    # Group predictor-specific eval scores
    eval_group = case_when(
      Eval_Score < 0 ~ "Negative",
      Eval_Score == 0 ~ "Zero",
      Eval_Score > 0 ~ "Positive",
      TRUE ~ NA_character_
    ),
    # Group grantham difference
    grantham_group = case_when(
      diff_aa_grantham_dis < 0 ~ "Negative",
      diff_aa_grantham_dis == 0 ~ "Zero",
      diff_aa_grantham_dis > 0 ~ "Positive",
      TRUE ~ NA_character_
    ),
    # Nice predictor labels
    Predictor = factor(Predictor,
                       levels = c("revel_eval", "bayesdel_eval", "mutpred_eval", "vest_eval"),
                       labels = c("REVEL", "BayesDel", "MutPred2", "VEST4")),
    # Factor groups so they appear Negative > Zero > Positive
    eval_group = factor(eval_group, levels = c("Negative", "Zero", "Positive")),
    grantham_group = factor(grantham_group, levels = c("Negative", "Zero", "Positive"))
  ) %>%
  filter(!is.na(eval_group), !is.na(grantham_group)) %>%
  count(Predictor, grantham_group, eval_group) %>%
  # Fill in missing combinations with 0
  complete(Predictor, grantham_group, eval_group, fill = list(n = 0)) %>%
  group_by(Predictor) %>%
  mutate(
    pct = n / sum(n) * 100,
    label_count = n,
    label_pct = paste0("(", round(pct, 1), "%)")
  ) %>%
  ungroup()

# Plot faceted heatmap
ggplot(heatmap_data, aes(x = eval_group, y = grantham_group, fill = pct)) +
  geom_tile(color = "white") +
  # First line: counts
  geom_text(aes(label = label_count), color = "black", size = 4, vjust = 0) +
  # Second line: percentages, pushed down a bit
  geom_text(aes(label = label_pct), color = "black", size = 3.5, vjust = 1.5) +
  scale_fill_gradient(low = "#f0f0f0", high = "#08519c") +
  facet_wrap(~ Predictor) +
  labs(
    x = "New VEP ACMG score (grouped)",
    y = "New PM5 ACMG score (grouped)",
    fill = "Percentage"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

