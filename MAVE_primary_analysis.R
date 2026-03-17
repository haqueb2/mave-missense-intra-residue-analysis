# Load libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(httr)
library(jsonlite)
library(purrr)

# Read in input dataset
maveData = read_excel("inputs/all_gene_missense_only_01052026.xlsx")
geneList = unique(maveData$gene) # 24 unique genes 

# Normalize protein_variant column 
aa1_to_aa3 = c(
  A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", 
  E = "Glu", Q = "Gln", G = "Gly", H = "His", I = "Ile", 
  L = "Leu", K = "Lys", M = "Met", F = "Phe", P = "Pro", 
  S = "Ser", T = "Thr", W = "Trp", Y = "Tyr", V = "Val"
)

# Remove rows with B or Z in protein_variant
maveData = maveData %>%
  filter(!str_detect(protein_variant, "B"),
         !str_detect(protein_variant, "Z"))

# Make function
convert_variant = function(variant) {
  if (grepl("^p\\.[A-Z]\\d+[A-Z]$", variant)) {
    # Extract reference, position, and alt amino acids
    ref = sub("^p\\.([A-Z]).*", "\\1", variant)
    pos = sub("^p\\.[A-Z](\\d+)[A-Z]$", "\\1", variant)
    alt = sub("^p\\.[A-Z]\\d+([A-Z])$", "\\1", variant)
    
    # Convert to 3-letter codes
    paste0("p.", aa1_to_aa3[[ref]], pos, aa1_to_aa3[[alt]])
  } else {
    # Already 3-letter or special format, leave unchanged
    variant
  }
}

# Apply function
maveData = maveData %>%
  rowwise() %>%
  mutate(protein_variant = convert_variant(protein_variant)) %>%
  ungroup()

# Identify any issues in protein_variant column
anomalies = maveData %>%
  filter(!grepl("^p\\.[A-Z][a-z]{2}\\d+[A-Z][a-z]{2}$", protein_variant))

# Create column for amino acid position
maveData = maveData %>%
  mutate(aa_pos = as.integer(sub("^p\\.[A-Z][a-z]{2}(\\d+)[A-Z][a-z]{2}$", "\\1", protein_variant)))

# Remove transcript_variant and score columns 
maveData_new = subset(maveData, select = -c(3,5))
maveData_new = unique(maveData_new)

# Find rows where gene & protein_variant are the same (duplicates)
duplicates = maveData_new %>%
  group_by(gene, protein_variant) %>%
  filter(n() > 1)

# Delete ALL rows where gene + protein_variant are duplicated
maveData_new = maveData_new %>%
  add_count(gene, protein_variant) %>%
  filter(n == 1) %>%   # keep only rows that appear once
  select(-n)           # drop the helper column

# Exclude INT variants 
maveData_new = maveData_new[maveData_new$func_classification != "INT",]

# Re-name
maveData = maveData_new

# Create column for variant count per residue & per gene
maveData = maveData %>%
  group_by(gene, aa_pos) %>%
  mutate(variant_count_per_residue = n()) %>%
  ungroup()

## PRIMARY ANALYSIS ##
# Group by residue position and only consider residues with >1 variant
filteredData = maveData[maveData$variant_count_per_residue>1, ]
filteredData$variant_count_per_residue = as.numeric(filteredData$variant_count_per_residue)

# Summarize data by residue 
varAnalysis = filteredData %>%
  group_by(gene, aa_pos) %>%
  summarise(
    # Get the variants and count information for each residue position
    protein_variants = paste(unique(protein_variant), collapse = ","),
    variant_count_per_residue = paste(unique(variant_count_per_residue), collapse = ","),
    
    # Collect gene_type and gene_only
    gene_type = paste(unique(gene_type), collapse = ","),
    gene_only = paste(unique(gene_only), collapse = ","),
    
    # Class counts
    FUNC_count = sum(func_classification == "FUNC", na.rm = TRUE),
    LOF_count = sum(func_classification == "LOF", na.rm = TRUE),
  ) 

# Create LOF_proportion column
varAnalysis$variant_count_per_residue = as.numeric(varAnalysis$variant_count_per_residue)
varAnalysis = varAnalysis %>%
  mutate(LOF_proportion = LOF_count / variant_count_per_residue)

# Create LOF_likelihood column
varAnalysis = varAnalysis %>%
  mutate(
    LOF_likelihood = ifelse(
      LOF_count >= 1 & variant_count_per_residue > 1,
      (LOF_count - 1) / (variant_count_per_residue - 1),
      NA_real_
    )
  )

# Create LR_pathogenicity column 
varAnalysis = varAnalysis %>%
  mutate(
    LR_pathogenicity = ifelse(
      (FUNC_count) > 0,
      as.character(LOF_count / (FUNC_count)),
      "Undefined"
    )
  )

# Create evidence_strength column 
varAnalysis = varAnalysis %>%
  mutate(
    evidence_strength = case_when(
      LR_pathogenicity == "Undefined" ~ "Very Strong",
      as.numeric(LR_pathogenicity) >= 350 ~ "Very Strong",
      as.numeric(LR_pathogenicity) >= 18.7 ~ "Strong",
      as.numeric(LR_pathogenicity) >= 4.33 ~ "Moderate",
      as.numeric(LR_pathogenicity) >= 2.08 ~ "Supporting",
      TRUE ~ NA_character_  
    )
  )

## INTERPRO ANNOTATION ##
# STEP 1: Get Uniprot protein ID from gene name 
idMap = read.table("annotations/HUMAN_9606_idmapping.dat", header = F, sep = "\t", quote = "", fill = TRUE)

# Filter idMap to keep only rows with Gene_Name
geneIdMap = idMap[idMap$V2 == "Gene_Name",]

# Create a new column 'uniprotID' in varAnalysis and initialize with NA
varAnalysis$uniprotID = NA

# Iterate through unique gene names in varAnalysis
for (gene in unique(varAnalysis$gene_only)) {
  # Check if the gene is present in geneIdMap
  if (gene %in% geneIdMap$V3) {
    # Extract the corresponding Uniprot IDs
    uniprotIDs = geneIdMap[geneIdMap$V3 == gene, "V1"]
    
    # Combine multiple Uniprot IDs into a single string separated by a semicolon
    uniprotID = paste(uniprotIDs, collapse = ";")
    
    # Update the uniprotID column in varAnalysis for all rows with the same gene
    varAnalysis$uniprotID[varAnalysis$gene_only %in% gene] = uniprotID
  }
}

# Only keep the first item in the semi-colon separated list of varAnalysis$uniprotID
varAnalysis$uniprotID = sapply(strsplit(as.character(varAnalysis$uniprotID), ";"), function(x) x[1])

# STEP 2: Annotate with InterPro
# Get the list of unique Uniprot IDs
proteinList = unique(varAnalysis$uniprotID)
combined_results = data.frame()  # Initialize an empty dataframe to store the combined results

for (protein_id in proteinList) {
  # Fetch the data from the API for the current protein ID
  response = GET(paste0("https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/reviewed/", protein_id))
  data = content(response, "text")
  
  # Skip to the next protein_id if data is empty
  if (data == "") {
    warning(paste("No data returned for protein ID:", protein_id))
    next
  }
  
  json_data = fromJSON(data)
  
  # Extract metadata and proteins
  metadata = json_data$results$metadata
  proteins = json_data$results$proteins
  
  # Skip to the next protein_id if metadata or proteins are NULL
  if (is.null(metadata) || is.null(proteins)) {
    warning(paste("No metadata or proteins data for protein ID:", protein_id))
    next
  }
  
  # Extract the entry_protein_locations data along with the associated accession IDs from metadata
  for (i in seq_along(proteins)) {
    protein = proteins[[i]]
    accession_id = metadata$accession[i]
    name = metadata$name[i]
    type = metadata$type[i]
    entry_protein_locations = protein$entry_protein_locations %>%
      map_df(~ as.data.frame(.)) %>%
      unnest(cols = fragments, names_sep = "_")
    
    if (nrow(entry_protein_locations) > 0) {
      # Add new columns for protein_id and accession_id
      entry_protein_locations = entry_protein_locations %>%
        mutate(protein_id = protein_id, accession_id = accession_id, name = name, type = type)
      
      # Append the results to the combined dataframe
      combined_results = bind_rows(combined_results, entry_protein_locations)
    }
  }
}

# Remove columns that are entirely NA
combined_results = combined_results %>% select_if(~ !all(is.na(.)))
combined_results$name = paste0(combined_results$name, " (", combined_results$type, ")")

# Create Interpro column 
varAnalysis$interpro = NA

for (i in 1:nrow(combined_results)) {
  # Get the start and end values from combined_results
  fragments_start = combined_results$fragments_start[i]
  fragments_end = combined_results$fragments_end[i]
  protein_id = combined_results$protein_id[i]
  
  # Find rows in varAnalysis that overlap with the current region
  overlap_rows = which(
    as.numeric(varAnalysis$aa_pos) <= fragments_end & 
      as.numeric(varAnalysis$aa_pos) >= fragments_start &
      varAnalysis$uniprotID == protein_id
  )
  
  # Populate the corresponding column in varAnalysis
  if (length(overlap_rows) > 0) {
    # Combine descriptions if multiple overlaps
    varAnalysis$interpro[overlap_rows] = ifelse(
      is.na(varAnalysis$interpro[overlap_rows]),
      combined_results$name[i],
      paste(varAnalysis$interpro[overlap_rows], combined_results$name[i], sep = ";")
    )
  }
}

# Extract families
# All known tags
tags = c("domain", "family", "homologous_superfamily",
          "binding_site", "active_site", "conserved_site", "ptm", "repeat")

# Create empty columns
for (tag in tags) {
  varAnalysis[[paste0("interpro_", tag)]] = NA
}

# Fill each column by extracting annotations of that type
for (tag in tags) {
  pattern = paste0("[^;]+\\(", tag, "\\)")
  
  varAnalysis[[paste0("interpro_", tag)]] = sapply(varAnalysis$interpro, function(x) {
    if (is.na(x)) return(NA)
    matches = str_extract_all(x, pattern)[[1]]
    if (length(matches) == 0) return(NA)
    paste(matches, collapse = ";")
  })
}

## MAKING PLOTS ## 
# 1. HEATMAP OF LOF_PROPORTION ALONG GENES
ggplot(varAnalysis, aes(
  x = aa_pos, 
  y = reorder(gene, aa_pos, FUN = median),
  fill = LOF_proportion
)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Amino Acid Residue Position", y = "Gene") + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),            # removes grid lines
    axis.text.y = element_text(face = "italic")  # italicizes y-axis labels
  )

# 2. SCATTERPLOT OF LOF_LIKELIHOOD VS. LOF_PROPORTION 
ggplot(varAnalysis, aes(x=LOF_proportion, y=LOF_likelihood, color=gene)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs (x = "Proportion of LOF Variants", y = "Likelihood of LOF") +
  theme_minimal() 

#3. PLOTTING EVIDENCE STRENGTH VS LOF PROPORTION
mavePlot = varAnalysis

# Ensure LR_pathogenicity is character before replacement
mavePlot$LR_pathogenicity = as.character(mavePlot$LR_pathogenicity)
mavePlot$LR_pathogenicity[mavePlot$LR_pathogenicity == "Undefined"] = "50"
mavePlot$LR_pathogenicity = as.numeric(mavePlot$LR_pathogenicity)

# Replace NA evidence_strength with "No pathogenic evidence strength"
mavePlot$evidence_strength = as.character(mavePlot$evidence_strength)
mavePlot$evidence_strength[is.na(mavePlot$evidence_strength)] = "None"
mavePlot$evidence_strength = factor(
  mavePlot$evidence_strength,
  levels = c("Very Strong", "Strong", "Moderate", "Supporting", "None")
)

# Find first LOF_proportion where evidence_strength is "Supporting"
supporting_threshold = min(mavePlot$LOF_proportion[mavePlot$evidence_strength == "Supporting"], na.rm = TRUE)

# Plot: LOF_proportion vs LR_pathogenicity
ggplot(mavePlot, aes(x = LOF_proportion, y = LR_pathogenicity, color = evidence_strength)) +
  geom_point() +
  {if (!is.na(supporting_threshold)) geom_vline(xintercept = supporting_threshold, linetype = "dashed", color = "gray")} +
  labs(x = "Proportion of LOF Variants", y = "Likelihood Ratio of Pathogenicity") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(
    breaks = c(1, 2, 5, 10, 20, 50),
    labels = function(x) ifelse(x == 50, "Undefined", x)
  ) +
  guides(color = guide_legend(title = "Evidence Strength", override.aes = list(alpha = 1)))

#4. PLOTTING EVIDENECE STRENGTH VS LOF LIKLIHOOD
# Find first LOF_likelihood where evidence_strength is "Supporting"
supporting_threshold = min(mavePlot$LOF_likelihood[mavePlot$evidence_strength == "Supporting"], na.rm = TRUE)

# Scatter plot: LOF_likelihood vs LR_pathogenicity
ggplot(mavePlot, aes(x = LOF_likelihood, y = LR_pathogenicity, color = evidence_strength)) +
  geom_point() +
  geom_vline(xintercept = supporting_threshold, linetype = "dashed", color = "gray") +
  labs(
    x = "Likelihood of LOF",
    y = "Likelihood Ratio of Pathogenicity"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(
    breaks = c(1, 2, 5, 10, 15, 20, 25, 50),
    labels = function(x) ifelse(x == 50, "Undefined", x)
  ) +
  guides(color = guide_legend(title = "Pathogenic Evidence Strength", override.aes = list(alpha = 1)))

# 4.1 PLOTTING EVIDENCE STRENGTH VS LOF LIKELIHOOD AT LOG SCALE (with jitter)
ggplot(mavePlot, aes(x = LOF_likelihood, y = LR_pathogenicity, color = evidence_strength)) +
  geom_count(alpha = 0.7) +   # count-based point size
  geom_vline(xintercept = supporting_threshold, linetype = "dashed", color = "gray") +
  labs(
    x = "Likelihood of LOF",
    y = "Likelihood Ratio of Pathogenicity (log10)"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  scale_y_log10(
    breaks = scales::log_breaks(n = 10), 
    labels = function(x) ifelse(x == 50, "Undefined", x)
  ) +
  guides(
    size = guide_legend(title = "Count"),  # show count in legend
    color = guide_legend(title = "Evidence Strength", override.aes = list(alpha = 1))
  )
