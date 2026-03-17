# Load libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(httr)
library(jsonlite)

# Read in input dataset
maveData = read_excel("inputs/all_gene_missense_only_01052026.xlsx")
geneList = unique(maveData$gene) # 21 unique genes 

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
    ref <- sub("^p\\.([A-Z]).*", "\\1", variant)
    pos <- sub("^p\\.[A-Z](\\d+)[A-Z]$", "\\1", variant)
    alt <- sub("^p\\.[A-Z]\\d+([A-Z])$", "\\1", variant)
    
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

# Re-name
maveData = maveData_new

# Create column with p. format
maveData$p_format = NA
maveData$p_format = paste0(maveData$gene_only, ":", maveData$protein_variant)

# Create column for variant count per residue & per gene
maveData = maveData %>%
  group_by(gene, aa_pos) %>%
  mutate(variant_count_per_residue = n()) %>%
  ungroup()

## ANNOTATE WITH GRANTHAM DISTANCE ## 
# Load in Grantham matrix
grantham_data = read.delim("annotations/grantham.matrix", header = T, row.names = 1)
grantham_matrix = as.matrix(grantham_data)

# Lookup table from 3-letter amino acid letters to 1-letter
aa3_to_aa1 = c(
  Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
  Glu = "E", Gln = "Q", Gly = "G", His = "H", Ile = "I",
  Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
  Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V"
)

# Create 1-letter amino acid columns
maveData = maveData %>%
  mutate(
    ref_aa3 = sub("^p\\.([A-Z][a-z]{2}).*$", "\\1", protein_variant),
    alt_aa3 = sub("^p\\.[A-Z][a-z]{2}\\d+([A-Z][a-z]{2})$", "\\1", protein_variant),
    ref_aa1 = aa3_to_aa1[ref_aa3],
    alt_aa1 = aa3_to_aa1[alt_aa3]
  )

# Add Grantham distance annotations
maveData = maveData %>%
  rowwise() %>%
  mutate(grantham_distance = grantham_matrix[ref_aa1, alt_aa1]) %>%
  ungroup()

## FIGURE OUT GENOMIC COORDINATES ##
# Path to VEP_out directory
vepDir = "/VEP_out"

# List all .txt files in VEP_out folder
txtFiles = list.files(vepDir, pattern = "\\.txt$", full.names = TRUE)

# Read in files into a list of data frames
vepDataList = lapply(txtFiles, read.delim, stringsAsFactors = FALSE)

# Only keep the first 32 columns in each dataframe 
vepDataList = lapply(vepDataList, function(df) df[, seq_len(min(32, ncol(df)))])

# Combine into one dataframe
vepDataCombined = do.call(rbind, vepDataList)

# Redefine gene list
geneList = unique(maveData$gene_only)

# Apply filters to data frame 
vepDataFilter = vepDataCombined[vepDataCombined$SYMBOL %in% geneList, ]

# Work on a copy
uv = vepDataFilter$X.Uploaded_variation

# Identify rows with 1-letter AA protein notation (e.g. p.G1820V)
is_one_letter = grepl(":p\\.[A-Z][0-9]+[A-Z]$", uv)

# Split only those rows
gene = sub(":p\\..*", "", uv[is_one_letter])
prot = sub(".*:p\\.", "", uv[is_one_letter])

parts = do.call(
  rbind,
  strsplit(prot, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)", perl = TRUE)
)

# Convert to 3-letter AA
parts[,1] = aa1_to_aa3[parts[,1]]
parts[,3] = aa1_to_aa3[parts[,3]]

# Recombine and replace ONLY those rows
uv[is_one_letter] = paste0(
  gene, ":p.", parts[,1], parts[,2], parts[,3]
)

# Put back into dataframe
vepDataFilter$X.Uploaded_variation <- uv

# Make sure the look-up columns have compatible types
vepDataFilter = vepDataFilter %>%
  mutate(Uploaded_variation = as.character(X.Uploaded_variation),
         Protein_position = as.numeric(Protein_position))

# Make sure the look-up columns have compatible types
maveData = maveData %>%
  mutate(p_format = as.character(p_format),
         aa_pos = as.numeric(aa_pos))

# Filter vepDataFilter by keeping only rows with matching pairs in maveData
vepDataFilter = vepDataFilter %>%
  semi_join(maveData,
            by = c("Uploaded_variation" = "p_format",
                   "Protein_position" = "aa_pos"))

vepDataFilter = unique(vepDataFilter)

# Modify vepDataFilter
# 1. Parse Location: chr, start, stop
vepDataFilter = vepDataFilter %>%
  separate(Location, into = c("chr", "coords"), sep = ":") %>%
  separate(coords, into = c("start", "stop"), sep = "-") %>%
  mutate(
    start = as.numeric(start),
    stop  = as.numeric(stop)
  )

# 2. Parse UPLOADED_ALLELE: ref, alt
vepDataFilter = vepDataFilter %>%
  separate(UPLOADED_ALLELE, into = c("ref", "alt"), sep = "/")

# 3. Join with maveData

# Prioritize MANE_Select and keep only one row per Uploaded_variation
vep_prioritized = vepDataFilter %>%
  group_by(Uploaded_variation) %>%
  slice_max(order_by = (MANE == "MANE_Select"), n = 1, with_ties = FALSE) %>%
  ungroup()

# Apply join
maveData_annotated = maveData %>%
  select(-any_of(c("chr", "start", "stop", "ref", "alt"))) %>%
  left_join(
    vep_prioritized %>% 
      select(Uploaded_variation, chr, start, stop, ref, alt, MANE_SELECT),
    by = c("p_format" = "Uploaded_variation")
  )

maveData_annotated = unique(maveData_annotated) # unique rows only

## SKIP IF RUN BEFORE ONCE, READ IN RESULTS INSTEAD ##
# 4. For missing information, use TransVar API 
missing_rows = maveData_annotated %>%
  filter(is.na(chr)) %>%
  select(p_format)

# Initialize empty dataframe before the loop
all_results = data.frame()

for (variant in missing_rows$p_format) {
  # Prepare POST body form-data
  body = list(
    app = "transvar",
    command_id = "transvar",
    task = "panno",
    refversion = "hg38",
    refseq = "refseq",
    typed_identifiers = variant
  )
  
  # POST request to TransVar API
  response = POST(
    url = "https://bioinformatics.mdanderson.org/transvar_dyce/",
    encode = "multipart",
    body = body
  )
  
  # Parse JSON content
  content_json = content(response, as = "text", encoding = "UTF-8")
  content_list = fromJSON(content_json)
  
  if (!is.null(content_list$results) && grepl("^https?://", content_list$results)) {
    # Download the results file (TSV)
    results_url = content_list$results
    results_txt = read.delim(results_url, stringsAsFactors = FALSE)
    
    # Add a column to identify which variant this came from
    results_txt$p_format = variant
    
    # Append to all_results
    all_results = rbind(all_results, results_txt)
  } else {
    # If no valid results, create a dummy row with NAs and append
    na_row = data.frame(matrix(NA, nrow=1, ncol=ncol(all_results)))
    colnames(na_row) = colnames(all_results)
    na_row$p_format = variant
    all_results = rbind(all_results, na_row)
  }
}

# Re-start script if any errors occur

# Track already processed variants
processed_variants = unique(all_results$p_format)

# Main loop over variants
for (variant in missing_rows$p_format) {
  if (variant %in% processed_variants) {
    cat("Skipping already processed:", variant, "\n")
    next
  }
  
  cat("Processing:", variant, "\n")
  
  # Prepare POST body
  body = list(
    app = "transvar",
    command_id = "transvar",
    task = "panno",
    refversion = "hg38",
    refseq = "refseq",
    typed_identifiers = variant
  )
  
  # API request
  response = POST(
    url = "https://bioinformatics.mdanderson.org/transvar_dyce/",
    encode = "multipart",
    body = body
  )
  
  # Parse JSON response
  content_json = content(response, as = "text", encoding = "UTF-8")
  content_list = fromJSON(content_json)
  
  if (!is.null(content_list$results) && grepl("^https?://", content_list$results)) {
    results_url = content_list$results
    
    tryCatch({
      results_txt = read.delim(results_url, stringsAsFactors = FALSE)
      results_txt$p_format = variant
      all_results = bind_rows(all_results, results_txt)
    }, error = function(e) {
      cat("Failed to read results for", variant, "\n")
      na_row = data.frame(p_format = variant)
      all_results = bind_rows(all_results, na_row)
    })
    
  } else {
    cat("Invalid response for", variant, "\n")
    na_row = data.frame(p_format = variant)
    all_results = bind_rows(all_results, na_row)
  }
}

# Save API request results
write.csv(all_results, "transvar/transvar_API_request_results.csv", row.names = F)

## READ IN TRANSVAR RESULTS INSTEAD ##
all_results = fread("transvar/transvar_API_request_results.csv")

# Replace the term LARGE1 in p_format to LARGE in maveData 
#maveData = maveData %>% mutate(p_format = str_replace(p_format, "LARGE1", "LARGE"))

# Modify all_results dataframe
all_results = all_results %>%
  separate(
    col = coordinates.gDNA.cDNA.protein.,
    into = c("gDNA", "cDNA", "protein_coord"),
    sep = "/",
    remove = FALSE
  ) # parse coordinates column

all_results = all_results %>%
  separate(coordinates.gDNA.cDNA.protein., into = c("gDNA", "cDNA", "protein"), sep = "/", remove = FALSE) %>%
  mutate(
    chr = str_extract(gDNA, "^chr[0-9XYM]+"),
    # Remove prefix for parsing
    pos = str_remove(gDNA, "^chr[0-9XYM]+:g\\."),
    # Handle start/stop
    start = str_extract(pos, "^[0-9]+"),
    stop = str_extract(pos, "(?<=_)[0-9]+"),
    stop = ifelse(is.na(stop), start, stop),
    # Extract ref and alt from different possible formats
    ref = case_when(
      str_detect(pos, "del[ACGT]+") ~ str_match(pos, "del([ACGT]+)")[,2],
      str_detect(pos, "[0-9]+[ACGT]>[ACGT]") ~ str_match(pos, "[0-9]+([ACGT])>")[,2],
      TRUE ~ NA_character_
    ),
    alt = case_when(
      str_detect(pos, "ins[ACGT]+") ~ str_match(pos, "ins([ACGT]+)")[,2],
      str_detect(pos, ">[ACGT]$") ~ str_match(pos, ">([ACGT])$")[,2],
      TRUE ~ NA_character_
    )
  ) # parse gDNA column

#5. Join with maveData_annotated to fill in missing rows 

# Make sure parsed columns exist in all_results
cols_to_merge = c("transcript", "p_format", "chr", "start", "stop", "ref", "alt")

all_results = all_results %>%
  mutate(
    ref  = as.character(ref),
    alt  = as.character(alt),
    chr  = as.character(chr),
    start = as.integer(start),
    stop  = as.integer(stop),
    transcript = as.character(transcript),
    p_format   = as.character(p_format)
  )


all_results_unique = all_results %>%
  select(all_of(cols_to_merge)) %>%
  group_by(p_format) %>%
  slice_head(n = 1) %>%   # safer than slice(1)
  ungroup()

# Join and fill only where missing
maveData_annotated = maveData_annotated %>%
  left_join(all_results_unique %>% select(all_of(cols_to_merge)),
            by = "p_format") %>%
  mutate(
    chr = ifelse(is.na(chr.x), chr.y, chr.x),
    start = ifelse(is.na(start.x), start.y, start.x),
    stop = ifelse(is.na(stop.x), stop.y, stop.x),
    ref = ifelse(is.na(ref.x), ref.y, ref.x),
    alt = ifelse(is.na(alt.x), alt.y, alt.x)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y"))

# Count how many rows still have missing chr/start/stop/ref/alt
missing_counts = maveData_annotated %>%
  summarise(across(c(chr, start, stop, ref, alt), ~ sum(is.na(.)), .names = "missing_{col}")) #760

# Modify chromosome column
maveData_annotated = maveData_annotated %>%
  mutate(chr = gsub("^chr", "", chr, ignore.case = TRUE))

## Get ref/alt alleles on same strand
# Make sure maveData is a data.table
maveData = maveData_annotated
maveData = as.data.table(maveData)

# Temporarily add 'chr' prefix if missing
maveData[, chr_orig := chr]  # backup original
maveData[, chr := ifelse(!startsWith(chr, "chr"), paste0("chr", chr), chr)]

# Filter rows that have no missing values in key columns
maveData_valid = maveData[!is.na(chr) & !is.na(start) & !is.na(stop) &
                            !is.na(ref) & !is.na(alt)]

# Load hg38 FASTA
hg38 = FaFile("Annotations/hg38.fa")
open(hg38)

# Function to normalize alleles to the plus strand
strand_normalize = function(chr, pos, ref, alt, fa) {
  gr = GRanges(chr, IRanges(pos, pos))
  ref_genome = as.character(getSeq(fa, gr))
  
  # Identify alleles that don't match the genome
  flip_idx = which(ref != ref_genome)
  
  if(length(flip_idx) > 0) {
    ref[flip_idx] = as.character(reverseComplement(DNAStringSet(ref[flip_idx])))
    alt[flip_idx] = as.character(reverseComplement(DNAStringSet(alt[flip_idx])))
  }
  
  list(ref_plus = ref, alt_plus = alt)
}

# Apply normalization to valid rows
norm_res = strand_normalize(
  maveData_valid$chr,
  maveData_valid$start,
  maveData_valid$ref,
  maveData_valid$alt,
  hg38
)

maveData_valid[, `:=`(ref_plus = norm_res$ref_plus,
                      alt_plus = norm_res$alt_plus)]

# Merge normalized alleles back into full maveData
maveData[maveData_valid, `:=`(ref_plus = i.ref_plus,
                              alt_plus = i.alt_plus),
         on = .(chr, start, stop, ref, alt)]

# Restore original chr column and remove backup
maveData[, chr := chr_orig]
maveData[, chr_orig := NULL]

# Write out results with genomic coordinates 
write.csv(maveData, "secondary_results/data.csv", row.names = F)

## READ IN ANNOVAR RESULTS ##
maveData = read.csv("secondary_results/data.csv")
invalidInput = fread("errors/myanno.refGene.invalid_input", header=FALSE)
setnames(invalidInput, c("chr", "start", "stop", "ref_plus", "alt_plus"))

# Calculate corrected stop coordinate
invalidInput[, corrected_stop := start + nchar(ref_plus) - 1]
invalidInput[, stop := corrected_stop]

# Deduplicate fixes by key columns
fixes_unique = unique(invalidInput[, .(chr, start, stop, ref_plus, alt_plus)], by = c("chr", "start", "ref_plus", "alt_plus"))
fixes_unique = as.data.frame(fixes_unique)

# Create keys for matching
maveData$key = paste(maveData$chr, maveData$start, maveData$ref_plus, maveData$alt_plus, sep = "_")
fixes_unique$key = paste(fixes_unique$chr, fixes_unique$start, fixes_unique$ref_plus, fixes_unique$alt_plus, sep = "_")

# Named vector for corrected stops
corrected_stops = setNames(fixes_unique$stop, fixes_unique$key)

# Replace stop coordinates in maveData where keys match
maveData$stop = ifelse(
  maveData$key %in% names(corrected_stops),
  corrected_stops[maveData$key],
  maveData$stop
)

# Remove helper key column
maveData$key = NULL

# Write out to CSV to reannotate with ANNOVAR
write.csv(maveData, "secondary_results/data.csv", row.names = F)

## CLINVAR ANNOTATIONS ##
# Read in annotated MAVE variants
annovarMAVE = read.csv("secondary_results/data.csv")
geneList = unique(annovarMAVE$gene_only) #39 genes

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

# Add StarRating column based on ClinVar ReviewStatus (SCV + VCV/RCV)
varSummary_filtered = varSummary_filtered %>%
  mutate(
    StarRating = case_when(
      # 4 stars
      ReviewStatus == "practice guideline" ~ 4,
      
      # 3 stars
      ReviewStatus == "reviewed by expert panel" ~ 3,
      
      # 2 stars (aggregate only)
      ReviewStatus == "criteria provided, multiple submitters, no conflicts" ~ 2,
      
      # 1 star
      ReviewStatus %in% c(
        "criteria provided, single submitter",
        "criteria provided, conflicting classifications"
      ) ~ 1,
      
      # 0 stars
      ReviewStatus %in% c(
        "no assertion criteria provided",
        "no classification provided",
        "no classification for the individual variant"
      ) ~ 0,
      
      # Anything unexpected
      TRUE ~ NA_real_
    )
  )

# Group p_format variants by classification 
clinvarGroups = varSummary_filtered %>%
  group_by(p_format) %>%
  summarise(
    clinical_significance = paste(unique(ClinicalSignificance), collapse = "; "),
    star_rating = paste(unique(StarRating), collapse = "; "),
    clinvar_entry_count = n(),
    .groups = "drop"
  )

# Apply join to maveData
annovarMAVE_clinvar = annovarMAVE %>%
  left_join(clinvarGroups, by = "p_format")

write.csv(annovarMAVE_clinvar, "secondary_results/data.csv", row.names = F) 

## SETTING THRESHOLD CLASS FOR IN SILICO PREDICTORS ## 
maveData = fread("secondary_results/data.csv")

# Convert "." to NA first
maveData$VEST4_score = as.numeric(ifelse(maveData$VEST4_score == ".", NA, maveData$VEST4_score))
maveData$REVEL_score = as.numeric(ifelse(maveData$REVEL_score == ".", NA, maveData$REVEL_score))
maveData$MutPred_score = as.numeric(ifelse(maveData$MutPred_score == ".", NA, maveData$MutPred_score))
maveData$BayesDel_noAF_score = as.numeric(ifelse(maveData$BayesDel_noAF_score == ".", NA, maveData$BayesDel_noAF_score))

# Create new columns 
maveData$VEST4_pred = NA
maveData$REVEL_pred = NA
maveData$MutPred_pred = NA

# Apply thresholds and give labels 
maveData$VEST4_pred = ifelse(
 maveData$VEST4_score > 0.764,  
 "D", "T"
)

maveData$REVEL_pred = ifelse(
  maveData$REVEL_score > 0.644,  
  "D", "T"
)

maveData$MutPred_pred = ifelse(
  maveData$MutPred_score > 0.737,  
  "D", "T"
)

write.csv(maveData, "secondary_results/data.csv", row.names = F) 

## ANNOTATING WITH ALPHAMISSENSE ##
## Read in MAVE dataset
maveData = fread("secondary_results/data.csv")

# Read in UniprotID file
idMap = read.table("Annotations/HUMAN_9606_idmapping.dat", header = F, sep = "\t", quote = "", fill = TRUE)

# Filter idMap to keep only rows with Gene_Name
geneIdMap = idMap[idMap$V2 == "Gene_Name",]

# Create a new column 'uniprotID'
maveData$uniprotID = NA

# Iterate through unique gene names 
for (gene in unique(maveData$gene_only)) {
  # Check if the gene is present in geneIdMap
  if (gene %in% geneIdMap$V3) {
    # Extract the corresponding Uniprot IDs
    uniprotIDs = geneIdMap[geneIdMap$V3 == gene, "V1"]
    
    # Combine multiple Uniprot IDs into a single string separated by a semicolon
    uniprotID = paste(uniprotIDs, collapse = ";")
    
    # Update the uniprotID column for all rows with the same gene
    maveData$uniprotID[maveData$gene_only %in% gene] <- uniprotID
  }
}

# Only keep the first item in the semi-colon separated list of maveData$uniprotID
maveData$uniprotID = sapply(strsplit(as.character(maveData$uniprotID), ";"), function(x) x[1])

# Read in AlphaMissense scores
amScores = fread("annotations/AlphaMissense_hg38.tsv")
colnames(amScores)[colnames(amScores) == "#CHROM"] <- "CHROM"

# Create new columns in amScores
amScores = amScores %>% mutate(ref_AA = substr(protein_variant, 1, 1))
amScores = amScores %>% mutate(alt_AA = substr(protein_variant, nchar(protein_variant), nchar(protein_variant)))
amScores = amScores %>%
  mutate(aa_pos = as.integer(substring(protein_variant, 2, nchar(protein_variant) - 1)))

# Remove 'chr' prefix from CHROM
amScores$CHROM = gsub("^chr", "", amScores$CHROM)

# Match rows and add only the new columns
amScores_subset = amScores %>%
  select(CHROM, POS, REF, ALT, am_pathogenicity, am_class)

maveData_annotated = maveData %>%
  left_join(
    amScores_subset,
    by = c("chr" = "CHROM", "start" = "POS", "ref_plus" = "REF", "alt_plus" = "ALT")
  )

# Identify rows with missing AlphaMissense data
missing_idx = which(is.na(maveData_annotated$am_pathogenicity) & is.na(maveData_annotated$am_class)) #124268
missing_annovar =  which((maveData_annotated$AlphaMissense_score == ".") & (maveData_annotated$AlphaMissense_pred == ".")) #124262

amScores_unique = amScores %>%
  distinct(uniprot_id, aa_pos, ref_AA, alt_AA, am_pathogenicity, am_class, .keep_all = TRUE)

# Second pass join — match by amino acids instead
maveData_aa_match = maveData_annotated %>%
  filter(row_number() %in% missing_idx) %>%
  left_join(
    amScores_unique %>%
      select(uniprot_id, aa_pos, ref_AA, alt_AA, am_pathogenicity, am_class),
    by = c("uniprotID" = "uniprot_id", "aa_pos" = "aa_pos", "ref_aa1" = "ref_AA", "alt_aa1" = "alt_AA")
  )

# Fill in missing values from amino acid match
maveData_annotated$am_pathogenicity[missing_idx] <- maveData_aa_match$am_pathogenicity.y
maveData_annotated$am_class[missing_idx] <- maveData_aa_match$am_class.y

# Identify rows with missing AlphaMissense data
missing_idx = which(is.na(maveData_annotated$am_pathogenicity) & is.na(maveData_annotated$am_class)) #123435

write.csv(maveData_annotated, "secondary_results/data.csv", row.names = F)

## ADD AM AND VARITY PRED COLUMNS ##
## SETTING THRESHOLD CLASS FOR IN SILICO PREDICTORS ## 
maveData = fread("secondary_results/data.csv")

# Convert "." to NA first
maveData$am_pathogenicity = as.numeric(ifelse(maveData$am_pathogenicity == ".", NA, maveData$am_pathogenicity))
maveData$VARITY_R_score = as.numeric(ifelse(maveData$VARITY_R_score == ".", NA, maveData$VARITY_R_score))

# Create new columns 
maveData$AM_pred = NA
maveData$VARITY_pred = NA

# Apply thresholds and give labels 
maveData$AM_pred = ifelse(
  maveData$am_pathogenicity > 0.564, "D",
  ifelse(maveData$am_pathogenicity < 0.34, "T", "A")
)

maveData$VARITY_pred = ifelse(
  maveData$VARITY_R_score > 0.75,  
  "D", "T"
)

write.csv(maveData, "secondary_results/final.csv", row.names = F) 
