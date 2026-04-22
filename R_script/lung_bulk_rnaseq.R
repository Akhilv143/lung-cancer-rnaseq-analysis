

# Step 1: Get raw count matrix BEFORE stripping versions
count_matrix_raw <- as.data.frame(counts(dds, normalized = FALSE))

cat("Dimensions before dedup:", nrow(count_matrix_raw), "genes x",
    ncol(count_matrix_raw), "samples\n")

# Step 2: Add clean Ensembl ID column (version stripped)
count_matrix_raw$ensembl_clean <- sub("\\..*", "", rownames(count_matrix_raw))

# Step 3: Check how many duplicates exist
dup_ids <- count_matrix_raw$ensembl_clean[
  duplicated(count_matrix_raw$ensembl_clean)]
cat("Duplicate base Ensembl IDs found:", length(dup_ids), "\n")
cat("Affected IDs:\n")
print(unique(dup_ids))

# Step 4: Aggregate duplicates by SUMMING counts per base Ensembl ID

count_matrix_ensembl <- count_matrix_raw %>%
  group_by(ensembl_clean) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  column_to_rownames("ensembl_clean")

cat("\nDimensions after dedup (summed):", nrow(count_matrix_ensembl),
    "genes x", ncol(count_matrix_ensembl), "samples\n")
cat("Example Ensembl IDs (first 5):\n")
print(rownames(count_matrix_ensembl)[1:5])


# biomaRt Ensembl -> Entrez conversion

ensembl_ids_clean <- rownames(count_matrix_ensembl)

cat("\nConnecting to biomaRt...\n")
mart <- tryCatch(
  useMart("ensembl", dataset = "hsapiens_gene_ensembl",
          host = "https://useast.ensembl.org"),
  error = function(e) {
    cat("Trying mirror...\n")
    useMart("ensembl", dataset = "hsapiens_gene_ensembl",
            host = "https://uswest.ensembl.org")
  }
)

anno <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids_clean,
  mart       = mart
)

cat("Mappings returned          :", nrow(anno), "\n")
cat("Ensembl IDs with Entrez    :",
    length(unique(anno$ensembl_gene_id[!is.na(anno$entrezgene_id)])), "\n")
cat("Lost (no Entrez ID)        :",
    length(ensembl_ids_clean) -
      length(unique(anno$ensembl_gene_id[!is.na(anno$entrezgene_id)])), "\n")

# Clean 1:1 mapping
anno_clean <- anno %>%
  filter(!is.na(entrezgene_id)) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%   # 1 entrez per ensembl
  distinct(entrezgene_id,   .keep_all = TRUE)        # no duplicate entrez

cat("Clean 1:1 mappings         :", nrow(anno_clean), "\n")

# Apply mapping
count_matrix_entrez <- count_matrix_ensembl %>%
  rownames_to_column("ensembl_gene_id") %>%
  inner_join(anno_clean %>% select(ensembl_gene_id, entrezgene_id),
             by = "ensembl_gene_id") %>%
  select(-ensembl_gene_id) %>%
  mutate(entrezgene_id = as.character(entrezgene_id)) %>%
  column_to_rownames("entrezgene_id")

cat("\nFinal Entrez count matrix:\n")
cat("  Genes  :", nrow(count_matrix_entrez), "\n")
cat("  Samples:", ncol(count_matrix_entrez), "\n")
cat("  First 5 Entrez IDs:", rownames(count_matrix_entrez)[1:5], "\n")
cat("  Preview (5 x 4):\n")
print(count_matrix_entrez[1:5, 1:4])




# Metadata

metadata <- file_df %>%
  mutate(
    Sample_Name = paste(PatientLabel,
                        ifelse(Condition == "Normal", "Normal", "Tumor")),
    Condition   = ifelse(Condition == "Normal", "Control", "Tumor"),
    GSE         = "GSE283245",
    Platform    = "Illumina NovaSeq 6000 (GPL24676)"
  ) %>%
  select(GSM_ID, Sample_Name, Condition, GSE, Platform) %>%
  arrange(Sample_Name)

cat("\n--- Metadata Preview (first 8 rows) ---\n")
print(head(metadata, 8))
cat("Total metadata rows:", nrow(metadata), "\n")


# Save outputs: Count matrix as TSV + Metadata as CSV

out_dir <- "C:/Users/akhil/comet_downlaod/lung_count/GSE283245_processed"
dir.create(out_dir, showWarnings = FALSE)

# Count matrix -> TSV
count_out <- file.path(out_dir, "GSE283245_count_matrix_entrez.tsv")
fwrite(count_matrix_entrez %>% rownames_to_column("GeneID"),
       file = count_out, sep = "\t", quote = FALSE)

# Metadata -> CSV
meta_out <- file.path(out_dir, "GSE283245_metadata.csv")
write.csv(metadata, file = meta_out, row.names = FALSE)

cat("\n============================================================\n")
cat("DONE! Saved to:", out_dir, "\n")
cat("  1. GSE283245_count_matrix_entrez.tsv —",
    nrow(count_matrix_entrez), "genes (Entrez) x",
    ncol(count_matrix_entrez), "samples\n")
cat("  2. GSE283245_metadata.csv            —",
    nrow(metadata), "samples\n")
cat("============================================================\n")



# Metadata Preparation
# Datasets: GSE283245 + GSE81089 + GSE159857

library(tidyverse)
library(data.table)

setwd("C:/Users/akhil/comet_downlaod/lung_count/count_matrix_unzipped")
cat("Working directory:", getwd(), "\n")
cat("Files present:\n")
print(list.files())


# PART 1: GSE283245 metadata — already exists, just load

cat("\n--- Loading GSE283245 metadata ---\n")

meta_283 <- read.csv("GSE283245_metadata.csv", stringsAsFactors = FALSE)

# Rename columns to match lihc pipeline style
meta_283 <- meta_283 %>%
  rename(
    GSM_ID    = GSM_ID,
    Condition = Condition,    # already "Control" / "Tumor"
    GSE       = GSE,
    Platform  = Platform
  ) %>%
  select(GSM_ID, Condition, GSE, Platform)

cat("GSE283245 samples:", nrow(meta_283), "\n")
print(table(meta_283$Condition))


# PART 2: GSE81089 metadata — matched pairs only
# Normal GSMs: GSM2142642–GSM2142660 (19 normals)

cat("\n--- Building GSE81089 metadata (matched pairs) ---\n")

# Normal samples — from GEO page (GSM ID -> Sample name)
normal_gsm <- c(
  "GSM2142642", "GSM2142643", "GSM2142644", "GSM2142645", "GSM2142646",
  "GSM2142647", "GSM2142648", "GSM2142649", "GSM2142650", "GSM2142651",
  "GSM2142652", "GSM2142653", "GSM2142654", "GSM2142655", "GSM2142656",
  "GSM2142657", "GSM2142658", "GSM2142659", "GSM2142660"
)

normal_name <- c(
  "L511N","L532N","L561N","L563N","L566N",
  "L572N","L606N","L616N","L644N","L656N",
  "L661N","L682N","L723N","L724N","L736N",
  "L738N","L809N","L831N","L881N"
)

# Corresponding tumor GSMs (matched by sample number from GEO listing)
# L511T, L532T, L561T ... paired with above normals
tumor_gsm <- c(
  "GSM2142481",  # L511T
  "GSM2142485",  # L532T
  "GSM2142497",  # L561T
  "GSM2142498",  # L563T
  "GSM2142500",  # L566T
  "GSM2142504",  # L572T
  "GSM2142520",  # L606T
  "GSM2142525",  # L616T
  "GSM2142537",  # L644T
  "GSM2142541",  # L656T
  "GSM2142543",  # L661T
  "GSM2142549",  # L682T
  "GSM2142560",  # L723T
  "GSM2142561",  # L724T
  "GSM2142566",  # L736T
  "GSM2142567",  # L738T
  "GSM2142596",  # L809T
  "GSM2142609",  # L831T
  "GSM2142636"   # L881T
)

tumor_name <- gsub("N$", "T", normal_name)

# Build metadata for GSE81089
meta_81089 <- data.frame(
  GSM_ID    = c(normal_gsm,   tumor_gsm),
  Condition = c(rep("Control", 19), rep("Tumor", 19)),
  SampleName= c(normal_name,  tumor_name),
  GSE       = "GSE81089",
  Platform  = "Illumina HiSeq 2500 (GPL16791)",
  stringsAsFactors = FALSE
)

cat("GSE81089 matched samples:", nrow(meta_81089), "\n")
print(table(meta_81089$Condition))

# PART 3: GSE159857 metadata — all 58 paired samples
# Pattern: A = Adenocarcinoma, S = Squamous

cat("\n--- Building GSE159857 metadata ---\n")

# All 58 GSM IDs in order (from GEO listing)
gsm_159857 <- paste0("GSM", 4848811:4848868)

# Sample names in same order as GEO listing
sample_names_159857 <- c(
  "A10N","A10T","A12N","A14N","A14T",
  "A15N","A15T","A18N","A18T","A1N",
  "A1T", "A2N", "A2T", "A3N", "A3T",
  "A4N", "A4T", "A5N", "A5T", "A64N",
  "A65T","A6N", "A6T", "A7N", "A7T",
  "A8N", "A8T", "A9N", "A9T", "S10N",
  "S10T","S12N","S12T","S13N","S13T",
  "S14N","S14T","S15N","S15T","S16N",
  "S16T","S17N","S17T","S2N", "S2T",
  "S3N", "S3T", "S4N", "S4T", "S5N",
  "S5T", "S6N", "S6T", "S7N", "S7T",
  "S8N", "S8T", "S9T"
)

meta_159857 <- data.frame(
  GSM_ID     = gsm_159857,
  SampleName = sample_names_159857,
  Condition  = ifelse(grepl("N$", sample_names_159857), "Control", "Tumor"),
  Histology  = ifelse(grepl("^A",  sample_names_159857),
                      "Adenocarcinoma", "Squamous"),
  GSE        = "GSE159857",
  Platform   = "Illumina NextSeq 500 (GPL18573)",
  stringsAsFactors = FALSE
)

cat("GSE159857 samples:", nrow(meta_159857), "\n")
print(table(meta_159857$Condition))
print(table(meta_159857$Histology, meta_159857$Condition))


# PART 4: Combine all metadata

cat("\n--- Combining all metadata ---\n")

# Harmonise columns before binding
meta_283_clean <- meta_283 %>%
  mutate(SampleName = NA_character_,
         Histology  = NA_character_) %>%
  select(GSM_ID, SampleName, Condition, Histology, GSE, Platform)

meta_81089_clean <- meta_81089 %>%
  mutate(Histology = NA_character_) %>%
  select(GSM_ID, SampleName, Condition, Histology, GSE, Platform)

meta_159857_clean <- meta_159857 %>%
  select(GSM_ID, SampleName, Condition, Histology, GSE, Platform)

# Bind all
combined_meta <- bind_rows(meta_283_clean,
                           meta_81089_clean,
                           meta_159857_clean) %>%
  mutate(
    condition = factor(Condition, levels = c("Control", "Tumor")),  # matches lihc style
    batch     = factor(GSE,
                       levels = c("GSE283245", "GSE81089", "GSE159857"))
  )

cat("\nCombined metadata summary:\n")
cat("Total samples:", nrow(combined_meta), "\n")
print(table(combined_meta$batch, combined_meta$condition))


# PART 5: Verify GSM IDs match count matrix columns

cat("\n--- Sanity checks against count matrices ---\n")

# Load just the header of each count matrix to check column names
cols_283   <- colnames(fread("GSE283245_count_matrix_entrez.tsv",
                             nrows = 0))[-1]
cols_81089  <- colnames(fread("GSE81089_raw_counts_GRCh38.p13_NCBI.tsv",
                              nrows = 0))[-1]
cols_159857 <- colnames(fread("GSE159857_raw_counts_GRCh38.p13_NCBI.tsv",
                              nrows = 0))[-1]

# Check GSE283245
gsm_283_meta <- combined_meta$GSM_ID[combined_meta$GSE == "GSE283245"]
cat("GSE283245 — meta samples:", length(gsm_283_meta),
    "| count cols:", length(cols_283), "\n")
cat("All GSE283245 GSMs in count matrix:",
    all(gsm_283_meta %in% cols_283), "\n")

# Check GSE81089 matched only
gsm_81089_meta <- combined_meta$GSM_ID[combined_meta$GSE == "GSE81089"]
cat("GSE81089  — meta samples:", length(gsm_81089_meta),
    "| count cols:", length(cols_81089), "\n")
cat("All GSE81089 GSMs in count matrix:",
    all(gsm_81089_meta %in% cols_81089), "\n")

# Check GSE159857
gsm_159857_meta <- combined_meta$GSM_ID[combined_meta$GSE == "GSE159857"]
cat("GSE159857 — meta samples:", length(gsm_159857_meta),
    "| count cols:", length(cols_159857), "\n")
cat("All GSE159857 GSMs in count matrix:",
    all(gsm_159857_meta %in% cols_159857), "\n")


# PART 6: Save metadata files

cat("\n--- Saving metadata ---\n")

# Individual dataset metadata
write.csv(meta_81089,  "GSE81089_metadata.csv",  row.names = FALSE)
write.csv(meta_159857, "GSE159857_metadata.csv", row.names = FALSE)

# Combined metadata (used directly in DESeq2 — like lihc pipeline)
write.csv(combined_meta, "LUNG_combined_metadata.csv", row.names = FALSE)

cat("Saved:\n")
cat("  GSE81089_metadata.csv       —", nrow(meta_81089),  "rows\n")
cat("  GSE159857_metadata.csv      —", nrow(meta_159857), "rows\n")
cat("  LUNG_combined_metadata.csv  —", nrow(combined_meta), "rows\n")


# PART 7: Preview 
cat("\n--- Combined Metadata Preview (first 8 rows) ---\n")
print(head(combined_meta, 8))

cat("\n--- Batch x Condition table ---\n")
print(table(combined_meta$batch, combined_meta$condition))



# LUNG CANCER - DGE


# Libraries 
library(data.table)
library(tidyverse)
library(DESeq2)
library(sva)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(BiocParallel)
library(grid)

register(SerialParam())
options(mc.cores = 1)

# Working directory & folder structure 
setwd("C:/Users/akhil/comet_downlaod/lung_count/count_matrix_unzipped")

dir.create("results_LUNG/01_objects",           recursive = TRUE, showWarnings = FALSE)
dir.create("results_LUNG/02_tables/DEG",        recursive = TRUE, showWarnings = FALSE)
dir.create("results_LUNG/02_tables/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("results_LUNG/03_plots/QC_PCA",      recursive = TRUE, showWarnings = FALSE)
dir.create("results_LUNG/03_plots/DEG",         recursive = TRUE, showWarnings = FALSE)
dir.create("results_LUNG/03_plots/Heatmap",     recursive = TRUE, showWarnings = FALSE)
dir.create("results_LUNG/03_plots/enrichment",  recursive = TRUE, showWarnings = FALSE)


# STEP 1: Load individual metadata files

cat("=== STEP 1: Load Metadata ===\n")

#  GSE283245: already created, load directly 
meta_283 <- read.csv("GSE283245_metadata.csv", stringsAsFactors = FALSE)

# Rename to match target format (like GSE109169_metadata-2.csv)
meta_283 <- meta_283 %>%
  mutate(
    sampleid    = GSM_ID,
    geoaccession= GSM_ID,
    title       = Sample_Name,
    condition   = ifelse(Condition == "Control", "control", "diseased"),
    patientid   = str_extract(Sample_Name, "\\d+"),
    batch       = "GSE283245",
    platform    = "GPL24676"
  ) %>%
  select(sampleid, geoaccession, title, condition, patientid, batch, platform)

#  GSE81089: matched pairs only 
normal_gsm  <- c("GSM2142642","GSM2142643","GSM2142644","GSM2142645","GSM2142646",
                 "GSM2142647","GSM2142648","GSM2142649","GSM2142650","GSM2142651",
                 "GSM2142652","GSM2142653","GSM2142654","GSM2142655","GSM2142656",
                 "GSM2142657","GSM2142658","GSM2142659","GSM2142660")

normal_name <- c("L511N","L532N","L561N","L563N","L566N","L572N","L606N",
                 "L616N","L644N","L656N","L661N","L682N","L723N","L724N",
                 "L736N","L738N","L809N","L831N","L881N")

tumor_gsm   <- c("GSM2142481","GSM2142485","GSM2142497","GSM2142498","GSM2142500",
                 "GSM2142504","GSM2142520","GSM2142525","GSM2142537","GSM2142541",
                 "GSM2142543","GSM2142549","GSM2142560","GSM2142561","GSM2142566",
                 "GSM2142567","GSM2142596","GSM2142609","GSM2142636")

tumor_name  <- gsub("N$", "T", normal_name)

# Patient IDs from sample names (e.g. L511N -> 511)
patient_ids_81089 <- gsub("[NT]$", "", gsub("^L", "", normal_name))

meta_81089 <- data.frame(
  sampleid     = c(normal_gsm, tumor_gsm),
  geoaccession = c(normal_gsm, tumor_gsm),
  title        = c(normal_name, tumor_name),
  condition    = c(rep("control", 19), rep("diseased", 19)),
  patientid    = c(patient_ids_81089, patient_ids_81089),
  batch        = "GSE81089",
  platform     = "GPL16791",
  stringsAsFactors = FALSE
)

# GSE159857: all 58 paired samples 
gsm_159857 <- paste0("GSM", 4848811:4848868)

sample_names_159857 <- c(
  "A10N","A10T","A12N","A14N","A14T","A15N","A15T","A18N","A18T","A1N",
  "A1T", "A2N", "A2T", "A3N", "A3T", "A4N", "A4T", "A5N", "A5T", "A64N",
  "A65T","A6N", "A6T", "A7N", "A7T", "A8N", "A8T", "A9N", "A9T", "S10N",
  "S10T","S12N","S12T","S13N","S13T","S14N","S14T","S15N","S15T","S16N",
  "S16T","S17N","S17T","S2N", "S2T", "S3N", "S3T", "S4N", "S4T", "S5N",
  "S5T", "S6N", "S6T", "S7N", "S7T", "S8N", "S8T", "S9T"
)

# Patient IDs: extract number from sample name (A10N -> 10, S10N -> 10)
patient_ids_159857 <- gsub("[ANT]", "", sample_names_159857)
patient_ids_159857 <- paste0(
  ifelse(grepl("^A", sample_names_159857), "A", "S"),
  patient_ids_159857
)

meta_159857 <- data.frame(
  sampleid     = gsm_159857,
  geoaccession = gsm_159857,
  title        = sample_names_159857,
  condition    = ifelse(grepl("N$", sample_names_159857), "control", "diseased"),
  patientid    = patient_ids_159857,
  batch        = "GSE159857",
  platform     = "GPL18573",
  stringsAsFactors = FALSE
)

#  Combine all metadata
combined_meta <- bind_rows(meta_283, meta_81089, meta_159857)

cat("Combined metadata — total samples:", nrow(combined_meta), "\n")
print(table(combined_meta$batch, combined_meta$condition))


# STEP 2: Load count matrices

cat("\n=== STEP 2: Load Count Matrices ===\n")

load_counts <- function(file) {
  df           <- fread(file, header = TRUE) %>% as.data.frame()
  rownames(df) <- df[[1]]
  df[[1]]      <- NULL
  as.matrix(df)
}

mat_283   <- load_counts("GSE283245_count_matrix_entrez.tsv")
mat_81089 <- load_counts("GSE81089_raw_counts_GRCh38.p13_NCBI.tsv")
mat_159857<- load_counts("GSE159857_raw_counts_GRCh38.p13_NCBI.tsv")

cat("GSE283245  matrix:", nrow(mat_283),    "genes x", ncol(mat_283),    "samples\n")
cat("GSE81089   matrix:", nrow(mat_81089),  "genes x", ncol(mat_81089),  "samples\n")
cat("GSE159857  matrix:", nrow(mat_159857), "genes x", ncol(mat_159857), "samples\n")

#  Subset GSE81089 to matched samples only 
keep_81089 <- c(normal_gsm, tumor_gsm)
mat_81089  <- mat_81089[, keep_81089]
cat("GSE81089 after subsetting matched pairs:",
    ncol(mat_81089), "samples\n")

# Find common genes across all 3 datasets 
common_genes <- Reduce(intersect, list(
  rownames(mat_283),
  rownames(mat_81089),
  rownames(mat_159857)
))
cat("Common Entrez genes across all 3 datasets:", length(common_genes), "\n")

mat_283    <- mat_283[common_genes, ]
mat_81089  <- mat_81089[common_genes, ]
mat_159857 <- mat_159857[common_genes, ]

# Verify gene order
stopifnot(
  all(rownames(mat_283)   == rownames(mat_81089)),
  all(rownames(mat_283)   == rownames(mat_159857))
)

#  Merge 
raw_counts <- cbind(mat_283, mat_81089, mat_159857)
raw_counts <- raw_counts[rowSums(raw_counts) > 0, ]
cat("Merged matrix:", nrow(raw_counts), "genes x", ncol(raw_counts), "samples\n")


# STEP 3: Align metadata to merged count matrix

cat("\n=== STEP 3: Align Metadata ===\n")

rownames(combined_meta) <- combined_meta$sampleid

# Keep only samples present in count matrix
combined_meta <- combined_meta[colnames(raw_counts), , drop = FALSE]

# Add DESeq2-style factor columns  (matching lihc pipeline exactly)
combined_meta$condition <- factor(combined_meta$condition,
                                  levels = c("control", "diseased"))
combined_meta$batch     <- factor(combined_meta$batch,
                                  levels = c("GSE283245","GSE81089","GSE159857"))

stopifnot(all(colnames(raw_counts) %in% rownames(combined_meta)))
stopifnot(all(rownames(combined_meta) %in% colnames(raw_counts)))
stopifnot(!any(is.na(combined_meta$condition)))

cat("Metadata aligned — samples:", nrow(combined_meta), "\n")
print(table(combined_meta$batch, combined_meta$condition))


# STEP 4: Map Entrez IDs -> Gene Symbols for merged matrix

cat("\n=== STEP 4: Map Gene Symbols ===\n")

entrez_ids <- rownames(raw_counts)

symbol_map <- mapIds(
  org.Hs.eg.db,
  keys      = entrez_ids,
  keytype   = "ENTREZID",
  column    = "SYMBOL",
  multiVals = "first"
)

gene_annotation <- data.frame(
  EntrezID = entrez_ids,
  Symbol   = unname(symbol_map),
  stringsAsFactors = FALSE
)

cat("Entrez IDs with gene symbols:",
    sum(!is.na(gene_annotation$Symbol)), "/", nrow(gene_annotation), "\n")


# STEP 5: Save merged count matrix + metadata

cat("\n=== STEP 5: Save Merged Files ===\n")

# Count matrix: GeneID (Entrez) + Symbol + sample columns
merged_counts_df <- as.data.frame(raw_counts) %>%
  rownames_to_column("EntrezID") %>%
  left_join(gene_annotation, by = "EntrezID") %>%
  select(EntrezID, Symbol, everything())   # EntrezID | Symbol | GSM... | GSM...

fwrite(merged_counts_df,
       "LUNG_merged_raw_counts.tsv",
       sep = "\t", quote = FALSE)

# Metadata: 
write.csv(combined_meta,
          "LUNG_combined_metadata.csv",
          row.names = FALSE)

cat("Saved:\n")
cat("  LUNG_merged_raw_counts.tsv  —",
    nrow(raw_counts), "genes x", ncol(raw_counts), "samples\n")
cat("  LUNG_combined_metadata.csv  —", nrow(combined_meta), "rows\n")

cat("\nMerged count matrix preview (first 5 rows, 5 cols):\n")
print(merged_counts_df[1:5, 1:5])

cat("\nMetadata preview (first 8 rows):\n")
print(head(combined_meta, 8))


# STEP 6: ComBat-seq batch correction

cat("\n=== STEP 6: ComBat-seq Batch Correction ===\n")

batch_vec <- as.integer(combined_meta$batch)
group_vec <- as.integer(combined_meta$condition)

cat("Batch vector table:\n")
print(table(combined_meta$batch))

corrected_counts <- ComBat_seq(
  counts = raw_counts,
  batch  = batch_vec,
  group  = group_vec
)

cat("ComBat-seq complete.\n")
cat("Corrected matrix:", nrow(corrected_counts), "genes x",
    ncol(corrected_counts), "samples\n")


# STEP 7: DESeq2

cat("\n=== STEP 7: DESeq2 ===\n")

dds <- DESeqDataSetFromMatrix(
  countData = corrected_counts,
  colData   = combined_meta,
  design    = ~ batch + condition
)

# Pre-filter: keep genes with >= 10 counts in >= 10 samples
keep <- rowSums(counts(dds) >= 10) >= 10
dds  <- dds[keep, ]
cat("Genes after pre-filtering:", nrow(dds), "\n")

dds <- DESeq(dds, parallel = FALSE)

summary(results(dds,
                contrast = c("condition", "diseased", "control"),
                alpha    = 0.05))



# STEP 8: Extract DESeq2 Results

cat("\n=== STEP 8: Extract DESeq2 Results ===\n")

res <- results(dds,
               contrast = c("condition", "diseased", "control"),
               alpha    = 0.05)

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  left_join(
    gene_annotation %>% select(EntrezID, Symbol),
    by = c("gene" = "EntrezID")
  ) %>%
  rename(symbol = Symbol) %>%
  mutate(
    regulation = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >  1 ~ "Up",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(padj)

message("Up: ",   sum(res_df$regulation == "Up",   na.rm = TRUE),
        " | Down: ", sum(res_df$regulation == "Down", na.rm = TRUE))

write.csv(res_df,
          "results_LUNG/02_tables/DEG/DESeq2_all_genes_with_symbols.csv",
          row.names = FALSE)


# PCA: BEFORE AFTER ComBat-seq

cat("\n=== PCA Before vs After ComBat-seq ===\n")

library(ggplot2)
library(patchwork)

# ---- Helper: PCA from raw integer count matrix 
pca_from_matrix <- function(mat, meta_df, title_label) {
  
  # VST-like: log1p normalization for visualization only
  log_mat <- log1p(mat)
  
  # PCA on transposed matrix (samples as rows)
  pca_res     <- prcomp(t(log_mat), scale. = TRUE, center = TRUE)
  percentVar  <- round(100 * summary(pca_res)$importance[2, 1:2], 1)
  
  pca_df <- data.frame(
    PC1       = pca_res$x[, 1],
    PC2       = pca_res$x[, 2],
    condition = meta_df[rownames(pca_res$x), "condition"],
    batch     = meta_df[rownames(pca_res$x), "batch"],
    row.names = rownames(pca_res$x)
  )
  
  # Condition plot
  p_cond <- ggplot(pca_df, aes(PC1, PC2,
                               colour = condition,
                               fill   = condition)) +
    geom_point(size = 2.5, alpha = 0.85) +
    stat_ellipse(geom = "polygon", level = 0.95,
                 alpha = 0.10, colour = NA) +
    stat_ellipse(level = 0.95, linewidth = 0.9) +
    scale_colour_manual(
      name   = "Condition",
      values = c("control"  = "#2166AC",
                 "diseased" = "#D6604D")
    ) +
    scale_fill_manual(
      name   = "Condition",
      values = c("control"  = "#2166AC",
                 "diseased" = "#D6604D")
    ) +
    labs(
      title = paste0(title_label, " — by Condition"),
      x     = paste0("PC1: ", percentVar[1], "% variance"),
      y     = paste0("PC2: ", percentVar[2], "% variance")
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title            = element_text(face = "bold", hjust = 0.5,
                                           size = 11),
      axis.title            = element_text(face = "bold"),
      axis.text             = element_text(face = "bold"),
      legend.title          = element_text(face = "bold"),
      legend.text           = element_text(face = "bold"),
      legend.position       = "right",
      legend.background     = element_rect(colour = "black", linewidth = 0.4,
                                           fill = "white"),
      panel.grid.minor      = element_blank()
    )
  
  # Batch plot
  p_batch <- ggplot(pca_df, aes(PC1, PC2, colour = batch)) +
    geom_point(size = 2.5, alpha = 0.85) +
    stat_ellipse(level = 0.95, linewidth = 0.9) +
    labs(
      title  = paste0(title_label, " — by Batch"),
      x      = paste0("PC1: ", percentVar[1], "% variance"),
      y      = paste0("PC2: ", percentVar[2], "% variance"),
      colour = "GEO Study"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title            = element_text(face = "bold", hjust = 0.5,
                                           size = 11),
      axis.title            = element_text(face = "bold"),
      axis.text             = element_text(face = "bold"),
      legend.title          = element_text(face = "bold"),
      legend.text           = element_text(face = "bold"),
      legend.position       = "right",
      legend.background     = element_rect(colour = "black", linewidth = 0.4,
                                           fill = "white"),
      panel.grid.minor      = element_blank()
    )
  
  list(condition = p_cond, batch = p_batch)
}

# ---- Generate PCA plots for BEFORE and AFTER 
cat("Computing PCA before ComBat-seq...\n")
pca_before <- pca_from_matrix(raw_counts,        combined_meta, "Before ComBat-seq")

cat("Computing PCA after ComBat-seq...\n")
pca_after  <- pca_from_matrix(corrected_counts,  combined_meta, "After ComBat-seq")

# ---- Combine: 2x2 grid (Condition row | Batch row) 
# Row 1: Before Condition | After Condition
# Row 2: Before Batch     | After Batch

combined_condition <- (pca_before$condition | pca_after$condition) +
  plot_annotation(
    title = "PCA by Condition — Before vs After ComBat-seq",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
    )
  )

combined_batch <- (pca_before$batch | pca_after$batch) +
  plot_annotation(
    title = "PCA by Batch — Before vs After ComBat-seq",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
    )
  )

combined_all <- (pca_before$condition | pca_after$condition) /
  (pca_before$batch     | pca_after$batch) +
  plot_annotation(
    title = "PCA Before vs After ComBat-seq Batch Correction",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
    )
  )

#  Save PNG 600 DPI 
dir.create("results_LUNG/03_plots/QC_PCA", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "results_LUNG/03_plots/QC_PCA/PCA_Condition_before_vs_after_ComBat.png",
  combined_condition,
  width = 14, height = 6, dpi = 600
)

ggsave(
  "results_LUNG/03_plots/QC_PCA/PCA_Batch_before_vs_after_ComBat.png",
  combined_batch,
  width = 14, height = 6, dpi = 600
)

ggsave(
  "results_LUNG/03_plots/QC_PCA/PCA_2x2_before_vs_after_ComBat.png",
  combined_all,
  width = 14, height = 12, dpi = 600
)


#  Non-coding filter, VST, Save Objects,
#               MA plot, Volcano plot


# step11_Non-coding filter 
noncoding_patterns <- c(
  "^LOC[0-9]", "^MIR[0-9]", "^LINC[0-9]", "^MT-",
  "^SNORD[0-9]", "^SNORA[0-9]", "^RNU[0-9]", "^RN7SL",
  "^SNHG[0-9]", "^SCARNA[0-9]", "^MALAT", "^NEAT",
  "^H19$", "^XIST$", "^MEG[0-9]", "^PEG[0-9]"
)

is_noncoding <- function(symbols) {
  pattern       <- paste(noncoding_patterns, collapse = "|")
  result        <- rep(FALSE, length(symbols))
  valid         <- !is.na(symbols)
  result[valid] <- grepl(pattern, symbols[valid], ignore.case = FALSE)
  result
}

res_df$is_noncoding <- is_noncoding(res_df$symbol)

sig_deg <- res_df %>%
  filter(regulation != "NS") %>%
  arrange(padj)

sig_deg_protein <- sig_deg %>%
  filter(!is_noncoding, !is.na(symbol))

message("Protein-coding DEGs: ", nrow(sig_deg_protein),
        " | Up: ",   sum(sig_deg_protein$regulation == "Up"),
        " | Down: ", sum(sig_deg_protein$regulation == "Down"))

write.csv(sig_deg,
          "results_LUNG/02_tables/DEG/DESeq2_significant_DEGs_all.csv",
          row.names = FALSE)
write.csv(sig_deg_protein,
          "results_LUNG/02_tables/DEG/DESeq2_significant_DEGs_proteinCoding.csv",
          row.names = FALSE)

# step12_VST normalization 
vst_mat  <- vst(dds, blind = FALSE)
vst_data <- assay(vst_mat)

write.csv(vst_data,
          "results_LUNG/02_tables/DEG/VST_normalized_expression.csv",
          row.names = TRUE)

cat("VST matrix dimensions:", nrow(vst_data), "genes x",
    ncol(vst_data), "samples\n")

# step13_Save R objects 
saveRDS(dds,             "results_LUNG/01_objects/LUNG_dds.rds")
saveRDS(res_df,          "results_LUNG/01_objects/LUNG_res_df.rds")
saveRDS(sig_deg_protein, "results_LUNG/01_objects/LUNG_sig_deg_protein.rds")
saveRDS(vst_mat,         "results_LUNG/01_objects/LUNG_vst_mat.rds")
saveRDS(vst_data,        "results_LUNG/01_objects/LUNG_vst_data.rds")

cat("R objects saved to: results_LUNG/01_objects/\n")

# step14_Plot helper: PNG + TIFF 600 DPI, no compression 
save_plot <- function(plot_obj, base_path, w, h) {
  ggsave(paste0(base_path, ".png"),
         plot_obj, width = w, height = h, dpi = 600)
  ggsave(paste0(base_path, ".tiff"),
         plot_obj, width = w, height = h, dpi = 600,
         compression = "none")
}

pal <- c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey75")

# step15_MA plot 
ma_plot <- ggplot(res_df,
                  aes(x      = log10(baseMean),
                      y      = log2FoldChange,
                      colour = regulation)) +
  geom_point(alpha = 0.4, size = 1.2, stroke = 0) +
  geom_hline(yintercept = c(-1, 1),
             linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0,
             linetype = "solid",  colour = "black", linewidth = 0.3) +
  scale_colour_manual(name = "Regulation", values = pal) +
  labs(
    x = "Mean expression (log\u2081\u2080)",
    y = "Log\u2082 fold change  (Diseased vs Control)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.title           = element_text(face = "bold"),
    axis.text            = element_text(face = "bold"),
    legend.title         = element_text(face = "bold"),
    legend.text          = element_text(face = "bold"),
    legend.position      = "right",
    legend.justification = "top",
    legend.background    = element_rect(colour = "black", linewidth = 0.4),
    panel.grid.minor     = element_blank()
  )

save_plot(ma_plot, "results_LUNG/03_plots/DEG/MA_plot", 8, 6)
message("MA plot saved.")


# step16_Volcano plot 
top10_up <- res_df %>%
  filter(regulation == "Up", !is.na(symbol), !is_noncoding) %>%
  arrange(padj) %>%
  slice_head(n = 10)

top10_dn <- res_df %>%
  filter(regulation == "Down", !is.na(symbol), !is_noncoding) %>%
  arrange(padj) %>%
  slice_head(n = 10)

top20_labeled <- bind_rows(top10_up, top10_dn)

message("Volcano labels — Up:   ", paste(top10_up$symbol, collapse = ", "))
message("Volcano labels — Down: ", paste(top10_dn$symbol, collapse = ", "))

volcano_plot <- ggplot(res_df,
                       aes(x      = log2FoldChange,
                           y      = -log10(padj),
                           colour = regulation)) +
  geom_point(alpha = 0.5, size = 1.5, stroke = 0) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "black", linewidth = 0.5) +
  ggrepel::geom_text_repel(
    data          = top20_labeled,
    aes(label     = symbol),
    size          = 3.5,
    fontface      = "bold.italic",
    max.overlaps  = 20,
    box.padding   = 0.5,
    point.padding = 0.3,
    show.legend   = FALSE
  ) +
  scale_colour_manual(
    name   = "Regulation",
    values = c("Up"   = "#D62728",
               "Down" = "#1F77B4",
               "NS"   = "grey70"),
    labels = c("Up"   = "Up-regulated",
               "Down" = "Down-regulated",
               "NS"   = "NS")
  ) +
  labs(
    x = expression(bold(Log[2]~Fold~Change~(Diseased~vs~Control))),
    y = expression(bold(-log[10](italic(p)[adj])))
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text             = element_text(face = "bold"),
    axis.title            = element_text(face = "bold"),
    legend.title          = element_text(face = "bold"),
    legend.text           = element_text(face = "bold"),
    legend.position       = "right",
    legend.justification  = "top",
    legend.key.size       = unit(0.5, "cm"),
    legend.background     = element_rect(colour = "black",
                                         linewidth = 0.5,
                                         fill = "white"),
    legend.box.background = element_rect(colour = "black",
                                         linewidth = 0.5),
    panel.grid.minor      = element_blank(),
    plot.title            = element_blank()
  )

save_plot(volcano_plot, "results_LUNG/03_plots/DEG/Volcano_plot", 8, 7)
message("Volcano plot saved.")



# step17_PCA plots

pca_df     <- plotPCA(vst_mat,
                      intgroup   = c("condition", "batch"),
                      returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

# PCA by Condition
pca_cond <- ggplot(pca_df,
                   aes(PC1, PC2,
                       colour = condition,
                       fill   = condition)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(geom = "polygon", level = 0.95,
               alpha = 0.12, colour = NA) +
  stat_ellipse(level = 0.95, linewidth = 1.1) +
  scale_colour_manual(
    name   = "Condition",
    values = c("control"  = "#2166AC",
               "diseased" = "#D6604D")
  ) +
  scale_fill_manual(
    name   = "Condition",
    values = c("control"  = "#2166AC",
               "diseased" = "#D6604D")
  ) +
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.title            = element_text(face = "bold"),
    axis.text             = element_text(face = "bold"),
    legend.title          = element_text(face = "bold"),
    legend.text           = element_text(face = "bold"),
    legend.position       = "right",
    legend.background     = element_rect(colour = "black", linewidth = 0.5,
                                         fill = "white"),
    legend.box.background = element_rect(colour = "black", linewidth = 0.5),
    panel.grid.minor      = element_blank()
  )

save_plot(pca_cond, "results_LUNG/03_plots/QC_PCA/PCA_by_Condition", 10, 7)

# PCA by Batch
pca_batch <- ggplot(pca_df,
                    aes(PC1, PC2, colour = batch)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    x      = paste0("PC1: ", percentVar[1], "% variance"),
    y      = paste0("PC2: ", percentVar[2], "% variance"),
    colour = "GEO Study"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.title            = element_text(face = "bold"),
    axis.text             = element_text(face = "bold"),
    legend.title          = element_text(face = "bold"),
    legend.text           = element_text(face = "bold"),
    legend.position       = "right",
    legend.background     = element_rect(colour = "black", linewidth = 0.5,
                                         fill = "white"),
    legend.box.background = element_rect(colour = "black", linewidth = 0.5),
    panel.grid.minor      = element_blank()
  )

save_plot(pca_batch, "results_LUNG/03_plots/QC_PCA/PCA_by_Batch", 10, 7)

message("PCA plots saved.")


#  Heatmap: Top 25 Up + Top 25 Down DEGs

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

#  Step A: Confirm res_df has all required columns 
message("res_df columns: ", paste(colnames(res_df), collapse = ", "))
message("Up: ",   sum(res_df$regulation == "Up",   na.rm = TRUE),
        " | Down: ", sum(res_df$regulation == "Down", na.rm = TRUE))

#  Step B: is_noncoding flag (if not already in res_df) 
if (!"is_noncoding" %in% colnames(res_df)) {
  noncoding_patterns <- c(
    "^LOC[0-9]", "^MIR[0-9]", "^LINC[0-9]", "^MT-",
    "^SNORD[0-9]", "^SNORA[0-9]", "^RNU[0-9]", "^RN7SL",
    "^SNHG[0-9]", "^SCARNA[0-9]", "^MALAT", "^NEAT",
    "^H19$", "^XIST$", "^MEG[0-9]", "^PEG[0-9]"
  )
  is_noncoding <- function(symbols) {
    pattern       <- paste(noncoding_patterns, collapse = "|")
    result        <- rep(FALSE, length(symbols))
    valid         <- !is.na(symbols)
    result[valid] <- grepl(pattern, symbols[valid], ignore.case = FALSE)
    result
  }
  res_df$is_noncoding <- is_noncoding(res_df$symbol)
  message("is_noncoding flag added.")
} else {
  message("is_noncoding flag already present.")
}

# Step C: Select top 25 Up + top 25 Down 
top25_up <- res_df %>%
  filter(!is.na(padj), !is.na(symbol),
         !is_noncoding, regulation == "Up") %>%
  arrange(padj) %>%
  slice_head(n = 25)

top25_dn <- res_df %>%
  filter(!is.na(padj), !is.na(symbol),
         !is_noncoding, regulation == "Down") %>%
  arrange(padj) %>%
  slice_head(n = 25)

top_degs <- bind_rows(top25_dn, top25_up)   
message("Heatmap genes — Up: ",   nrow(top25_up),
        " | Down: ", nrow(top25_dn),
        " | Total: ", nrow(top_degs))

# ---- Step D: Build heatmap matrix 
in_vst   <- top_degs$gene %in% rownames(vst_data)
top_filt <- top_degs[in_vst, ]

stopifnot("No genes found in vst_data" = nrow(top_filt) > 0)
message("Genes found in vst_data: ", nrow(top_filt), " / ", nrow(top_degs))


meta <- combined_meta


col_order <- c(
  rownames(meta)[meta$condition == "control"],
  rownames(meta)[meta$condition == "diseased"]
)

# Subset vst_data
hmat           <- vst_data[top_filt$gene, col_order, drop = FALSE]
rownames(hmat) <- top_filt$symbol

# Z-score per gene (row-wise)
hmat_z <- t(scale(t(hmat)))

# Remove rows with NA or Inf
valid  <- apply(hmat_z, 1,
                function(x) !all(is.na(x)) & !any(is.infinite(x)))
hmat_z <- hmat_z[valid, , drop = FALSE]
message("Genes retained after Z-score filter: ", nrow(hmat_z))

# Sync top_filt to valid rows only
top_filt_valid <- top_filt[top_filt$symbol %in% rownames(hmat_z), ]

# Row order: Down on top, Up on bottom
reg_vec <- top_filt_valid$regulation[
  match(rownames(hmat_z), top_filt_valid$symbol)]
row_ord <- c(which(reg_vec == "Down"), which(reg_vec == "Up"))
hmat_z  <- hmat_z[row_ord, ]
reg_vec <- reg_vec[row_ord]

# Column split by condition
col_split <- factor(
  meta[col_order, "condition"],
  levels = c("control", "diseased")
)

# Step E: Colours 

# Condition annotation colors
condition_colors <- c(
  "control"  = "#4575B4",   # blue
  "diseased" = "#FFD700"    # yellow
)

# Batch / Dataset colors
n_batch      <- length(levels(meta$batch))
batch_colors <- setNames(
  brewer.pal(max(3, n_batch), "Set2")[1:n_batch],
  levels(meta$batch)
)

# Regulation colors — NEW
reg_colors <- c(
  "Up"   = "#8BC34A",   # lime green
  "Down" = "#CE93D8"    # lavender
)

# Heatmap Z-score color scale
col_fun <- colorRamp2(c(-2, 0, 2), c("green3", "black", "red"))

# ---- Step F: Annotations 
ann_df <- data.frame(
  Condition = meta[col_order, "condition"],
  Dataset   = as.character(meta[col_order, "batch"]),
  row.names = col_order
)

# Top annotation: Condition + Dataset
ha_top <- HeatmapAnnotation(
  Condition = ann_df$Condition,
  Dataset   = ann_df$Dataset,
  col = list(
    Condition = condition_colors,
    Dataset   = batch_colors
  ),
  annotation_name_gp   = gpar(fontface = "bold", fontsize = 10),
  show_annotation_name = TRUE,
  show_legend          = FALSE
)

# Left annotation: Regulation — NEW COLORS
ha_left <- rowAnnotation(
  Regulation = reg_vec,
  col = list(
    Regulation = c(
      "Up"   = "#8BC34A",   # lime green
      "Down" = "#CE93D8"    # lavender
    )
  ),
  width                = unit(0.4, "cm"),
  annotation_name_gp   = gpar(fontface = "bold", fontsize = 10),
  annotation_name_side = "top",
  show_annotation_name = TRUE,
  show_legend          = FALSE
)

# ---- Step G: Heatmap object 
ht <- Heatmap(
  hmat_z,
  name                  = "Z-score",
  col                   = col_fun,
  top_annotation        = ha_top,
  left_annotation       = ha_left,
  cluster_rows          = FALSE,
  cluster_columns       = TRUE,
  cluster_column_slices = FALSE,
  column_split          = col_split,
  column_gap            = unit(0, "mm"),
  rect_gp               = gpar(col = NA),
  show_row_names        = TRUE,
  show_column_names     = FALSE,
  row_names_side        = "right",
  row_names_gp          = gpar(fontface = "bold.italic", fontsize = 9),
  show_heatmap_legend   = FALSE,
  row_split             = factor(reg_vec, levels = c("Down", "Up")),
  row_gap               = unit(0, "mm"),
  row_title             = NULL,
  column_title          = NULL,
  border                = TRUE
)

#  Step H: Manual legends 

# Condition legend
lgd_condition <- Legend(
  title     = "Condition",
  at        = c("control", "diseased"),
  labels    = c("Control (Normal)", "Diseased (Tumor)"),
  legend_gp = gpar(fill = c("#4575B4", "#FFD700")),
  title_gp  = gpar(fontface = "bold", fontsize = 10),
  labels_gp = gpar(fontface = "bold", fontsize = 9)
)

# Regulation legend — NEW COLORS
lgd_regulation <- Legend(
  title     = "Regulation",
  at        = c("Down", "Up"),
  labels    = c("Down-regulated", "Up-regulated"),
  legend_gp = gpar(fill = c("#CE93D8", "#8BC34A")),  # Down=lavender, Up=green
  title_gp  = gpar(fontface = "bold", fontsize = 10),
  labels_gp = gpar(fontface = "bold", fontsize = 9)
)

# Dataset legend
lgd_dataset <- Legend(
  title     = "Dataset",
  at        = levels(meta$batch),
  legend_gp = gpar(fill = batch_colors),
  title_gp  = gpar(fontface = "bold", fontsize = 10),
  labels_gp = gpar(fontface = "bold", fontsize = 9)
)

# Z-score legend
lgd_zscore <- Legend(
  title         = "Z-score",
  col_fun       = col_fun,
  at            = c(-2, 0, 2),
  labels        = c("-2", "0", "2"),
  legend_height = unit(3.5, "cm"),
  title_gp      = gpar(fontface = "bold", fontsize = 10),
  labels_gp     = gpar(fontface = "bold", fontsize = 9),
  direction     = "vertical"
)

# Pack all legends vertically
all_legends <- packLegend(
  lgd_condition,
  lgd_regulation,
  lgd_dataset,
  lgd_zscore,
  direction = "vertical",
  gap       = unit(4, "mm")
)

#  Step I: Draw function 
draw_ht <- function() {
  draw(
    ht,
    annotation_legend_list  = all_legends,
    heatmap_legend_side     = "right",
    annotation_legend_side  = "right",
    align_annotation_legend = "heatmap_top",
    merge_legend            = FALSE,
    padding                 = unit(c(2, 2, 2, 2), "mm")
  )
}

# Step J: Save PNG + TIFF 
dir.create("results_LUNG/03_plots/Heatmap",
           recursive = TRUE, showWarnings = FALSE)

# PNG
png("results_LUNG/03_plots/Heatmap/Heatmap_top25up_top25dn.png",
    width = 14, height = 11, units = "in", res = 600)
draw_ht()
dev.off()

# TIFF
tiff("results_LUNG/03_plots/Heatmap/Heatmap_top25up_top25dn.tiff",
     width = 14, height = 11, units = "in", res = 600,
     compression = "none")
draw_ht()
dev.off()

message("Heatmap saved to results_LUNG/03_plots/Heatmap/")
message("  Heatmap_top25up_top25dn.png  (14x11 in, 600 DPI)")
message("  Heatmap_top25up_top25dn.tiff (14x11 in, 600 DPI)")

#go and kegg

library(tidyverse)
library(ggplot2)
library(scales)

setwd("C:/Users/akhil/comet_downlaod/lung_count/count_matrix_unzipped")


# STEP 1: Create organised folder structure

cat("Creating folder structure...\n")

folders <- c(
  "results_LUNG/03_plots/enrichment/png",
  "results_LUNG/03_plots/enrichment/svg",
  "results_LUNG/03_plots/enrichment/tiff"
)

for (d in folders) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  cat("  Created:", d, "\n")
}

cat("Folder structure ready.\n\n")


# STEP 2: Verify enrichment objects in environment

cat("=== CHECKING OBJECTS IN ENVIRONMENT ===\n")

obj_check <- data.frame(
  Object  = c("go_bp_simp", "go_cc_simp", "go_mf_simp", "kegg"),
  Exists  = c(exists("go_bp_simp"), exists("go_cc_simp"),
              exists("go_mf_simp"), exists("kegg")),
  Rows    = c(
    if (exists("go_bp_simp") && !is.null(go_bp_simp)) nrow(go_bp_simp) else 0,
    if (exists("go_cc_simp") && !is.null(go_cc_simp)) nrow(go_cc_simp) else 0,
    if (exists("go_mf_simp") && !is.null(go_mf_simp)) nrow(go_mf_simp) else 0,
    if (exists("kegg")       && !is.null(kegg))        nrow(kegg)       else 0
  )
)

print(obj_check)

# Stop if objects missing
if (any(!obj_check$Exists)) {
  stop("Some enrichment objects are missing from environment!\n",
       "Please re-run the enrichment script first.")
}
cat("All objects confirmed in environment.\n\n")


save_3fmt <- function(p, fname_base, width = 11, height = 10) {
  cat("  Saving", fname_base, "...\n")
  
  # PNG
  ggsave(
    file.path("results_LUNG/03_plots/enrichment/png",
              paste0(fname_base, ".png")),
    plot   = p,
    width  = width,
    height = height,
    dpi    = 600,
    units  = "in"
  )
  
  # SVG
  ggsave(
    file.path("results_LUNG/03_plots/enrichment/svg",
              paste0(fname_base, ".svg")),
    plot   = p,
    width  = width,
    height = height
  )
  
  # TIFF (uncompressed)
  ggsave(
    file.path("results_LUNG/03_plots/enrichment/tiff",
              paste0(fname_base, ".tiff")),
    plot        = p,
    width       = width,
    height      = height,
    dpi         = 600,
    compression = "none"
  )
  
  # Confirm with file size
  png_size  <- round(file.info(file.path("results_LUNG/03_plots/enrichment/png",
                                         paste0(fname_base,".png")))$size / 1024, 1)
  tiff_size <- round(file.info(file.path("results_LUNG/03_plots/enrichment/tiff",
                                         paste0(fname_base,".tiff")))$size / 1024, 1)
  
  cat("    PNG :", png_size,  "KB |",
      "TIFF:", tiff_size, "KB |",
      "saved at", format(Sys.time(), "%H:%M:%S"), "\n")
}


# STEP 4: Dot plot function

make_dotplot <- function(enrich_obj, title_str, top_n = 20) {
  
  df <- as.data.frame(enrich_obj) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    arrange(Count) %>%
    mutate(
      Description = str_wrap(Description, width = 45),
      Description = factor(Description, levels = unique(Description))
    )
  
  if (nrow(df) == 0) {
    message("No terms to plot for: ", title_str)
    return(NULL)
  }
  
  p_min <- min(df$p.adjust, na.rm = TRUE)
  p_max <- max(df$p.adjust, na.rm = TRUE)
  if (p_min == p_max) p_max <- p_min * 1.01
  
  ggplot(df, aes(x = Count, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.95) +
    
    scale_color_gradientn(
      name   = "p.adjust",
      colors = c("#2D004B","#6A0573","#C2185B","#F06292","#FFCCBC"),
      values = scales::rescale(c(
        p_min,
        p_min + (p_max - p_min) * 0.25,
        p_min + (p_max - p_min) * 0.50,
        p_min + (p_max - p_min) * 0.75,
        p_max
      )),
      limits = c(p_min, p_max),
      guide  = guide_colorbar(
        barheight      = unit(5,    "cm"),
        barwidth       = unit(0.55, "cm"),
        title.position = "top",
        title.hjust    = 0.5,
        frame.colour   = "black",
        ticks.colour   = "black"
      )
    ) +
    
    scale_size_continuous(
      name   = "GeneCount",
      range  = c(3, 14),
      breaks = c(5, 10, 15, 20, 25)
    ) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.04, 0.25))
    ) +
    
    labs(x = "Count", y = NULL, title = title_str) +
    
    theme_bw(base_size = 14) +
    theme(
      plot.title       = element_text(face  = "bold", hjust  = 0.5,
                                      size  = 16,    color  = "black"),
      axis.text.x      = element_text(face  = "bold", size   = 13,
                                      color = "black"),
      axis.text.y      = element_text(face  = "bold", size   = 12,
                                      color = "black"),
      axis.title.x     = element_text(face  = "bold", size   = 14,
                                      color = "black"),
      axis.title.y     = element_text(face  = "bold", size   = 14,
                                      color = "black"),
      legend.title     = element_text(face  = "bold", size   = 11,
                                      color = "black"),
      legend.text      = element_text(face  = "bold", size   = 10,
                                      color = "black"),
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color     = "grey90",
                                      linewidth = 0.3),
      panel.border     = element_rect(color     = "black",
                                      linewidth = 0.8,
                                      fill      = NA)
    )
}


# STEP 5: Generate + Save all plots

cat("\n=== GENERATING PLOTS ===\n")
cat("Start time:", format(Sys.time(), "%H:%M:%S"), "\n\n")

# ---- GO Biological Process 
cat("1/4  GO Biological Process (", nrow(go_bp_simp), "terms )\n")
p_bp <- make_dotplot(go_bp_simp, "GO Biological Process", top_n = 20)
if (!is.null(p_bp)) save_3fmt(p_bp, "GO_BP_LUNG", width = 11, height = 10)

# ---- GO Cellular Component 
cat("\n2/4  GO Cellular Component (", nrow(go_cc_simp), "terms )\n")
p_cc <- make_dotplot(go_cc_simp, "GO Cellular Component", top_n = 20)
if (!is.null(p_cc)) save_3fmt(p_cc, "GO_CC_LUNG", width = 11, height = 10)

# ---- GO Molecular Function 
cat("\n3/4  GO Molecular Function (", nrow(go_mf_simp), "terms )\n")
p_mf <- make_dotplot(go_mf_simp, "GO Molecular Function", top_n = 20)
if (!is.null(p_mf)) save_3fmt(p_mf, "GO_MF_LUNG", width = 11, height = 10)

# ---- KEGG Pathway 
cat("\n4/4  KEGG Pathway (", nrow(kegg), "terms )\n")
p_kegg <- make_dotplot(kegg, "KEGG Pathway Enrichment", top_n = 20)
if (!is.null(p_kegg)) save_3fmt(p_kegg, "KEGG_LUNG", width = 11, height = 9)

