library(dbplyr)
library(biomaRt)
library(tidyverse)
library(GEOquery)
library(data.table)   # fast file reading
library(stringr) 

#### 1. Set folder path ####
folder_path2 <- "/home/user12345/Documents/BHU_Dessertation/DATA"
#Match all .genes.results.gz files
files2 <- list.files(folder_path2, pattern = ".gz$", full.names = TRUE)


data_list <- lapply(files2, function(file) {
  dat <- fread(file)   # fread works with .gz directly
  
  # Extract patient ID from filename (remove extension)
  patient_id <- str_remove(basename(file), "\\..*")
  
  # Check column names and adjust if needed
  if (!("gene_id" %in% colnames(dat))) {
    stop(paste("Column 'gene_id' not found in file:", file))
  }
  if (!("expected_count" %in% colnames(dat))) {
    stop(paste("Column 'expected_count' not found in file:", file))
  }
  # Keep only gene_id and expected_count columns
  dat<- dat[, .(gene_id, expected_count)]
  setnames(dat, "expected_count", patient_id)  # rename expected_count column with patient ID
  
  return(dat)
})

##### 3. Merge all patients by gene_id ##
raw_data <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), data_list)  #60671

#### 4. Convert to expression matrix ####
rownames(raw_data) <- raw_data$gene_id
expr_matrix <- as.matrix(raw_data[, -1, with = FALSE])

#### 5. Quick check ####
dim(expr_matrix)         # genes x patients   60671    94
expr_matrix[1:5, 1:5]    # first 5 genes × first 5 patients

library(org.Hs.eg.db)
library(AnnotationDbi)

EnsembleIds <- raw_data$gene_id

# Map Ensembl IDs to HGNC symbols
hgnc.list <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = EnsembleIds,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENSEMBL")
)

raw_symbol <- dplyr::left_join(
  as.data.frame(raw_data),
  as.data.frame(hgnc.list),
  by = c("gene_id" = "ENSEMBL")
) 

# Reorder columns: gene_id | SYMBOL | rest
raw_symbol <- raw_symbol[,c("gene_id", "SYMBOL", setdiff(names(raw_symbol), c("gene_id", "gene_symbol")))]
raw_symbol <- raw_symbol[, !names(raw_symbol) %in% "SYMBOL.1"]


head(raw_symbol)
sum(is.na(raw_symbol)) #23911

# Remove rows with NA in gene symbol column
raw_symbol_noNA <- raw_symbol[!is.na(raw_symbol$SYMBOL), ]

library(janitor)
# Convert all columns except the first two to numeric
# Create a copy of your data frame to work on
data_numeric <- raw_symbol_noNA

# 2 Get the correct column range (from 3rd column to last)
numeric_cols <- 3:ncol(data_numeric)

# 3 Convert those columns to numeric
data_numeric[, numeric_cols] <- lapply(
  raw_symbol_noNA[, numeric_cols],
  function(x) as.numeric(as.character(x))
)

str(data_numeric)

# Check dimensions before and after
dim(raw_symbol) #  64785    96
dim(data_numeric) # 40874    96

# Verify that no NA remains in gene_id column
sum(is.na(data_numeric))

dup <- get_dupes(data_numeric,gene_id) 

raw_symbol_unique <- data_numeric %>%
  as.data.frame() %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

dim(raw_symbol_unique)# 37550    96

#rename columns
colnames(raw_symbol_unique)[3:96] <- c(
  "GSM8813321", "GSM8813322", "GSM8813323", "GSM8813324", "GSM8813325", "GSM8813326",
  "GSM8813327", "GSM8813328", "GSM8813329", "GSM8813330", "GSM8813331", "GSM8813332",
  "GSM8813333", "GSM8813334", "GSM8813335", "GSM8813336", "GSM8813337", "GSM8813338",
  "GSM8813339", "GSM8813340", "GSM8813341", "GSM8813342", "GSM8813343", "GSM8813344",
  "GSM8813345", "GSM8813346", "GSM8813347", "GSM8813348", "GSM8813349", "GSM8813350",
  "GSM8813351", "GSM8813352", "GSM8813353", "GSM8813354", "GSM8813355", "GSM8813356",
  "GSM8813357", "GSM8813358", "GSM8813359", "GSM8813360", "GSM8813361", "GSM8813362",
  "GSM8813363", "GSM8813364", "GSM8813365", "GSM8813366", "GSM8813367", "GSM8813368",
  "GSM8813369", "GSM8813370", "GSM8813371", "GSM8813372", "GSM8813373", "GSM8813374",
  "GSM8813375", "GSM8813376", "GSM8813377", "GSM8813378", "GSM8813379", "GSM8813380",
  "GSM8813381", "GSM8813382", "GSM8813383", "GSM8813384", "GSM8813385", "GSM8813386",
  "GSM8813387", "GSM8813388", "GSM8813389", "GSM8813390", "GSM8813391", "GSM8813392",
  "GSM8813393", "GSM8813394", "GSM8813395", "GSM8813396", "GSM8813397", "GSM8813398",
  "GSM8816432", "GSM8816433", "GSM8816434", "GSM8816435", "GSM8816436", "GSM8816437",
  "GSM8816438", "GSM8816439", "GSM8816440", "GSM8816441", "GSM8816442", "GSM8816443",
  "GSM8816444", "GSM8816445", "GSM8816446", "GSM8816447"
)


# Convert to matrix 
raw_matrix <- as.matrix(raw_symbol_unique[, 3:96])

# Convert matrix to data frame first
raw_matrix_df <- as.data.frame(raw_matrix)

# Add SYMBOL column safely
raw_matrix_df$SYMBOL <- raw_symbol_unique$SYMBOL
raw_matrix_df <- raw_matrix_df[, c("SYMBOL", setdiff(names(raw_matrix_df), "SYMBOL"))]

numeric_data <- raw_matrix_df[, -(1)]

# Calculate fraction of zeros per gene (row)
zero_fraction <- rowMeans(numeric_data == 0)

# Define threshold (30%)
threshold <- 0.30

# Keep only genes where less than 30% of patients have zero expression
data_filtered <- raw_matrix_df[zero_fraction < threshold, ]
# Check results
dim(raw_symbol_unique) # 37550    96
dim(raw_matrix_df) # 37550    95
dim(data_filtered) # [1] 22985    95















                   
