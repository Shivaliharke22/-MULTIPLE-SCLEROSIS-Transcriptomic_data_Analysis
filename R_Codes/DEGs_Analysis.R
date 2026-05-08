###DGE Analyis######### 
# Convert to matrix for normalization
data_filtered_matrix <- as.matrix(data_filtered[, 2:95])
data_filtered_df <- as.data.frame(data_filtered)
rownames(data_filtered_df) <- make.unique(data_filtered$SYMBOL)

sample_info <- data.frame(
  row.names = c("GSM8813321",  "GSM8813322", "GSM8813323",  "GSM8813324",  "GSM8813325",  "GSM8813326",  "GSM8813327",  "GSM8813328" , "GSM8813329" ,
                "GSM8813330",  "GSM8813331", "GSM8813332",  "GSM8813333",  "GSM8813334",  "GSM8813335" , "GSM8813336",  "GSM8813337" , "GSM8813338",  "GSM8813339" , "GSM8813340" ,
                "GSM8813341", "GSM8813342",  "GSM8813343",  "GSM8813344" , "GSM8813345",  "GSM8813346" , "GSM8813347",  "GSM8813348" , "GSM8813349" , "GSM8813350" , "GSM8813351" ,
                "GSM8813352",  "GSM8813353",  "GSM8813354",  "GSM8813355" , "GSM8813356",  "GSM8813357",  "GSM8813358",  "GSM8813359" , "GSM8813360" , "GSM8813361",  "GSM8813362" ,
                "GSM8813363",  "GSM8813364",  "GSM8813365",  "GSM8813366" , "GSM8813367",  "GSM8813368",  "GSM8813369",  "GSM8813370",  "GSM8813371" , "GSM8813372" , "GSM8813373" ,
                "GSM8813374",  "GSM8813375",  "GSM8813376", "GSM8813377" , "GSM8813378",  "GSM8813379" , "GSM8813380",  "GSM8813381" , "GSM8813382",  "GSM8813383" , "GSM8813384" ,
                "GSM8813385",  "GSM8813386",  "GSM8813387",  "GSM8813388",  "GSM8813389",  "GSM8813390",  "GSM8813391",  "GSM8813392",  "GSM8813393",  "GSM8813394" , "GSM8813395" ,
                "GSM8813396",  "GSM8813397",  "GSM8813398","GSM8816432","GSM8816433","GSM8816434","GSM8816435", "GSM8816436","GSM8816437","GSM8816438","GSM8816439","GSM8816440",
                "GSM8816441","GSM8816442","GSM8816443","GSM8816444","GSM8816445","GSM8816446","GSM8816447"),
  
  Condition = c("MS", "control","MS","MS","MS","MS","control","control","MS","control","MS","control",
                "control","control","control","control","MS", "control","control","control","MS","MS", "control",
                "MS","MS","MS","control","MS","control","MS","MS","control","MS","MS", "control","MS","MS","MS","MS",
                "control","MS","MS","MS","MS","MS","control","MS","MS","MS","control","control", "MS","MS","MS","MS","MS","MS",
                "MS","control","control", "control","MS","control","control","MS","MS","MS","MS","MS","MS","MS","MS",
                "control","control","MS","MS","control","control","control","control","control","MS","MS","MS","MS","MS",
                "MS","MS","MS","MS","MS","MS","MS","control"))



# Assume first column = gene IDs, rest = counts
gene_symbol <- data_filtered_df[, 1]
# Convert to matrix
filtered_matrix_int <- as.matrix(data_filtered_matrix)
# 3. Force storage mode to integer (important for DESeq2/edgeR)
mode(filtered_matrix_int) <- "integer"
# 4. Add gene IDs as rownames
rownames(filtered_matrix_int) <- gene_symbol

class(filtered_matrix_int)

filtered_matrix_df <- as.data.frame(filtered_matrix_int)

# Compute total expression per gene require dataframe
filtered_matrix_df$Total <- rowSums(filtered_matrix_df)

# Or if gene symbols are in rownames (not a column)
filtered_matrix_int <- filtered_matrix_int[!is.na(rownames(filtered_matrix_int)), ]
filtered_matrix_int <- filtered_matrix_int[order(rowSums(filtered_matrix_int), decreasing = TRUE), ]
filtered_matrix_int <- filtered_matrix_int[!duplicated(rownames(filtered_matrix_int)), ]

# Check unique genes now
length(unique(rownames(filtered_matrix_int))) == nrow(filtered_matrix_int)

library(DESeq2)

sample_info$Condition <- as.factor(sample_info$Condition)
DGE_Data <- DESeqDataSetFromMatrix(countData = filtered_matrix_int,
                                   colData = sample_info,
                                   design = ~ Condition)
#Run DESeq2 Analysis
DGE_Data <- DESeq(DGE_Data) 


vsd <- vst(DGE_Data)
plotPCA(vsd, intgroup="Condition")

# Variance stabilized data (recommended)
boxplot(assay(vsd),
        outline = FALSE,   # removes dots
        las = 2,
        main = "VST Normalized Data",
        col = "lightblue")


#Create the results object
DGE_result <- results(DGE_Data)

result_clean <- DGE_result[!is.na(DGE_result$padj), ]
sum(is.na(result_clean$`padj`))  

# Set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1  # means 2-fold change

all_DEG <- subset(
  result_clean,
  padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff 
)

write.csv(all_DEG, "DESeq2_all_genes.csv", row.names = TRUE) #179

# Convert to dataframe for easy handling
result_clean <- as.data.frame(DGE_result)

# Add regulation column
result_clean$regulation <- "NotSig"
result_clean$regulation[result_clean$padj < padj_cutoff & result_clean$log2FoldChange > log2fc_cutoff] <- "Up"
result_clean$regulation[result_clean$padj < padj_cutoff & result_clean$log2FoldChange < -log2fc_cutoff] <- "Down"


# Upregulated genes
up_genes <- subset(result_clean, padj < padj_cutoff & log2FoldChange > log2fc_cutoff)
up_genes <- subset(result_clean, regulation == "Up") #116
# Downregulated genes
down_genes <- subset(result_clean, padj < padj_cutoff & log2FoldChange < -log2fc_cutoff)
#Extract all rows where regulation column is "down"
down_genes <- subset(result_clean, regulation == "Down") #63

combined_df_up <- rbind(up_genes, down_genes)

######VOLCANO PLOT############

# VOLCANO Plots
ggplot(result_clean, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
  labs(title = "DEG_RESULT",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value") +
  theme_minimal()

