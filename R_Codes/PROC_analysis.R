#PRoc Analysis
PRoc_data <- t(data_filtered)
#Convert matrix → data frame
PRoc_data <- as.data.frame(PRoc_data, stringsAsFactors = FALSE)

#Use first row as gene names
# Extract gene names
Gene_name <- as.character(PRoc_data[1, ])

# Remove first row
PRoc_data <- PRoc_data[-1, ]

# Assign gene names as column names
colnames(PRoc_data) <- Gene_name

#Add SampleID column from rownames
PRoc_data$SampleID <- rownames(PRoc_data)

#Convert expression values to numeric
gene_cols <- setdiff(colnames(PRoc_data), "SampleID")

PRoc_data[gene_cols] <- lapply(PRoc_data[gene_cols], function(x) {
  as.numeric(trimws(x))
})


#Convert condition to binary (REQUIRED for pROC)
sample_info$SampleID <- rownames(PRoc_data)

sample_info$Condition_binary <- ifelse(
  sample_info$Condition == "MS", 1, 0
)

#Merge condition into proc_data
PRoc_data <- merge(
  PRoc_data,
  sample_info[, c("SampleID", "Condition_binary")],
  by = "SampleID"
)

# Put SampleID back as rownames
rownames(PRoc_data) <- PRoc_data$SampleID
PRoc_data$SampleID <- NULL

#Define hubgene
hub_genes <- c(
  "VCAM1",
  "ICAM1",
  "CXCL8",
  "IL6",
  "CXCL3",
  "IL1B"
)

#Run ROC analysis for ALL hub genes
library(pROC)

roc_list <- list()
auc_results <- data.frame(
  Gene = character(),
  AUC  = numeric(),
  stringsAsFactors = FALSE
)


for (gene in hub_genes) {
  
  # Skip gene if not present
  if (!gene %in% colnames(PRoc_data)) {
    warning(paste("Gene not found:", gene))
    next
  }
  
  roc_obj <- roc(
    response  = PRoc_data$Condition_binary,
    predictor = PRoc_data[[gene]],
    quiet = TRUE
  )
  
  roc_list[[gene]] <- roc_obj
  
  auc_results <- rbind(
    auc_results,
    data.frame(
      Gene = gene,
      AUC  = as.numeric(auc(roc_obj))
    )
  )
}

#Plot ROC

dir.create("ROC_Plots", showWarnings = FALSE)

for (gene in hub_genes) {
  
  if (!gene %in% colnames(PRoc_data)) next
  
  roc_obj <- roc(
    response  = PRoc_data$Condition_binary,
    predictor = PRoc_data[[gene]],
    quiet = TRUE
  )
  
  png(
    filename = paste0("ROC_Plots/ROC_", gene, ".png"),
    width = 800,
    height = 800,
    res = 150
  )
  
  # ---- Plot ROC ----
  plot(
    roc_obj,
    col = "darkred",
    lwd = 2,
    legacy.axes = TRUE,
    xlim = c(1, 0),
    ylim = c(0, 1),
    xlab = "1 − Specificity",
    ylab = "Sensitivity",
    main = paste(
      "ROC Curve for", gene,
      "\nAUC =", round(as.numeric(auc(roc_obj)), 3)
    )
  )
  
  
  # ---- Legend (INSIDE png device) ----
  legend(
    "bottomright",
    legend = c(
      paste("Empirical ROC (AUC =", round(as.numeric(auc(roc_obj)), 3), ")"),
      "Random classifier"
    ),
    col  = c("darkred", "gray50"),
    lwd  = c(2, 2),
    lty  = c(1, 2),
    cex  = 0.9,
    bty  = "n"
  )
  
  dev.off()
}
