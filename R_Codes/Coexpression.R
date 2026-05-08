#Now run WGCNA:
library(WGCNA)


# Extract expression matrix
expr_matrix <- assay(vsd)

# Transpose for WGCNA (samples = rows, genes = columns)
datExpr <- t(expr_matrix)

# Convert to data frame
datExpr <- as.data.frame(datExpr)

#: Check for Bad Genes / Samples
gsg <- goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}


powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers)
sft.data <- sft$fitIndices


library(ggplot2)
plot1 <- ggplot(sft.data, aes(Power,SFT.R.sq, label = Power)) +
  geom_point()+
  geom_text(nudge_y =0.1)+
  geom_hline(yintercept = 0.8, color = 'red')+
  labs( x = 'Power', y = "free topology model fit(signed R^2)")+
  theme_classic()

plot2 <- ggplot(sft.data, aes(Power,mean.k., label = Power)) +
  geom_point()+
  geom_text(nudge_y = 0.1)+
  labs( x = 'Power', y = 'mean connectivity')+
  theme_classic()

library(gridExtra)
grid.arrange(plot1, plot2, nrow= 2)

#Hierarchial clustering
sft_power <- 9
Temp_cor <- cor
cor <- WGCNA::cor

net <- blockwiseModules(
  datExpr,
  maxBlockSize = 10000,
  power = sft_power,
  TOMType = "signed",
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  randomSeed = 1234,
  verbose = 3)

cor <- Temp_cor
moduleColors <- net$colors

# Module eigengenes
MEs <- net$MEs
# Order module eigengenes for consistency
MEs <- orderMEs(MEs)
table(net$colors)

#get no. of gene for each module
library(pheatmap)

#Dendrogram
block <- 1
blockGenes <- net$blockGenes[[block]]
plotDendroAndColors(net$dendrograms[[block]],
                    cbind(
                      net$unmergedColors[blockGenes],
                      net$colors[blockGenes]),
                    c("unmerged", "merged",
                      "module colors",
                      dendroLabels = FALSE,
                      addGuide = TRUE,
                      hang = 0.03,
                      guideHang = 0.05)
)

# Relate module trait associations

#import data
my_data <- read.csv("/home/user12345/Documents/BHU_Dessertation/clinical_data.csv")
library(dplyr)

my_data <- my_data %>%
  dplyr::select(-starts_with("data."))


traits <- data.frame(
  PPMS = ifelse(my_data$classification_of_ms == "PPMS", 1, 0),
  SPMS = ifelse(my_data$classification_of_ms == "SPMS", 1, 0),
  RRMS = ifelse(my_data$classification_of_ms == "RRMS", 1, 0),
  Control = ifelse(my_data$classification_of_ms == "Non-MS control", 1, 0)
)

rownames(traits) <- rownames(MEs)
all(rownames(traits) == rownames(MEs))
all(rownames(datExpr) == rownames(traits))


#Compute P-values
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,
                                      nSamples = nrow(datExpr))


# Define numbers of genes and samples
nSamples <- nrow(datExpr)
nGenes <- ncol(datExpr)

traits$Disease <- as.numeric(
  rowSums(traits[, c("PPMS", "SPMS", "RRMS")]) > 0
)

#check module–trait correlation:
moduleTraitCor <- cor(MEs,
                      traits,
                      use = "pairwise.complete.obs",
                      method = "pearson")


# visualize module-trait association as a heatmap
heatmap.data <- merge(MEs, traits, by = 'row.names')

head(heatmap.data)
library(tibble)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

moduleTraitCor <- cor(MEs, traits, use = "p")

pdf("Module_Trait_Heatmap.pdf", width = 14, height = 10)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  zlim = c(-1, 1),
  cex.text = 0.6,
  main = "Module-Trait Relationships"  # normal hyphen
)

dev.off()

module.gene.mapping <- data.frame(
  gene   = colnames(datExpr),
  module = net$colors,
  stringsAsFactors = FALSE
)

# Extract genes of one module (turquoise) 
turquoise_genes <- module.gene.mapping$gene[
  module.gene.mapping$module == "turquoise"
]

#Create list of genes per module
module.gene.table <- data.frame(
  gene   = colnames(datExpr),
  module = net$colors,
  stringsAsFactors = FALSE
)

write.table(turquoise_genes,
            file = "turquoise",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
