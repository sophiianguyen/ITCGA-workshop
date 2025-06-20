#''''''''''''''''''''''''''''''''''''
#' Read Count Analysis with edgeR
#' @author Cooper Kimball-Rhines
#' @date 2025-06-10
#''''''''''''''''''''''''''''''''''''

# Install package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
install.packages("tidyverse")
install.packages("ggrepel")


# Load library
library(edgeR)
library(tidyverse)
library(ggrepel)

# Set working directory manually: Session -> Set Working Directory -> To Source File Location
# Load data from counts.txt file
counts <- read_tsv("~/Desktop/Job/Research/counts.txt", skip = 1) %>%
  rename_with(~str_remove(., '/itcgastorage/share_home/kenneth.chen001/project/results/bam')) %>%
  rename_with(~str_remove(., '.bam'))

head(counts)

# Provide treatment metadata
condition <- c("A", "A", "A", "P", "P", "P")
dim(counts)

# Filter out genes expressed below cutoff
totalexp <- rowSums(counts[,7:12])
hist(totalexp)

counts <- filter(counts, totalexp > 10)

# Move annotation info to separate object
ann <- counts[,1:6]

# Combine gene counts and metadata into edgeR object
d <- DGEList(counts = counts[,7:12],
             group = factor(condition),
             genes = counts[,1:6])
str(d)
dim(d) # 12,558 genes across 9 samples

# Normalize the data
d <- calcNormFactors(d, method = "TMM")
d$samples

# Plot MDS
samples <- c("A1","A2", "A3", "P1", "P2", "P3")

plotMDS(d, col = as.numeric(d$samples$group), label = d$samples$group)

# Generate design matrix
design <- model.matrix(~ 0 + d$samples$group)
design
colnames(design) <- levels(d$samples$group)

# Estimate dispersion
dg <- estimateGLMCommonDisp(d, design)
plotBCV(dg)

# Fit the model
fit <- glmFit(dg, design)

# This compares group 1 (Cs, 1) to group 2 (Ts, -1), ignoring group 3 (Vs, 0)
# and does a likelihood ratio test for each gene
fitCT <- glmLRT(fit, contrast=c(1, -1))

# Sort out the differentially expressed genes
tabCT <- topTags(fitCT,n=Inf,adjust.method="BH", sort.by = "PValue")$table

# Make a significance column
tabCT <- tabCT %>% 
  mutate(significance = case_when((FDR < 0.05 & logFC > 1) ~ "Upregulated", (FDR < 0.05 & logFC < -1) ~ "Downregulated", .default = "Not significant"))

# Pull out the top differentially expressed genes
top_genes <- filter(tabCT, -log10(PValue) > 5)

# Visualize the expression data!
volcano_plot <- ggplot(tabCT, aes(x = logFC, y = -log10(PValue), color = significance)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Not significant" = "grey","Upregulated" = "red", "Downregulated" = "blue")) +
  geom_text_repel(data = top_genes, aes(label = Geneid), size = 3.5, fontface = 'bold') +
  geom_hline(yintercept = 5, linetype = "dashed") +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  theme_classic() +
  theme(legend.position = 'None',
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank()) +
  labs(x = "Log2 Fold Change", y = "-Log10 p-value")

volcano_plot

# Write out the most significantly differential genes
tabCT |>
  filter(significance != "Not significant") |>
  write_tsv(file = "topGenes.tsv")

# Let's investigate the top differentially expressed gene
topCount <- counts |>
  filter(Geneid == "ENSG00000041982") |>
  select(-c(Geneid, Chr, Start, End, Strand,Length)) |>
  t()

topGene <- data.frame(cbind(as.character(d$samples$group), topCount))
colnames(topGene) <- c("Treatment","Reads")
topGene$Reads <- as.numeric(topGene$Reads)

genePlot <- ggplot(data = topGene, 
                   mapping = aes(x = Treatment, y = Reads, color = Treatment)) +
  geom_jitter(width = 0.25) +
  labs(y = "Read counts at HSPA6") +
  theme_bw()

genePlot # We can see that the counts of this gene are way higher in the TGF-b treatment

