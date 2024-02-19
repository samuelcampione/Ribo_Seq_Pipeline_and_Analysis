
library(GenomicFeatures)
library(DESeq2)
library(tidyverse)
library(limma)
library(edgeR)
library(sva)


# assumes samples are in same directory and named like "SRA123456_ribo_cond_sample_batch_series.tsv"

####################################
# GIR2 Knockout vs WT
####################################

# Set the working directory to where the count files are located
setwd("/Users/scampione/Projects/Buck_Institute/Ribo_Profiling/PMID34900236/Proof_of_Prinicpal")
count_files <- list.files(pattern = "\\.tsv$")


# Processing sample metadata
sample_info <- data.frame(
  sample = gsub("_b[[:digit:]]+_s[[:digit:]]+tsv", "", count_files),
  condition = gsub("SRR[[:digit:]]+_ribo_(.*)_[[:digit:]]_b[[:digit:]]+_s[[:digit:]]+\\.tsv", "\\1", count_files),
  batch = str_extract(count_files, "_b[[:digit:]]+"),
  stringsAsFactors = FALSE
) %>% mutate(across(everything(), as.factor))

col_names <- c("gene_id", "count")


# Read the files and create a list of data frames, capturing potential problems
count_list <- lapply(count_files, function(f) {
  read_tsv(f, 
           col_names = col_names, 
           col_types = cols(
             gene_id = col_character(),
             count = col_double()), 
           skip = 4,
           col_select = c("gene_id", "count"))
})


# Extract gene identifiers from the first file
genes <- count_list[[1]] %>% pull(gene_id)


# Validate that all files have the same gene identifiers
for(df in count_list) {
  if(!all(genes == df %>% pull(gene_id))) {
    stop("Gene identifiers do not match across all files.")
  }
}


# Combine count data into a single matrix
count_matrix <- do.call(cbind, lapply(count_list, function(df) df %>% pull(count)))
rownames(count_matrix) <- genes # row names as gene names
colnames(count_matrix) <- factor(sample_info$sample)


# Create a DGEList object (edgeR)
y <- DGEList(counts=count_matrix, group=sample_info$condition)


# Filter for only genes with a minimum of 10 in at least 4 samples
cpm_values <- cpm(y) 
keep <- rowSums(cpm_values >= 10) >=4 
y <- y[keep, ] 


# Batch Correction with ComBat-Seq
log_counts <- cpm(y, log=TRUE, prior.count=2) # log the counts
batch_corrected_log_counts <- ComBat_seq(log_counts, batch=sample_info$batch)
batch_corrected_counts <- exp(batch_corrected_log_counts) 


# UpdateDGEList object with batch-corrected counts
y$counts <- batch_corrected_counts



####################################
# Differential translation analysis
####################################

# Ensure WT is control and gir2 knockout is treatment
sample_info$condition <- relevel(factor(sample_info$condition), ref = "wt")


# Create a design matrix
design <- model.matrix(~ condition, data=sample_info) 


# Estimate dispersion
y <- estimateGLMRobustDisp(y, design) 
#y <- estimateDisp(y, design) 

# Fit the model
fit <- glmQLFit(y, design) 
qlf <- glmQLFTest(fit, coef=2) # Run the likelihood ratio tests

# Extract differentially expressed genes
results <- topTags(qlf, n=Inf)
sig_genes <- results$table[results$table$FDR < 0.05, ]


# Upregulated genes
sig_genes["YDL194W",]$logFC #aka SNF3
sig_genes["YGL184C",]$logFC #aka STR3

# Downregulated genes
sig_genes["YNL011C",]$logFC 
sig_genes["snR86",]$logFC 

# Non-differentially expressed genes
# Is SPO24 differentially expressed?
spo24 <- "YPR036W-A" %in% sig_genes
spo24
 
# Is tma10 differentially expressed?
tma10 <- "YLR327C" %in% sig_genes
tma10



summary(decideTests(qlf))

plotMD(qlf, main="Differential Translation: gir2Δ vs WT")
abline(h=c(-1,1), col="blue")

annotate_genes <- c("YDL194W", "YGL184C", "YNL011C", "snR86", "YPR036W-A", "YLR327C", "YDR152W")
annotations <- c("SNF3", "STR3", "YNL011C", "snR86", "SPO24", "TMA10", "gir2")

results$table["YLR327C",]

logFC <- qlf$table[annotate_genes, "logFC"]
logCPM <- qlf$table[annotate_genes, "logCPM"]

cex_size <- 1.4  # Adjust as needed for your plot
text_color <- "blue"
border_color <- "blue"
plotting_character <- 13
text_font <- 2

text(logCPM, logFC, labels=annotations, pos=3, cex=cex_size, col=text_color,
     pch=plotting_character, bg="gray", font=text_font)
points(logCPM, logFC, pch=plotting_character, bg="white", col=border_color, cex=cex_size)



