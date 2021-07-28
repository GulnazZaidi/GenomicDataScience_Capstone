library(SummarizedExperiment)
library(edgeR)
library(annotationDbi)

#import feature counts table and phenotype table
counts <- read.table("Featurecounts.tabular", header = TRUE)
pheno <- read.table("Phenotypes_Capstone.txt", header = TRUE)

#create summarized experiment object from counts and pheno
se <- SummarizedExperiment(assays=list(counts=counts), colData=pheno)
gene_count <- assay(se)

#plot boxplot of gene counts
boxplot(log2(gene_count + 1), 
ylab = "log2 geneCount +1", 
col=as.factor(pheno$Group), 
medcol = "yellow", 
range=0, las =2, 
cex.axis = 0.6, 
main = "Boxplot of log transformed counts")
legend(x = "topright", 
legend = c("Adult", "fetal"), 
fill = c(1,2), 
box.lty = 0, 
cex = 0.5)

# Calculate cpm 
cpm_gene_count <- cpm(gene_count)

#plot boxplot for cpm+0.1 of gene counts
boxplot(log2(cpm_gene_count + 0.1), 
ylab = "log2 cpm_geneCount + 0.1",  
col=as.factor(pheno$Group), 
medcol = "yellow", range=0, las =2, 
cex.axis = 0.6, 
main = "Boxplot of log transformed cpm")
legend(x = "topright", 
legend = c("Adult", "fetal"), 
fill = c(1,2), 
box.lty = 0, 
cex = 0.5)

#for PCA
#remove lowly expressed genes for cpm < 0.5
#cpm of 0.5 corresponds to 10-12 reads. rowSums of (cpm>0.5) > 1 means atleast 1 sample in the row shows exp higher than cpm of 0.5

keep <- rowSums(cpm_gene_count > 0.5) >= 1
keep_gene_count <- gene_count[keep,]
dim(keep_gene_count)
dim(gene_count)

#log of filtered gene cpms
log_keep <- log2(keep_gene_count + 0.1)

#PCA of log transformed genes
pca_gene <- prcomp(log_keep, center = TRUE, scale. = TRUE)
pca_Rot <- data.frame(pca_gene$rotation)

#plot PC1 and PC2
plot(pca_Rot$PC1, pca_Rot$PC2, col= as.factor(pheno$Group), 
xlab = "PC1", ylab = "PC2", 
main = "Scatterplot of PC1 and PC2", pch=19)
legend(x = "topright", 
legend = c("Adult", "fetal"), 
fill = c(1,2), 
box.lty = 0, 
cex = 0.5)





