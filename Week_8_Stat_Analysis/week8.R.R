library(SummarizedExperiment)
library(AnnotationDbi)
library(tidyverse)
library(org.Hs.eg.db)
library(DESeq2)

#load the featureCount result and phenotype table
counts <- as.matrix(read.table("Featurecounts.tabular", header = TRUE))
pheno <- read.table("Phenotypes_Capstone.txt", header = TRUE)

#create summarized experiment
se <- SummarizedExperiment(assays=list(counts=counts), colData=pheno)
seDESeq <- DESeqDataSet(se, design= ~ Group)
seDESeq <- DESeq(seDESeq)
result_DES <- as.data.frame(results(seDESeq,name="Group_Fetal_vs_Adult")) %>% rownames_to_column("ENTREZID")
result_ord <- result_DES[order(result_DES$pvalue),]

#make dataframe for annotation and final result
anno <- AnnotationDbi::select(org.Hs.eg.db,keys=result_ord$ENTREZID, columns=c("SYMBOL","GENENAME"),keytype="ENTREZID")

result_final <- data.frame(GeneName = anno$GENENAME, FoldChange = result_ord$log2FoldChange, Pvalue = result_ord$pvalue, Padj = result_ord$padj)

#print result to tab delimited and csv file
write.table(result_final,file ="Fetal_vs_Adult_DESeq_annotated.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
write.csv(result_final,file="Fetal_vs_Adult_DESeq_annotated.csv",row.names=FALSE)

#for volcano plot
#add filtering criteria, sort, annotate with gene symbols
res01_lfc <- as.data.frame(results(seDESeq, alpha=0.05, lfcThreshold = 0.58))%>% rownames_to_column("ENTREZID")
res01_lfc <- res01_lfc[order(res01_lfc$pvalue),]
res01_lfc <- bind_cols(res01_lfc, Symbol=anno$SYMBOL[match(res01_lfc$ENTREZID, anno$ENTREZID)])
res01_lfc$top10 <- " "
res01_lfc$top10[1:10] <- as.character(res01_lfc$Symbol[1:10])

#volcano plot with top10 DE genes labeled
ggplot(res01_lfc, aes(x = log2FoldChange, y = -log10(padj))) + 
geom_point() + 
geom_text(aes(label= top10, nudge_x= 0.25, nudge_y= 0.25, check_overlap= T )) +
ggtitle("Volcano Plot of Fetal vs Adult") +
xlab("log2 fold change") + 
ylab("-log10 adjusted p-value")
