library(AnnotationHub)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)

# rework the DESEQ2 results to get list of DE genes
result_DES <- as.data.frame(results(seDESeq,name="Group_Fetal_vs_Adult", alpha = 0.01, lfcThreshold = 1)) %>% rownames_to_column("ENTREZID")
result_DES <- result_DES[order(result_DES$pvalue),]
ann_df <- AnnotationDbi::select(org.Hs.eg.db,keys=result_DES$ENTREZID, columns=c("SYMBOL","GENENAME"),keytype="ENTREZID")
result_DES <- bind_cols(result_DES, Symbol=ann_df$SYMBOL[match(result_DES$ENTREZID, anno$ENTREZID)])
res.na <- na.omit(result_DES)
result_DES_1 = result_DES[result_DES$padj < 0.01,]
result_DES_2 <- result_DES_1[abs(result_DES_1$log2FoldChange) >= 1,]

#extract the up- and down-regulated genes list
de_up <- result_DES_2$ENTREZID[result_DES_2$log2FoldChange>0]
de_up <- de_up[!is.na(de_up)]
de_down <- result_DES_2$Symbol[result_DES_2$log2FoldChange<0]
de_down <- de_down[!is.na(de_down)]

#extract the promoter regions in Upregulated & downregulated genes
UCSC_promoters = promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, columns=c("gene_id"))
DE_promUp = UCSC_promoters[as.character(UCSC_promoters$gene_id) %in% de_up,]
length(DE_promUp)
DE_promDown = UCSC_promoters[as.character(UCSC_promoters$gene_id) %in% de_down,]
length(DE_promDown)

#query the Annotation hub cell types and Refseq
ah <- AnnotationHub()
qLiver <- query(ah, c("H3K4me3", "E066"))
qBrainF <- query(ah, c("H3K4me3", "E081"))
qBrainA <- query(ah, c("H3K4me3", "E073"))
qSeq <- query(ah, "RefSeq")

#download the narrowPeak file for liver, fetal brain and adult brain
grN_liver <- qLiver[[2]]
grN_BrainF <- qBrainF[[2]]
grN_BrainA <- qBrainA[[2]]
#check the width of the peaks
summary(width(grN_liver))

#check overlaps between Up/down promoters & narrowPeaks of Fetal, adult, liver
BrainA_Up <- subsetByOverlaps(grN_BrainA, DE_promUp)
length(BrainA_Up)
BrainA_Down <- subsetByOverlaps(grN_BrainA, DE_promDown)
length(BrainA_Down)
BrainF_Up <- subsetByOverlaps(grN_BrainF, DE_promUp)
length(BrainF_Up)
BrainF_Down <- subsetByOverlaps(grN_BrainF, DE_promDown)
length(BrainF_Down)
liver_Up <- subsetByOverlaps(grN_liver, DE_promUp)
length(liver_Up)
liver_Down <- subsetByOverlaps(grN_liver, DE_promDown)
length(liver_Down)
