#We have to make a new counts table for RNA

BiocManager::install("DESeq2")

#Load (this must be done every time)
library(DESeq2)

#Load additional libraries that might be useful (these may need to be installed)
library(ggplot2)

counts_file<- read.csv("e910_raw_counts.csv", row.names = 1)
c2<- read.csv("e12_raw_counts.csv", row.names = 1)
c3<- read.csv("e14_raw_counts.csv", row.names = 1)
c4<- read.csv("e16_raw_counts.csv", row.names = 1)


counts_file<- cbind(counts_file, c2[,4:7])
counts_file<- cbind(counts_file, c3[,4:6])
counts_file<- cbind(counts_file, c4[,4:7])

row.names(counts_file)<- make.names(counts_file$`symbol`, unique = TRUE)

counts_protein_coding <- counts_file

metadata<- read.csv("metadata.csv", row.names = 1)
all(names(counts_protein_coding) %in% rownames(metadata))

#DDS by "exp"
dds<- DESeqDataSetFromMatrix(countData = counts_protein_coding, colData = metadata, design = ~ timepoint )

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds<- estimateSizeFactors(dds)
sizeFactors(dds)
combined_normalized_counts<- counts(dds, normalized =TRUE)
combined_normalized_counts<- as.data.frame(combined_normalized_counts)
write.table(combined_normalized_counts, file= "dds_combined_normalized_counts.txt", sep = "\t", quote = FALSE)

#Change Ensembl to normal gene names
write.table(row.names(combined_normalized_counts), file = "ensemb.txt", row.names = F, quote = F)
genes<- read.table("ensemb.txt", sep = "\t")
combined_normalized_counts$ensembl<- rownames(combined_normalized_counts)
colnames(genes)<- c("ensembl", "symbol")
x<- merge(combined_normalized_counts,genes, by = "ensembl")
combined_normalized_counts<- x
combined_normalized_counts<- combined_normalized_counts[,-1]
row.names(combined_normalized_counts)<- make.names(combined_normalized_counts$`symbol`, unique = TRUE)
combined_normalized_counts<- combined_normalized_counts[,-19]

write.table(rnavalues, file= "dds_combined_normalized_counts_proteingenes.txt", sep = "\t", quote = FALSE)

