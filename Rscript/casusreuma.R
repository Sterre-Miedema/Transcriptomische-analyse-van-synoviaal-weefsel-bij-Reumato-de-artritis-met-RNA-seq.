setwd("/Users/sterremiedema/OneDrive - NHL Stenden/transcriptomics/")
install.packages("BiocManager")
BiocManager::install('Rsubread')
library(Rsubread)
buildindex(
  basename = 'ref_humaan',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)

align.nor1 <- align(index = "ref_humaan", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2 = "SRR4785819_2_subset40k.fastq" , output_file = "nor1.BAM")
align.nor2 <- align(index = "ref_humaan", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2 = "SRR4785820_2_subset40k.fastq" , output_file = "nor2.BAM")
align.nor3 <- align(index = "ref_humaan", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2 = "SRR4785828_2_subset40k.fastq" , output_file = "nor3.BAM")
align.nor4 <- align(index = "ref_humaan", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2 = "SRR4785831_2_subset40k.fastq" , output_file = "nor4.BAM")
align.ra1 <- align(index = "ref_humaan", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2 = "SRR4785979_2_subset40k.fastq" , output_file = "ra1.BAM")
align.ra2 <- align(index = "ref_humaan", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2 = "SRR4785980_2_subset40k.fastq" , output_file = "ra2.BAM")
align.ra3 <- align(index = "ref_humaan", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2 = "SRR4785986_2_subset40k.fastq" , output_file = "ra3.BAM")
align.ra4 <- align(index = "ref_humaan", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2 = "SRR4785988_2_subset40k.fastq" , output_file = "ra4.BAM")
# Laad Rsamtools voor sorteren en indexeren (dowloaden indien nodig)
BiocManager::install('Rsamtools')
library(Rsamtools)

# Bestandsnamen van de monsters
samples <- c('eth1', 'eth2', 'eth3', 'con1', 'con2', 'con3')

# Voor elk monster: sorteer en indexeer de BAM-file
# Sorteer BAM-bestanden
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})
# Indexeer de gesorteerde BAM-file
lapply(samples, function(s) {indexBam(file = paste0(s, '.sorted.bam'))
})

allsamples <- c("nor1.BAM", "nor2.BAM", "nor3.BAM", "nor4.BAM","ra1.BAM" , "ra2.BAM" , "ra3.BAM" , "ra4.BAM")

count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

str(count_matrix)

countsgithub <- count_matrix$counts
head(countsgithub)

colnames(countsgithub) <- c("nor1" , "nor2" , "nor3" , "nor4" , "ra1" , "ra2", "ra3" , "ra4")
head(countsgithub)

write.csv(countsgithub , "githubreuma.csv")
 
# werkcollege 3

read.csv("githubreuma.csv")

echtedata <- read.table('/Users/sterremiedema/OneDrive - NHL Stenden/transcriptomicsgithub/count_matrix_RA.txt')

treatment <- c("control" , "control" , "control" , "control" , "reuma" , "reuma" ,"reuma" ,"reuma")
treatment_table <- data.frame(treatment)

head(treatment_table)

dds <- DESeqDataSetFromMatrix(countData = echtedata,
                              colData = treatment_table,
                              design = ~ treatment)
dds=DESeq(dds)
resultatenR=results(dds)
resultatenR

write.table(resultaten , file = 'Resultaten.csv' , row.names = TRUE , col.names = TRUE)

sum(resultatenR$padj < 0.05 & resultatenR$log2FoldChange > 1, na.rm = TRUE)
sum(resultatenR$padj < 0.05 & resultatenR$log2FoldChange < -1, na.rm = TRUE)

hoogste_fold_change <- resultatenR[order(resultatenR$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultatenR[order(resultatenR$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultatenR[order(resultatenR$padj, decreasing = FALSE), ]
hoogste_fold_change
laagste_fold_change
laagste_p_waarde

EnhancedVolcano(resultatenR,
                lab = rownames(resultatenR),
                x = 'log2FoldChange',
                y = 'padj')
dev.copy(png, 'VolcanoplotReuma.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

# je doet eerst GO analyse

BiocManager::install("goseq")
BiocManager::install("geneLenDataBase")
BiocManager::install("org.Hs.eg.db") # Hs inplaats van Dm
library("goseq")
library("geneLenDataBase")
library("org.Hs.eg.db")

ALL = rownames(resultatenR)
head(ALL)

res <- as.data.frame(resultatenR)
DEG = res %>%
  filter(padj <0.05)
DEG1 = rownames(DEG)

gene.vector=as.integer(ALL%in%DEG1)
names(gene.vector)=ALL
#lets explore this new vector a bit
head(gene.vector)
tail(gene.vector)

pwf=nullp(gene.vector,"hg19","geneSymbol")

GO.wall=goseq(pwf,"hg19","geneSymbol")

#How many enriched GO terms do we have
class(GO.wall)
head(GO.wall)
nrow(GO.wall)

enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]
#NOTE: They recommend using a more stringent multiple testing corrected p value here

#How many GO terms do we have now?
class(enriched.GO)
head(enriched.GO)
length(enriched.GO)

top10<- GO.wall %>% arrange(over_represented_pvalue) %>%
  slice(1:10)
# ggplot 1.0
ggplot(top10, aes(x = reorder(category, -over_represented_pvalue),
                  y = -log10(over_represented_pvalue))) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 meest verrijkte GO‑categorieën",
    x = "GO‑categorie",
    y = "-log10(p‑waarde)"
  ) +
  theme_minimal()
#ggplot 2.0
ggplot(top10, aes(x = -log10(over_represented_pvalue),
                  y = reorder(category, over_represented_pvalue))) +
  geom_point(size = 4, color = "darkred") +
  labs(
    title = "Top 10 GO‑categorieën (dotplot)",
    x = "-log10(p‑waarde)",
    y = "GO‑categorie"
  ) +
  theme_minimal()
#ggplot 3.0
ggplot(top10,
       aes(x = reorder(category, over_represented_pvalue),
           y = -log10(over_represented_pvalue),
           fill = -log10(over_represented_pvalue))) +
  
  geom_col(width = 0.8) +
  
  coord_flip() +
  
  scale_fill_gradient(
    low = "skyblue",
    high = "darkblue"
  ) +
  
  labs(
    title = "Top 10 meest verrijkte GO-categorieën",
    subtitle = "Over-representation analyse",
    x = "GO-categorie",
    y = expression(-log[10](pwaarde)),
    fill = expression(-log[10](pwaarde))
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

gene_conversion <- bitr(
  res_df$ENSEMBL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db )
res_df <- as.data.frame(resultaten)

res_df$SYMBOL <- rownames(res_df)

res_annotated <- merge(
  res_df,
  gene_conversion,
  by = "SYMBOL"
)
gene_list <- res_annotated$log2FoldChange

names(gene_list) <- res_annotated$ENTREZID 
pathview(
  gene.data = gene_list,
  pathway.id = "hsa05323",
  species = "hsa"
)

pathview(
  gene.data = resultatenR,
  pathway.id = "04010",  
  species = "hsa",          
  gene.idtype = "KEGG",     
  limit = list(gene = 5)    
)
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

# Convert DESeq2 results to dataframe
res_df <- as.data.frame(resultaten)

# Add Ensembl IDs from rownames
res_df$SYMBOL <- rownames(res_df)

# Convert ENSEMBL -> ENTREZID
gene_conversion <- bitr(
  res_df$ENSEMBL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Merge annotations
res_annotated <- merge(
  res_df,
  gene_conversion,
  by = "SYMBOL"
)


resultatenR[1] <- NULL
resultatenR[2:5] <- NULL

# Run pathview
pathview(
  gene.data = resultatenR,
  pathway.id = "hsa05323",
  species = "hsa",
  gene.idtype = "SYMBOL"
)


# -> extra boven resultaten # Create named vector for pathview
# gene_list <- res_annotated$log2FoldChange
# names(gene_list) <- res_annotated$ENTREZID