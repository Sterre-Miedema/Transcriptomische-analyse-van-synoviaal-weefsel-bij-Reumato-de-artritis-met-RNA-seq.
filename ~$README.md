## Inleiding

Reumatoïde artritis (RA) is een chronische, systemische auto-immuunziekte die wereldwijd ongeveer 0,3/0,5-1% van de bevolking treft. (1, 2) Het treft vrouwen vaker van mannen (1,2,3). Deze ontstekingsziekte ontstaat doordat het afweersysteem de eigen gewrichten aan, waardoor zowel grote als kleine gewrichten ontstoken raken. (3) De exacte oorzaak van RA is nog niet volledig bekend maar de ziekte ontstaat vermoedelijk door een complexe interactie tussen genetische aanleg, verandering in genexpressie, omgevingsfactoren en een ontregeld immuunsysteem (1). Omdat blijvende schade al vroeg kan ontstaan tijdens de ziekte door aanhoudende ontsteking (1), is een snelle diagnose en behandeling essentieel. Hoewel verschillende geneesmiddelen de ziekteactiviteit kunnen remmen, bestaat er momenteel geen genezende therapie. 
Door de grote variatie in symptomen en behandelrespons wordt aangenomen dat RA uit meerdere subtypen bestaat (1). Transcriptomics kan hierin waardevolle mogelijkheden bieden. Met methoden zoals RNA-sequencing kan de expressie van duizenden genen tegelijk worden onderzocht waardoor er een overzicht ontstaat in moleculaire mechanismen en signaalroutes (4). Aangezien RA een ontstekingsreactie geeft in het weefsel, kan transcriptomics een mogelijke bijdrage leveren aan een betere diagnose.


## Materiaal en methode
Voor deze transcriptomics analyse werden RNA-sequenties gebruikt afkomstig uit de studie van platzer et al. (2019). In deze studie werden synoviumbiopten verzameld van 4 gezonde individuen, en 4 patiënten met vastgestelde Reumatoïde artritis (RA) uit het gewrichtsslijmvlies. Alle deelnemers waren tussen de 15 en 66 jaar en vrouwelijk. De gebruikte samples zijn:
SRR4785819
31VrouwNormaalSRR4785820
15VrouwNormaalSRR4785828
31VrouwNormaalSRR4785831
42VrouwNormaalSRR4785979
54VrouwRASRR4785980
66VrouwRASRR4785986
60VrouwRASRR4785988
59VrouwRA
De ruwe RNA-seq FASTQ bestanden werden vanuit de Sequence Read Archive (SRA) gedownload. De verdere analyse werd uitgevoerd in R. 

Alle analyses werden uitgevoerd in R (versie ? 4.0) op macOS. De volgende Bioconductor- packages werden gebruikt: Rsubread, Rsamtools, DESeg2, EnhancedVolcano, goseq, geneLenDataBase, org.Hs.eg.db, clusterProfiler, pathview en ggplot2. Alle packages werden geïnstalleerd via BiocManager. 

Het humane referentiegenoom (GCF_000001405.40) en de GTF-annotatie werden gedownload vanaf NCBI. Met de functie buildindex() uit Rsubread werd een genoomindex aangemaakt: memory = 4000, indexSlit = TRUE en basename = “ref_humaan”.

Voor elk sample werden de paired-end FASTQ-bestanden uitgelijnd tegen het humane referentiegenoom met de functie align(). Dit resulteerde in 8 BAM-bestanden (nor1,2,3,4 en ra1,2,3,4).

Hierna werden de 8 BAM-bestanden gesorteerd via sortBam() en geïndexeerd met indexBam()

Hierna werden de genexpressieniveau’s bepaald door featureCounts:
isPairedEnd = TRUE, isGTFAnnotationFile = TRUE, GTF.attrType = “gene_id”, use MetaFeatures = TRUE. De count matrix werd hierna opgeslagen als csv bestand. 

De counts werden ingelezen en gecombineerd in een metadata-tabel waarin de condities waren vastgelegd: control (NOR) en reuma (RA). Significante genen werden gedefinieerd als padj < 0.05 en log2FoldChange > 1.

Voor de visualisatie van differentiële expressie werd een volcanoplot gegenereerd met EnhancedVolcano waarin log2FoldChange (x-as) en padj (y-as) werden weergegeven en opgeslagen als PNG. 

Omdat RNA-seq gevoelig is voor genlengte-bias (ff ander woord voor zoeken) werd goseq gebruikt. Hierin werd een binaire vector gemaakt met onderscheid tussen wel of niet differentieel tot expressie komt.

Voor een KEGG-analyse moesten de genen omgezet worden naar Entrez-ID’s. met bitr() werden SYMBOL -> ENTREZID conversies uitgevoerd (wat betekend conversies?). Met pathview werden genexpressiewijzigingen geprojecteerd op KEGG-pathways (hsa05323 – RA, hsa04010 – MAPK signaling pathway).


