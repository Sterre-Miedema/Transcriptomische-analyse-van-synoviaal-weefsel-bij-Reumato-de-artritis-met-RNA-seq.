# Transcriptomische analyse van synoviaal weefsel bij *Reumatoïde artritis* met RNA-seq.





## Inleiding

*Reumatoïde artritis* (RA) is een chronische, systemische auto-immuunziekte die wereldwijd ongeveer 0,3/0,5-1% van de bevolking treft [(1,2)](./Bronnen). Het treft vrouwen vaker van mannen [(1,2,3)](./Bronnen). Deze ontstekingsziekte ontstaat doordat het afweersysteem de eigen gewrichten aan, waardoor zowel grote als kleine gewrichten ontstoken raken [(3)](./Bronnen). De exacte oorzaak van RA is nog niet volledig bekend maar de ziekte ontstaat vermoedelijk door een complexe interactie tussen genetische aanleg, verandering in genexpressie, omgevingsfactoren en een ontregeld immuunsysteem [(1)](./Bronnen). Omdat blijvende schade al vroeg kan ontstaan tijdens de ziekte door aanhoudende ontsteking [(1)](./Bronnen), is een snelle diagnose en behandeling essentieel. Hoewel verschillende geneesmiddelen de ziekteactiviteit kunnen remmen, bestaat er momenteel geen genezende therapie. 
Door de grote variatie in symptomen en behandelrespons wordt aangenomen dat RA uit meerdere subtypen bestaat [(1)](./Bronnen). Transcriptomics kan hierin waardevolle mogelijkheden bieden. Met methoden zoals RNA-sequencing kan de expressie van duizenden genen tegelijk worden onderzocht waardoor er een overzicht ontstaat in moleculaire mechanismen en signaalroutes [(4)](./Bronnen). Aangezien RA een ontstekingsreactie geeft in het weefsel, kan transcriptomics een mogelijke bijdrage leveren aan een betere diagnose. Het doel van deze studie is om met behulp van transcriptomics te onderzoeken welke genen en biologische routes verschillend tot expressie komen in synoviaal weefsel van patiënten met RA ten opzichte van gezonde controles.

Deelvragen:

- Welke genen komen differentieel tot expressie in synoviaal weefsel van patiënten met *Reumatoïde artritis* ten opzichte van de gezonde controle?

- Welke biologische processen en signaalroutes tonen een verschil in synoviaal weefsel van patiënten met *Reumatoïde artritis* ten opzichte van de gezonde controle door middel van de GO- en KEGG-analyse?

## Materiaal en methode
Voor deze transcriptomics analyse werden RNA-sequenties gebruikt en Fastq bestanden gedownload afkomstig uit de studie van platzer et al. (2019) [(1)](./Bronnen). In de  studie werden synoviumbiopten verzameld van 4 gezonde individuen, en 4 patiënten met vastgestelde Reumatoïde artritis (RA) uit het gewrichtsslijmvlies. 

Alle analyses werden uitgevoerd in R (versie > 4.0) op macOS. De volgende Bioconductor- packages werden gebruikt: Rsubread(2.24.0) [(5)](./Bronnen), Rsamtools (2.26.0) [(6)](./Bronnen), DESeg2 (1.50.2) [(7)](./Bronnen), EnhancedVolcano [(8)](./Bronnen), goseq (1.62.0) [(9)](./Bronnen), geneLenDataBase (1.28.2) [(10)](./Bronnen), org.Hs.eg.db (3.22.0) [(11)](./Bronnen), clusterProfiler (4.18.4) [(11)](./Bronnen), pathview (1.50.0) [(11)](./Bronnen), dplyr (1.2.1) [(12)](./Bronnen) en ggplot2 (4.0.3) [(13)](./Bronnen). Alle packages werden geïnstalleerd via BiocManager. 

Het humane referentiegenoom [(GCF_000001405.40)](Bronnen/GCF_000001405.40.pdf) en de GTF-annotatie werden gedownload vanaf NCBI. 
Vervolgens weden de BAM-bestanden gesorteerd en geïndexeerd. Met featureCounts werd een count matrix opgesteld. Deze count matrix werd gebruikt voor de differentiële expressieanalyse met DESeq2 (pad < 0,05;Log2FoldChange > 1). De RNA-seq analyse identificeerde 29407 genen in synovia weefsel

De resultaten van de expressieanalyse werden gevisualiseerd in een volcanoplot met EnhancedVolcano. Vervolgens werd een GO-analyse uitgevoerd met goseq om de biologische processen te identificeren. Hieruit is een KEGG-pathway (hsa04660; T cell receptor signaling pathway) gekozen. Met bitr() zijn de gen-ID's omgezet naar Entrez-ID's.



## Resultaten

De RNA-seq analyse identificeerde 29407 genen. Met DESeq2 werden tot expressie gekomen genen bepaald (pad < 0,05; log2FoldChange > 1). Vooral immuunglobulinegenen (IGHV3-53 [(14)](./Bronnen), IGHV1-69 [(15)](./Bronnen) en IGHV4-31 [(16)](./Bronnen)) en ontstekingsgerelateerde genen (CXCR1 [(17)](./Bronnen) en PTGFR [(18)](./Bronnen)) waren verhoogd tot expressie gebracht.

**Verhoogde expressie van immuunglobuline- en ontsteking gerelateerde genen in Reumatoïde artritis vergeleken met de controlegroep.**

De RNA-seq analyse van synoviaal weefsel lieten duidelijke verschillen zien in genexpressie. In totaal werden er 29407 genen getest op differentiële expressie en deze resultaten zijn in een volcano plot (figuur 1) uitgezet. Op basis van de drempels padj < 0.05 en log2FC > 1 is er een significante toe- of afname zichtbaar in verschillende immuunglobuline-gerelateerde genen(IGHV3-53 [(14)](./Bronnen), IGHV1-69 [(15)](./Bronnen) en IGHV4-31 [(16)](./Bronnen)) en ontsteking gerelateerde genen (CXCR1 [(17)](./Bronnen) en PTGFR [(18)](./Bronnen)). [Figuur 1](Resultaten/VolcanoplotReuma.png) toont 3 groepen genen: Grijs: niet-significante genen, groen: genen met grote fold-change maar niet significant en rood: genen die een grote fold change hebben en significant zijn.


**Verhoogde activatie van immuunprocessen in RA.**

De GO-analyse laat een verijking zien van biologische processen betrokken bij immuungerelateerde processen [(figuur 2)](Resultaten/GO-analysereuma.png). De meest significante GO-categorie is het immunoglobulin complex, gevolgd door adaptive immune respons. Deze resultaten wijzen op een sterke activatie van immuunprocessen in het synoviale weefsel van patiënten met RA. Voor de KEGG-analyse wordt het immuunsysteem bekeken.

**Verhoogde T-celactivatie en MAPK-signaaltransductie in RA.**

De KEGG-analyse [(figuur 3)](Resultaten/hsa04660.pathview.png) van de T-celreceptor signaalroute (hsa04660) laat zien dat T-celactivatie verhoogd tot expressie komt. De pathway toont dat zowel co-stimulatoire moleculen zoals CD28, CTLA4 en CD86 en de MAPK-route verhoogde expressie tonen. MAPK is een belangrijke route voor signaaltransductieroutes zoals celprofilatie, celdiffenrentiatie en celdood [(19)](./Bronnen) en dus belangrijke factor in ontstekingsprocessen [(20)](./Bronnen).


## Conclusie
Op basis van de uitgevoerde RNA-seq analyse van synoviaal weefsel uit *reumatoïde artritis* (RA) en gezonde controles kan geconcludeerd worden dat er een verhoogde ontstekingsactiviteit en immuun activatie is. Vooral immuunglobulinegenen en T-cel co-stimulatie moleculen (CD28, CTLA4 en CD86).De analyse van de KEGG-pathway laat zien dat er een belangrijke signaalroute zoals MAPK  en T-celreceptorroute sterk geactiveerd zijn. Dit draagt bij aan celprofilatie, celdifferentiatie en celdood [(19)](./Bronnen) maar ook een grote rol spelen in ontstekingsprocessen [(20)](./Bronnen). Dit komt overeen met de al bekende kenmerken van RA. Ook sluiten de resultaten  goed aan bij de bevindingen uit de studie van Platzer et al. (2019). Hieruit is de gebruikte dataset afkomstig. De combinatie van de DE-analyse, GO-verrijking en de KEGG-pathway visualisatie geeft een goed beeld van RA door ontstekingen en auto-immuun activatie. Het doel is dus behaald.

