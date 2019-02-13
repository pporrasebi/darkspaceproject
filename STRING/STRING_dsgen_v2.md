---
output: 
  html_document: 
    keep_md: yes
---
STRING dataset generator, version 2
========================================================

#### Load STRING data



I download the latest version of the STRING data for human from STRING website: MITAB2.5. Current file was downloaded on Wed Feb 13 11:23:15 2019. 


```r
STRING_mitab_25_url = "http://string-db.org/download/psicquic-mitab_2.5.v10/9606.psicquic-mitab_2.5.v10.txt.gz"
STRING_mitab_25_file_gz = "./source_files/9606.psicquic-mitab_2.5.v10.txt.gz"
STRING_mitab_25_file = "./source_files/9606.psicquic-mitab_2.5.v10.txt"
if(!file.exists(STRING_mitab_25_file_gz)) download(STRING_mitab_25_url,STRING_mitab_25_file_gz)
gunzip(STRING_mitab_25_file_gz, remove = F, overwrite = T)
# Reading the table, cleaning
STRING_mitab = fread(STRING_mitab_25_file, colClasses = "character", header = F)
unlink(STRING_mitab_25_file)
STRING_proc = STRING_mitab[,.(ida = V1, idb = V2, method = V7, taxida = V10, taxidb = V11, interaction_type = V12, database = V13, score = V15)]
STRING_proc[, ida_ENSP := gsub("^string:9606\\.","", ida)]
STRING_proc[, ida_ENSP := gsub("\\|uniprotkb:.+$","", ida_ENSP)]
STRING_proc[, idb_ENSP := gsub("^string:9606\\.","", idb)]
STRING_proc[, idb_ENSP := gsub("\\|uniprotkb:.+$","", idb_ENSP)]
STRING_proc[, score := gsub("^score:","", score)]
```

#### STRING PMID-interactions text-mined pairs

I obtain the STRING text mining data from Felix Langer, who has made the full SQL installation of the STRING database to obtain the interaction-PMID links we lack. 


```r
pmid_pairs <- unique(fread(
        "./source_files/flanger_files/interactions_with_pmids.csv",
        sep = ",",
        header = T,
        colClasses = "character",
        stringsAsFactors = F
))
```

```
## Avoidable 10.177 seconds. This file is very unusual: it ends abruptly without a final newline, and also its size is a multiple of 4096 bytes. Please properly end the last row with a newline using for example 'echo >> file' to avoid this  time to copy.
```

```r
# Refs come from PMID and OMIM, need to clean them up and discard OMIM.

pmid_pairs_clean <- unique(pmid_pairs[grepl("PMID:",pmids),.(
        ida_ENSP = gsub(".*[.]","",protein1),
        taxID.a = gsub("[.].*","",protein1),
        idb_ENSP = gsub(".*[.]","",protein2),
        taxID.b = gsub("[.].*","",protein2),
        textmining,
        textmining_transferred,
        combined_score,
        pmids = gsub(", ",",",gsub("\'","",gsub("\\]","",gsub("\\[","",gsub("OMIM:[0-9]+\', ","",pmids)))))
)])
```


```r
string_proteins = STRING_proc[,unique(c(ida_ENSP, idb_ENSP))]
# map ENSEMBL_PROTEIN to UNIPROTKB
uniprot = UniProt.ws(taxId=9606)
# getting names of the keys and columns from Uniprot
# keytypes(uniprot)
# mapping IDs
string_proteins_up = data.table()
for(i in seq(1, length(string_proteins),500)){
ENSP2uniprot_temp <- as.data.table(select(uniprot, 
               keys = string_proteins[i:(i+500)], 
               columns = c("ENSEMBL_PROTEIN", "UNIPROTKB"),
               keytype = "ENSEMBL_PROTEIN"))
ENSP2uniprot = copy(ENSP2uniprot_temp[complete.cases(ENSP2uniprot_temp),])
string_proteins_up = rbind(string_proteins_up, ENSP2uniprot)
}

# Cleaning up mapping column

string_proteins_up <- unique(string_proteins_up[,.(
        ENSEMBL_PROTEIN,
        UNIPROTKB = gsub("-[0-9]+","",gsub(".* -> ","",UNIPROTKB))
)])
colnames(string_proteins_up)[2] = "id_up"

# Mapping STRING data obtained through their website
STRING_proc <- STRING_proc[, `:=` (ida_up = string_proteins_up$id_up[match(ida_ENSP,string_proteins_up$ENSEMBL_PROTEIN)], idb_up = string_proteins_up$id_up[match(idb_ENSP,string_proteins_up$ENSEMBL_PROTEIN)])]

STRING_proc = STRING_proc[complete.cases(STRING_proc)]

# Mapping STRING data from Felix pipeline
pmid_pairs_mapped <- pmid_pairs_clean[taxID.a=="9606" & taxID.b=="9606", `:=` (ida_up = string_proteins_up$id_up[match(ida_ENSP,string_proteins_up$ENSEMBL_PROTEIN)], idb_up = string_proteins_up$id_up[match(idb_ENSP,string_proteins_up$ENSEMBL_PROTEIN)])]

pmid_pairs_mapped = pmid_pairs_mapped[complete.cases(pmid_pairs_mapped)]
```


```r
STRING_proc[, sort(table(method), decreasing = T)]
```

```
## method
##                 psi-mi:"MI:0064"(interologs mapping) 
##                                               240998 
##                          psi-mi:"MI:0362"(inference) 
##                                               188579 
##             psi-mi:"MI:0087"(predictive text mining) 
##                                               166720 
## psi-mi:"MI:0045"(experimental interaction detection) 
##                                                89136 
##               psi-mi:"MI:0085"(phylogenetic profile) 
##                                                 2866 
##                      psi-mi:"MI:0036"(domain fusion) 
##                                                  114
```

```r
# generating interacting pairs
STRING_proc[, pair_id_clean := apply(data.table(ida_up,idb_up,stringsAsFactors = F), 1,
                                               function(a) { z = sort(a)
                                               paste0(z[1],"_",z[2]) })]

pmid_pairs_mapped[, pair_id_clean := apply(data.table(ida_up,idb_up,stringsAsFactors = F), 1,
                                               function(a) { z = sort(a)
                                               paste0(z[1],"_",z[2]) })]

text_mining_pmids <- unique(pmid_pairs_mapped[,.(
        pmids = paste(unique(pmids),collapse = ",")
),
by = c("pair_id_clean")])

# Bring PMIDs to STRING data 

STRING_proc_v2 <- STRING_proc[, `:=` (pmids = text_mining_pmids$pmids[match(pair_id_clean,text_mining_pmids$pair_id_clean)])]

STRING_textmining = STRING_proc_v2[method == "psi-mi:\"MI:0087\"(predictive text mining)",]
STRING_textmining[, sort(table(database), decreasing = T)]
```

```
## psi-mi:"MI:1014"(string) 
##                   166720
```

```r
STRING_textmining_pmids <- STRING_textmining[!is.na(pmids)]
STRING_textmining_pmid <- cSplit(
        STRING_textmining_pmids,
        splitCols = "pmids",
        sep = ",",
        direction = "long",
        drop = F
)

STRING_pathway_inference = STRING_proc_v2[method %in% c("psi-mi:\"MI:0362\"(inference)"),]
STRING_pathway_inference[, sort(table(database), decreasing = T)][1:20]
```

```
## database
##                                                 psi-mi:"MI:0467"(reactome)|psi-mi:"MI:0470"(kegg_pathways) 
##                                                                                                      62631 
##                                                                                 psi-mi:"MI:0467"(reactome) 
##                                                                                                      52614 
##                                                                            psi-mi:"MI:0470"(kegg_pathways) 
##                                                                                                      23874 
##                                                                            psi-mi:"MI:0448"(gene_ontology) 
##                                                                                                       6956 
##                                                                                      psi-mi:"MI:1107"(pid) 
##                                                                                                       4214 
##                                                      psi-mi:"MI:1107"(pid)|psi-mi:"MI:0470"(kegg_pathways) 
##                                                                                                       3002 
##                                                                                 psi-mi:"MI:1108"(biocarta) 
##                                                                                                       2205 
##                                                        psi-mi:"MI:0469"(intact)|psi-mi:"MI:0467"(reactome) 
##                                                                                                       1897 
##                                 psi-mi:"MI:0469"(intact)|psi-mi:"MI:0467"(reactome)|psi-mi:"MI:0471"(mint) 
##                                                                                                       1753 
##                                            psi-mi:"MI:0470"(kegg_pathways)|psi-mi:"MI:0448"(gene_ontology) 
##                                                                                                       1435 
##                                                 psi-mi:"MI:1108"(biocarta)|psi-mi:"MI:0470"(kegg_pathways) 
##                                                                                                       1211 
##                                                       psi-mi:"MI:0467"(reactome)|psi-mi:"MI:0463"(biogrid) 
##                                                                                                       1064 
##                 psi-mi:"MI:0467"(reactome)|psi-mi:"MI:0470"(kegg_pathways)|psi-mi:"MI:0448"(gene_ontology) 
##                                                                                                       1063 
##                        psi-mi:"MI:0469"(intact)|psi-mi:"MI:0467"(reactome)|psi-mi:"MI:0470"(kegg_pathways) 
##                                                                                                       1040 
##                           psi-mi:"MI:1107"(pid)|psi-mi:"MI:0467"(reactome)|psi-mi:"MI:0470"(kegg_pathways) 
##                                                                                                        795 
##          psi-mi:"MI:0469"(intact)|psi-mi:"MI:0738"(pride)|psi-mi:"MI:0467"(reactome)|psi-mi:"MI:1042"(pmc) 
##                                                                                                        704 
## psi-mi:"MI:0469"(intact)|psi-mi:"MI:0467"(reactome)|psi-mi:"MI:0470"(kegg_pathways)|psi-mi:"MI:0471"(mint) 
##                                                                                                        660 
##                                                   psi-mi:"MI:0469"(intact)|psi-mi:"MI:0448"(gene_ontology) 
##                                                                                                        628 
##                       psi-mi:"MI:0467"(reactome)|psi-mi:"MI:0470"(kegg_pathways)|psi-mi:"MI:0463"(biogrid) 
##                                                                                                        534 
##                                psi-mi:"MI:0469"(intact)|psi-mi:"MI:0738"(pride)|psi-mi:"MI:0467"(reactome) 
##                                                                                                        514
```

```r
STRING_phylo_predictions = STRING_proc_v2[method %in% c("psi-mi:\"MI:0064\"(interologs mapping)","psi-mi:\"MI:0085\"(phylogenetic profile)", "psi-mi:\"MI:0036\"(domain fusion)"),]
STRING_phylo_predictions[, sort(table(database), decreasing = T)]
```

```
## psi-mi:"MI:1014"(string) 
##                   243978
```

```r
# STRING_proc = STRING_proc[-grep("intact", database),]
# STRING_proc = STRING_proc[-grep("biogrid", database),]
# STRING_proc = STRING_proc[-grep("reactome", database),]
# STRING_proc = STRING_proc[-grep("hprd", database),]
# STRING_proc = STRING_proc[-grep("mint", database),]
# STRING_proc = STRING_proc[-grep("dip", database),]
# STRING_proc = STRING_proc[-grep("bind", database),]
```



```r
# saving STRING_textmining table with standard columns
STRING_textmining = STRING_textmining[,.(pair_id_clean, ida_clean = ida_up, idb_clean = idb_up, taxon = "9606", STRING_score = score, database, STRING_textmining = 1)]
fwrite(x = unique(STRING_textmining), 
       file = "./results/pairs_STRING_textmining.txt", sep = "\t")
# saving STRING_textmining_pmid table with standard columns
STRING_textmining_pmid = STRING_textmining_pmid[,.(pair_id_clean, ida_clean = ida_up, idb_clean = idb_up, taxon = "9606", pmid = pmids, STRING_score = score, database, STRING_textmining = 1)]
fwrite(x = unique(STRING_textmining_pmid), 
       file = "./results/pairs_STRING_textmining_pmid.txt", sep = "\t")
system("zip -rm ./results/pairs_STRING_textmining_pmid.txt.zip ./results/pairs_STRING_textmining_pmid.txt")
# saving STRING_pathway_inference table with standard columns
STRING_pathway_inference = STRING_pathway_inference[,.(pair_id_clean, ida_clean = ida_up, idb_clean = idb_up, taxon = "9606", STRING_score = score, database, STRING_pathway_inference = 1)]
fwrite(x = unique(STRING_pathway_inference), 
       file = "./results/pairs_STRING_pathway_inference.txt", sep = "\t")
# saving STRING_phylo_predictions table with standard columns
STRING_phylo_predictions = STRING_phylo_predictions[,.(pair_id_clean, ida_clean = ida_up, idb_clean = idb_up, taxon = "9606", STRING_score = score, database, STRING_phylo_predictions = 1)]
fwrite(x = unique(STRING_phylo_predictions), 
       file = "./results/pairs_STRING_phylo_predictions.txt", sep = "\t")
```

The total number of interacting pairs in the filtered by database STRING_textmining dataset: 166397

The total number of interacting pairs in the filtered by database STRING_pathway_inference dataset: 188526 

The total number of interacting pairs in the filtered by database STRING_phylo_predictions dataset: 240920 

#### Compare STRING "textmining", "pathway inference", "phylogeny and homology-based predictions" datasets

I calculate how many interactions overlap between different STRING datasets.


```r
area1 = STRING_textmining[STRING_textmining == 1,.N]
area2 = STRING_pathway_inference[STRING_pathway_inference == 1,.N]
area3 = STRING_phylo_predictions[STRING_phylo_predictions == 1,.N]
n12 = sum(!is.na(match(STRING_textmining[,unique(pair_id_clean)], STRING_pathway_inference[,unique(pair_id_clean)])))
n23 = sum(!is.na(match(STRING_pathway_inference[,unique(pair_id_clean)], STRING_phylo_predictions[,unique(pair_id_clean)])))
n13 = sum(!is.na(match(STRING_textmining[,unique(pair_id_clean)], STRING_phylo_predictions[,unique(pair_id_clean)])))
n123 = length(intersect(intersect(STRING_textmining$pair_id_clean,STRING_pathway_inference$pair_id_clean),STRING_phylo_predictions$pair_id_clean))

venn_Interactions_1 = draw.triple.venn(area1 = area1, 
                          area2 = area2, 
                          area3 = area3, 
                          n12 = n12, 
                          n23 = n23,
                          n13 = n13,
                          n123 = n123, 
                          category = c("STRING_textmining", 
                                       "STRING_pathway_inference", 
                                       "STRING_phylo_predictions"), 
                          lty = rep("blank", 3), 
                          fill = c("blue", "red", "green"), 
                          alpha = rep(0.5, 3), cat.pos = c(350, 25, 160), 
                          cat.dist = c(0.08,0.035, 0.08), 
                          cat.cex = c(0.9,0.9, 0.9), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 0.8
)
```

![](STRING_dsgen_v2_files/figure-html/STRING_comparison -1.png)<!-- -->


#### Compare STRING datasets and IMEx 

I calculate how many interactions in STRING datasets match to IMEx.  


```r
imex = fread("https://raw.githubusercontent.com/pporrasebi/darkspaceproject/master/IMEx/results/imex_full.txt", header = T, sep = "\t", colClasses = "character")
imex_human = imex[taxid_a == "9606" | taxid_b == "9606",]

area1 = STRING_textmining[STRING_textmining == 1,.N]
area2 = STRING_pathway_inference[STRING_pathway_inference == 1,.N]
area3 = STRING_phylo_predictions[STRING_phylo_predictions == 1,.N]
area4 = imex_human[,length(unique(pair_id))]
n12 = sum(!is.na(match(STRING_textmining[,unique(pair_id_clean)], STRING_pathway_inference[,unique(pair_id_clean)])))
n13 = sum(!is.na(match(STRING_textmining[,unique(pair_id_clean)], STRING_phylo_predictions[,unique(pair_id_clean)])))
n14 = sum(!is.na(match(STRING_textmining[,unique(pair_id_clean)], imex_human[,unique(pair_id)])))
n23 = sum(!is.na(match(STRING_pathway_inference[,unique(pair_id_clean)], STRING_phylo_predictions[,unique(pair_id_clean)])))
n24 = sum(!is.na(match(STRING_pathway_inference[,unique(pair_id_clean)], imex_human[,unique(pair_id)])))
n34 = sum(!is.na(match(STRING_phylo_predictions[,unique(pair_id_clean)], imex_human[,unique(pair_id)])))
n123 = length(intersect(intersect(unique(STRING_textmining$pair_id_clean),unique(STRING_pathway_inference$pair_id_clean)),unique(STRING_phylo_predictions$pair_id_clean)))
n124 = length(intersect(intersect(unique(STRING_textmining$pair_id_clean),unique(STRING_pathway_inference$pair_id_clean)),unique(imex_human$pair_id)))
n134 = length(intersect(intersect(unique(STRING_textmining$pair_id_clean),unique(STRING_phylo_predictions$pair_id_clean)),unique(imex_human$pair_id)))
n234 = length(intersect(intersect(unique(STRING_pathway_inference$pair_id_clean),unique(STRING_phylo_predictions$pair_id_clean)),unique(imex_human$pair_id)))
n1234 = length(intersect(intersect(intersect(unique(STRING_textmining$pair_id_clean),unique(STRING_pathway_inference$pair_id_clean)),unique(STRING_phylo_predictions$pair_id_clean)), unique(imex_human$pair_id)))


venn_Interactions_1 = draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24,
    n34, n123, n124, n134, n234, n1234, 
                          category = c("STRING_textmining", 
                                       "STRING_pathway_inference", 
                                       "STRING_phylo_predictions",
                                       "IMEx"), 
                          lty = rep("blank", 4), 
                          fill = c("blue", "red", "green", "grey"), 
                          alpha = rep(0.5, 4), cat.pos = c(350, 25, 160,45), 
                          cat.dist = c(0.08,0.035, -0.08, 0.08), 
                          cat.cex = c(0.9,0.9, 0.9, 0.9), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 0.8
)
```

![](STRING_dsgen_v2_files/figure-html/STRING_imex_comparison -1.png)<!-- -->
