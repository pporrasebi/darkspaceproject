---
output: 
  html_document: 
    keep_md: yes
---
GO IPI, Gene Ontology - Inferred from Physical Interaction - dataset generator
========================================================

#### Load Gene Ontology IPI dataset - EBI_GOA_nonIntAct


```
## Warning: replacing previous import 'IRanges::desc' by 'plyr::desc' when loading
## 'PSICQUIC'
```


I download the latest version of Gene Ontology IPI dataset - EBI_GOA_nonIntAct - using PSICQUIC service. Current file was downloaded on Fri Jan 17 11:32:04 2020. 



```r
## Load PSICQUIC functionality
    psicquic <- PSICQUIC()
    providers <- providers(psicquic)
# query EBI-GOA-nonIntAct database		
database = "EBI-GOA-nonIntAct"
SPECIES_ID = "9606"
MITAB = "tab25"
EBI_GOA_nonIntAct = data.table()
if(database %in% providers){
        ## Query for the number of interactions
        PSICQUIC_query = paste("taxidA:",SPECIES_ID," AND taxidB:",SPECIES_ID, sep = "")
        PSICQUIC_query1 = paste0(PSICQUIC_query, "?format=count") 
        N_interactions <- unlist(rawQuery(psicquic, database, PSICQUIC_query1))
        
        if(N_interactions > 0){
          N_start = 1
          N_nrows = 2500
          for(n_starts in seq(from = N_start, to = N_interactions, by = N_nrows)){
            PSICQUIC_query2 = paste0(PSICQUIC_query, "?format=",MITAB,"&firstResult=", n_starts,"&maxResults=", N_nrows) 
            SPECIES_ID_interactome_d <- as.data.table(rawQuery(psicquic, database, PSICQUIC_query2))
            EBI_GOA_nonIntAct <- rbind(EBI_GOA_nonIntAct, SPECIES_ID_interactome_d)
          }
        }
}
fwrite(x = unique(EBI_GOA_nonIntAct), 
       file = "./source_files/EBI_GOA_nonIntAct_mitab25.txt", sep = "\t")

system("perl ./scripts/MITAB25extractor_v12.pl ./source_files/EBI_GOA_nonIntAct_mitab25.txt ./processed_files/EBI_GOA_nonIntAct_pairs.txt")
```


Saving a table of interacting pairs, publication IDs and EBI_GOA_nonIntAct tag.



```r
EBI_GOA_nonIntAct = fread("./processed_files/EBI_GOA_nonIntAct_pairs.txt", header = T, sep = "\t", colClasses = "character")
fwrite(x = unique(EBI_GOA_nonIntAct[, .(pair_id_clean, pubid, EBI_GOA_nonIntAct = rep(1, .N))]), 
       file = "./results/pairs_pmids_EBI_GOA_nonIntAct.txt", sep = "\t")
N_EBI_GOA_nonIntAct = length(EBI_GOA_nonIntAct[,unique(pair_id_clean)])
```


The EBI_GOA_nonIntAct dataset contains 10456 human interacting pairs. 


Creating a list of PMIDs that have been curated into GO annotations of human proteins (EBI_GOA_nonIntAct dataset)



```r
EBI_GOA_nonIntAct_pmids <- data.frame(unique(EBI_GOA_nonIntAct$pubid))

write.table(EBI_GOA_nonIntAct_pmids, "./results/biogrid_pmids.txt", quote=F, sep ="\t", row.names = F, col.names = T)
```


5827 publications (human) are curated into EBI_GOA_nonIntAct dataset 

#### Compare human EBI_GOA_nonIntAct interactions and publications to IMEx 

I calculate how many interactions in EBI_GOA_nonIntAct dataset match to IMEx.


```r
imex = fread("https://raw.githubusercontent.com/pporrasebi/darkspaceproject/master/IMEx/results/imex_full.txt", header = T, sep = "\t", colClasses = "character")
imex_human = imex[taxid_a == "9606" & taxid_b == "9606",]
N_imex = length(imex_human[, unique(pair_id_clean)])
N_EBI_GOA_nonIntAct = length(EBI_GOA_nonIntAct[,unique(pair_id_clean)])
N_overlap = sum(!is.na(match(EBI_GOA_nonIntAct[,unique(pair_id_clean)], imex_human[, unique(pair_id_clean)])))

venn.d = draw.pairwise.venn(area1 = N_imex, area2 = N_EBI_GOA_nonIntAct, cross.area = N_overlap, category = c("IMEx", "EBI_GOA_nonIntAct"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(0, 135), 
                          cat.dist = rep(0.035, 2), 
                          cat.cex = c(1,1), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          direct.area = TRUE,
                          cex = 1)
```

![](GO_IPI_dsgen_files/figure-html/biogrid_vs_imex-1.png)<!-- -->

I calculate how many publications in EBI_GOA_nonIntAct dataset match to IMEx.


```r
N_pubid_imex = length(imex_human[, unique(pubid)])
N_pubid_EBI_GOA_nonIntAct = length(EBI_GOA_nonIntAct[,unique(pubid)])
N_pubid_overlap = sum(!is.na(match(EBI_GOA_nonIntAct[,unique(pubid)], imex_human[, unique(pubid)])))

venn.d = draw.pairwise.venn(area1 = N_pubid_imex, area2 = N_pubid_EBI_GOA_nonIntAct, cross.area = N_pubid_overlap, category = c("IMEx", "EBI_GOA_nonIntAct"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(0, 135), 
                          cat.dist = rep(0.035, 2), 
                          cat.cex = c(1,1), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          direct.area = TRUE,
                          cex = 1)
```

![](GO_IPI_dsgen_files/figure-html/biogrid_vs_imex_pub-1.png)<!-- -->

I calculate how many interactions published in specific articles (the same interaction can come from different publications) in EBI_GOA_nonIntAct dataset match to IMEx.


```r
N_pub_int_imex = length(imex_human[, unique(paste0(pubid,"_",pair_id_clean))])
N_pub_int_EBI_GOA_nonIntAct = length(EBI_GOA_nonIntAct[,unique(paste0(pubid,"_",pair_id_clean))])
N_pub_int_overlap = sum(!is.na(match(EBI_GOA_nonIntAct[,unique(paste0(pubid,"_",pair_id_clean))], imex_human[, unique(paste0(pubid,"_",pair_id_clean))])))

venn.d = draw.pairwise.venn(area1 = N_pub_int_imex, area2 = N_pub_int_EBI_GOA_nonIntAct, cross.area = N_pub_int_overlap, category = c("IMEx", "EBI_GOA_nonIntAct"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(0, 135), 
                          cat.dist = rep(0.035, 2), 
                          cat.cex = c(1,1), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          direct.area = TRUE,
                          cex = 1)
```

![](GO_IPI_dsgen_files/figure-html/biogrid_vs_imex_pub_inter-1.png)<!-- -->
