---
output: 
  html_document: 
    keep_md: yes
---
BioGRID dataset generator
========================================================

#### Load BioGRID data



I download the latest version of the BioGRID data from mentha website. I have tried using PSICQUIC, however, downloading whole mentha using PSICQUIC takes too much time. Now, however, the link should be updated manually to get the latest mentha release. Current file was downloaded on Fri Jan 17 11:19:28 2020. 


```r
if(!file.exists("./source_files/2020-01-14_MITAB-2.5.zip")){
        mentha_url = "http://mentha.uniroma2.it/dumps/2020-01-14_MITAB-2.5.zip"
        # httr::GET(http://mentha.uniroma2.it/doDownload.php?file=2017-02-27_MITAB-2.5.zip) - this is how I got the direct link
        download(mentha_url, destfile = "./source_files/2020-01-14_MITAB-2.5.zip")
}
unzip("./source_files/2020-01-14_MITAB-2.5.zip", exdir = "./source_files/")
```

Filtering BioGRID only interactions and cleaning up the BioGRID file.


```r
mentha = fread("./source_files/2020-01-14_MITAB-2.5", 
               header = F, 
               sep = "\t",
               colClasses = "character")
biogrid_from_mentha = mentha[V13 == "psi-mi:\"MI:0463\"(biogrid)",]
biogrid_from_mentha[,table(V12)]
```

```
## V12
##       psi-mi:"MI:0403"(colocalization)   psi-mi:"MI:0407"(direct interaction) 
##                                  32713                                 167508 
##          psi-mi:"MI:0914"(association) psi-mi:"MI:0915"(physical association) 
##                                   9412                                 349991
```

```r
fwrite(x = biogrid_from_mentha, file = "./processed_files/biogrid_mitab25.txt",sep = "\t")
system("perl ./scripts/MITAB25extractor_v12.pl ./processed_files/biogrid_mitab25.txt ./processed_files/biogrid_pairs.txt")
unlink(c("./source_files/2020-01-14_MITAB-2.5","./processed_files/biogrid_mitab25.txt"))
```

Saving a table of interacting pairs, publication IDs and BioGRID tag. Selecting human protein-human protein only interactions.


```r
biogrid_from_mentha = fread("./processed_files/biogrid_pairs.txt", header = T, sep = "\t", colClasses = "character")
biogrid_from_mentha_human = biogrid_from_mentha[taxid_a == "9606" & taxid_b == "9606",]
fwrite(x = unique(biogrid_from_mentha_human[, .(pair_id_clean, pubid, biogrid = rep(1, .N))]), 
       file = "./results/pairs_pmids_biogrid.txt", sep = "\t")
N_biogrid = length(biogrid_from_mentha_human[,unique(pair_id_clean)])
```

The BioGRID dataset contains 247094 human interacting pairs. 

Files in the /processed_files/ contain interactions for all species.
Files in the /results/ contain interactions for human only.

Creating a list of PMIDs that have been curated into human-human, human-other_species, human-compound interaction


```r
biogrid_pmids <- data.frame(unique(biogrid_from_mentha_human$pubid))

write.table(biogrid_pmids, "./results/biogrid_pmids.txt", quote=F, sep ="\t", row.names = F, col.names = T)
```

25213 publications (human) are curated into BioGRID database. 

#### Compare human BioGRID interactions and publications to IMEx 

I calculate how many interactions in BioGRID match to IMEx.


```r
imex = fread("https://raw.githubusercontent.com/pporrasebi/darkspaceproject/master/IMEx/results/imex_full.txt", header = T, sep = "\t", colClasses = "character")

imex_human = imex[taxid_a == "9606" & taxid_b == "9606",]
N_imex = length(imex_human[, unique(pair_id_clean)])
N_biogrid = length(biogrid_from_mentha_human[,unique(pair_id_clean)])
N_overlap = sum(!is.na(match(biogrid_from_mentha_human[,unique(pair_id_clean)], imex_human[, unique(pair_id_clean)])))

venn.d = draw.pairwise.venn(area1 = N_imex, area2 = N_biogrid, cross.area = N_overlap, category = c("IMEx", "BioGRID"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(0, 135), 
                          cat.dist = rep(0.035, 2), 
                          cat.cex = c(1,1), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          direct.area = TRUE,
                          cex = 1)
```

![](BioGRID_dsgen_files/figure-html/biogrid_vs_imex-1.png)<!-- -->

I calculate how many publications in BioGRID match to IMEx.


```r
N_pubid_imex = length(imex_human[, unique(pubid)])
N_pubid_biogrid = length(biogrid_from_mentha_human[,unique(pubid)])
N_pubid_overlap = sum(!is.na(match(biogrid_from_mentha_human[,unique(pubid)], imex_human[, unique(pubid)])))

venn.d = draw.pairwise.venn(area1 = N_pubid_imex, area2 = N_pubid_biogrid, cross.area = N_pubid_overlap, category = c("IMEx", "BioGRID"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(0, 135), 
                          cat.dist = rep(0.035, 2), 
                          cat.cex = c(1,1), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          direct.area = TRUE,
                          cex = 1)
```

![](BioGRID_dsgen_files/figure-html/biogrid_vs_imex_pub-1.png)<!-- -->

I calculate how many interactions published in specific articles (the same interaction can come from different publications) in BioGRID match to IMEx.


```r
N_pub_int_imex = length(imex_human[, unique(paste0(pubid,"_",pair_id_clean))])
N_pub_int_biogrid = length(biogrid_from_mentha_human[,unique(paste0(pubid,"_",pair_id_clean))])
N_pub_int_overlap = sum(!is.na(match(biogrid_from_mentha_human[,unique(paste0(pubid,"_",pair_id_clean))], imex_human[, unique(paste0(pubid,"_",pair_id_clean))])))

venn.d = draw.pairwise.venn(area1 = N_pub_int_imex, area2 = N_pub_int_biogrid, cross.area = N_pub_int_overlap, category = c("IMEx", "BioGRID"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(0, 135), 
                          cat.dist = rep(0.035, 2), 
                          cat.cex = c(1,1), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          direct.area = TRUE,
                          cex = 1)
```

![](BioGRID_dsgen_files/figure-html/biogrid_vs_imex_pub_inter-1.png)<!-- -->
