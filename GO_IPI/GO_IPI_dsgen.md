GO IPI, Gene Ontology - Inferred from Physical Interaction - dataset generator
========================================================

#### Load BioGRID data



I download the latest version of the GO annotations (UniProt KnowledgeBase (UniProtKB), IntAct protein complexes, and RNAcentral identifiers - multispecies) data from geneontology.org. Current file was downloaded on Wed Mar  1 17:20:39 2017. 


```r
if(!file.exists("./source_files/goa_uniprot_all_noiea.gaf.gz")){
        GO_url = "http://geneontology.org/gene-associations/goa_uniprot_all_noiea.gaf.gz"
        download(GO_url, destfile = "./source_files/goa_uniprot_all_noiea.gaf.gz")
}
gunzip(filename = "./source_files/goa_uniprot_all_noiea.gaf.gz", destname = "./source_files/goa_uniprot_all_noiea.gaf", remove = F)
readLines("./source_files/goa_uniprot_all_noiea.gaf", n = 10)
```

```
##  [1] "!gaf-version: 2.0"                                                                                                         
##  [2] "!Project_name: UniProt GO Annotation program (UniProt-GOA)"                                                                
##  [3] "!URL: http://www.ebi.ac.uk/GOA"                                                                                            
##  [4] "!Contact Email: goa@ebi.ac.uk"                                                                                             
##  [5] "!"                                                                                                                         
##  [6] "!This file contains all GO annotations and gene product information for proteins in the UniProt KnowledgeBase (UniProtKB),"
##  [7] "!IntAct protein complexes, and RNAcentral identifiers."                                                                    
##  [8] "!"                                                                                                                         
##  [9] "!Generated: 2017-02-13 16:27"                                                                                              
## [10] "!GO-version: http://purl.obolibrary.org/obo/go/releases/2017-02-11/go.owl"
```

Reading and filtering GO by evidence code "IPI" - Inferred from Physical Interaction (http://www.geneontology.org/page/ipi-inferred-physical-interaction).


```r
GO = fread("./source_files/goa_uniprot_all_noiea.gaf", skip = 11, header = F, sep = "\t", colClasses = "character")
```

```
## 
Read 45.6% of 1338386 rows
Read 85.9% of 1338386 rows
Read 1338386 rows and 17 (of 17) columns from 0.215 GB file in 00:00:04
```

```r
GO_IPI = unique(GO[V7 == "IPI",])
unlink(c("./source_files/goa_uniprot_all_noiea.gaf"))
```

Below you can see how many of GO IPI annotations originate from each source database.


```r
GO_IPI[,table(V15)]
```

```
## V15
##           AgBase            AspGD          BHF-UCL            CACAO 
##              777               16               74                1 
##              CGD            DFLAT          FlyBase           GeneDB 
##               22                4                7              126 
##       GO_Central               GR             HGNC           IntAct 
##               15                3                3             4516 
##              MGI          MTBBASE ParkinsonsUK-UCL        PseudoCAP 
##               26              504               14                7 
##          UniProt         WormBase 
##             1385                3
```

Generating a list of PMIDs that have been curated into GO IPI annotations


```r
GO_IPI[,V6 := gsub("^PMID:", "", V6)]
```

```
##               V1                  V2                  V3 V4         V5
##    1:  UniProtKB          A0A024A2C9                 lph    GO:0005515
##    2:  UniProtKB          A0A024A2C9                 lph    GO:0005515
##    3:  UniProtKB          A0A024A2C9                 lph    GO:0005515
##    4:  UniProtKB          A0A024A2C9                 lph    GO:0005515
##    5:  UniProtKB          A0A024A2C9                 lph    GO:0005515
##   ---                                                                 
## 7499:  UniProtKB              X2KN52              X2KN52    GO:0097677
## 7500:  UniProtKB              X2KN52              X2KN52    GO:0097677
## 7501:  UniProtKB              X2KN52              X2KN52    GO:0097677
## 7502:  UniProtKB              X2KN52              X2KN52    GO:1990782
## 7503: RNAcentral URS00006F9B91_11676 URS00006F9B91_11676    GO:0008134
##             V6  V7
##    1: 24835392 IPI
##    2: 26538390 IPI
##    3: 26538390 IPI
##    4: 26538390 IPI
##    5: 26538390 IPI
##   ---             
## 7499: 25379383 IPI
## 7500: 25379383 IPI
## 7501: 25379383 IPI
## 7502: 25379383 IPI
## 7503: 25116364 IPI
##                                                                                         V8
##    1:                                                                     UniProtKB:P08603
##    2:                                                                     UniProtKB:P04004
##    3:                                                                     UniProtKB:P04004
##    4:                                                                     UniProtKB:P08603
##    5:                                                                     UniProtKB:P08603
##   ---                                                                                     
## 7499:                                                                     UniProtKB:A7LH49
## 7500:                                                                     UniProtKB:B5X3E8
## 7501:                                                                     UniProtKB:C7EY83
## 7502: UniProtKB:S5R2L7|UniProtKB:S5R5H7|UniProtKB:S5RA51|UniProtKB:S5RNB2|UniProtKB:S5RRC4
## 7503:                                                                     UniProtKB:P04608
##       V9
##    1:  F
##    2:  F
##    3:  F
##    4:  F
##    5:  F
##   ---   
## 7499:  F
## 7500:  F
## 7501:  F
## 7502:  F
## 7503:  F
##                                                                          V10
##    1:                                                 Lipoprotein binding FH
##    2:                                                 Lipoprotein binding FH
##    3:                                                 Lipoprotein binding FH
##    4:                                                 Lipoprotein binding FH
##    5:                                                 Lipoprotein binding FH
##   ---                                                                       
## 7499:                       Signal transducer and activator of transcription
## 7500:                       Signal transducer and activator of transcription
## 7501:                       Signal transducer and activator of transcription
## 7502:                       Signal transducer and activator of transcription
## 7503: Human immunodeficiency virus 1 Trans-activation response element (TAR)
##                                             V11     V12         V13
##    1: A0A024A2C9_HAEIF|lph|tbp2|ERS515279_00376 protein   taxon:727
##    2: A0A024A2C9_HAEIF|lph|tbp2|ERS515279_00376 protein   taxon:727
##    3: A0A024A2C9_HAEIF|lph|tbp2|ERS515279_00376 protein   taxon:727
##    4: A0A024A2C9_HAEIF|lph|tbp2|ERS515279_00376 protein   taxon:727
##    5: A0A024A2C9_HAEIF|lph|tbp2|ERS515279_00376 protein   taxon:727
##   ---                                                              
## 7499:                              X2KN52_SALSA protein  taxon:8030
## 7500:                              X2KN52_SALSA protein  taxon:8030
## 7501:                              X2KN52_SALSA protein  taxon:8030
## 7502:                              X2KN52_SALSA protein  taxon:8030
## 7503:                                             miRNA taxon:11676
##            V14              V15                              V16 V17
##    1: 20170213           IntAct                                     
##    2: 20160817           AgBase positively_regulates(GO:0050840)    
##    3: 20170213           IntAct                                     
##    4: 20160817           AgBase                                     
##    5: 20170213           IntAct                                     
##   ---                                                               
## 7499: 20160315           AgBase                                     
## 7500: 20160315           AgBase                                     
## 7501: 20160315           AgBase                                     
## 7502: 20160315           AgBase                                     
## 7503: 20160607 ParkinsonsUK-UCL
```

```r
GO_IPI_publications = GO_IPI[,.(GO_IPI_pmids = unique(V6))]

write.table(GO_IPI_publications, "./results/GO_IPI_pmids.txt", quote=F, sep ="\t", row.names = F, col.names = T)
```

2288 publications are curated into GO IPI annotations 

#### Compare GO IPI publications to IMEx 


```r
imex = fread("https://raw.githubusercontent.com/pporrasebi/darkspaceproject/master/IMEx/results/imex_full.txt", header = T, sep = "\t", colClasses = "character")
N_pubid_imex = length(imex[,unique(pubid)])
N_pubid_GO_IPI = length(GO_IPI_publications$GO_IPI_pmids)
N_pubid_overlap = sum(!is.na(match(GO_IPI_publications$GO_IPI_pmids, imex[,unique(pubid)])))

venn.d = draw.pairwise.venn(area1 = N_pubid_imex, area2 = N_pubid_GO_IPI, cross.area = N_pubid_overlap, category = c("IMEx", "GO_IPI"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(0, 135), 
                          cat.dist = rep(0.035, 2), 
                          cat.cex = c(1,1), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          direct.area = TRUE,
                          cex = 1)
```

![](GO_IPI_dsgen_files/figure-html/GO_IPI_vs_imex_pmids-1.png)<!-- -->


I save a table of interacting pairs, publication IDs and BioGRID tag.


```r
#biogrid_from_mentha = fread("./processed_files/biogrid_pairs.txt", header = T, sep = "\t", colClasses = "character")
#fwrite(x = unique(biogrid_from_mentha[, .(pair_id_clean, pubid, biogrid = rep(1, .N))]), 
#       file = "./results/pairs_pmids_biogrid.txt", sep = "\t")
#N_biogrid = length(biogrid_from_mentha[,unique(pair_id_clean)])
```

The BioGRID dataset contains  interacting pairs. 

#### Compare BioGRID interactions and publications to IMEx 

I calculate how many interactions in BioGRID match to IMEx.


```r
#imex = fread("https://raw.githubusercontent.com/pporrasebi/darkspaceproject/master/IMEx/results/imex_full.txt", header = T, sep = "\t", colClasses = "character")
#N_imex = length(imex[,unique(pair_id_clean)])
#N_biogrid = length(biogrid_from_mentha[,unique(pair_id_clean)])
#N_overlap = sum(!is.na(match(biogrid_from_mentha[,unique(pair_id_clean)], imex[,unique(pair_id_clean)])))

#venn.d = draw.pairwise.venn(area1 = N_imex, area2 = N_biogrid, cross.area = N_overlap, category = c("IMEx", "BioGRID"), 
#                          lty = rep("blank", 2), 
#                          fill = c("blue", "red"), 
#                          alpha = rep(0.5, 2), cat.pos = c(0, 135), 
#                          cat.dist = rep(0.035, 2), 
#                          cat.cex = c(1,1), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
#                          direct.area = TRUE,
#                          cex = 1)
```

I calculate how many publications in BioGRID match to IMEx.


```r
#N_pubid_imex = length(imex[,unique(pubid)])
#N_pubid_biogrid = length(biogrid_from_mentha[,unique(pubid)])
#N_pubid_overlap = sum(!is.na(match(biogrid_from_mentha[,unique(pubid)], imex[,unique(pubid)])))

#venn.d = draw.pairwise.venn(area1 = N_pubid_imex, area2 = N_pubid_biogrid, cross.area = N_pubid_overlap, category = c("IMEx", "BioGRID"), 
#                          lty = rep("blank", 2), 
#                          fill = c("blue", "red"), 
#                          alpha = rep(0.5, 2), cat.pos = c(0, 135), 
#                          cat.dist = rep(0.035, 2), 
#                          cat.cex = c(1,1), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
#                          direct.area = TRUE,
#                          cex = 1)
```

I calculate how many interactions published in specific articles (the same interaction can come from different publications) in BioGRID match to IMEx.


```r
#N_pub_int_imex = length(imex[,unique(paste0(pubid,"_",pair_id_clean))])
#N_pub_int_biogrid = length(biogrid_from_mentha[,unique(paste0(pubid,"_",pair_id_clean))])
#N_pub_int_overlap = sum(!is.na(match(biogrid_from_mentha[,unique(paste0(pubid,"_",pair_id_clean))], imex[,unique(paste0(pubid,"_",pair_id_clean))])))

#venn.d = draw.pairwise.venn(area1 = N_pub_int_imex, area2 = N_pub_int_biogrid, cross.area = N_pub_int_overlap, category = c("IMEx", "BioGRID"), 
#                          lty = rep("blank", 2), 
#                          fill = c("blue", "red"), 
#                          alpha = rep(0.5, 2), cat.pos = c(0, 135), 
#                          cat.dist = rep(0.035, 2), 
#                          cat.cex = c(1,1), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
#                          direct.area = TRUE,
#                          cex = 1)
```
