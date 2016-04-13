---
title: Evaluation text mining predictions EPMC
date: 2016-04-11
author: Pablo Porras
---


IID interaction predictions: Technical report
========================================================

### Synopsis

I take the fraction of IID (http://dcv.uhnres.utoronto.ca/iid/) data that is either predicted using the FpClass algorithm (http://www.ncbi.nlm.nih.gov/pubmed/25402006) or derived from orthologous interactions as a predictor of potentially interacting pairs we have missed. 

### Part 1: Download and process information from IID

I get the data from the IID web site.


```r
if(!file.exists("./source_files/iid.human.2016-03.txt.gz")){
        download.file("http://dcv.uhnres.utoronto.ca/downloads/iid.human.2016-03.txt.gz", destfile="./source_files/iid.human.2016-03.txt.gz", method= "curl")
}

iidhumanall <- read.table(gzfile("./source_files/iid.human.2016-03.txt.gz"), header = T, sep = "\t")
```

I have a look at the data to check how much of it is 100% predicted. 


```r
table(iidhumanall$evidence.type)
```

```

             exp       exp, ortho exp, ortho, pred        exp, pred 
          203529             6389            10535            47571 
           ortho      ortho, pred             pred 
           36519             8291           598612 
```

I will discard all the interactions that have been only detected experimentally. 


```r
iid_pred <- subset(iidhumanall, evidence.type != "exp")
table(iid_pred$evidence.type)
```

```

             exp       exp, ortho exp, ortho, pred        exp, pred 
               0             6389            10535            47571 
           ortho      ortho, pred             pred 
           36519             8291           598612 
```

Now I select the fields I will need for the comparison. I keep the field in which IID codes the database of origin to assess if there are mapping issues in the comparison. 


```r
library(dplyr)
iid_pred_sel <- unique(select(iid_pred,uniprot1,uniprot2,dbs))
```

I generate pair identifiers for each interaction reported in the dataset.


```r
iid_pred_sel$pair_id <- apply(iid_pred_sel[,1:2], 1,function(i){
        paste(sort(i),collapse = "_")
})
iid_pred_sel$iid <- "yes"
```

I save this as a text file for further comparison with other datasets.


```r
iid_pred_pairs <- unique(select(iid_pred_sel,pair_id))
iid_pred_pairs$iid_pred <- 1
write.table(iid_pred_pairs,"./results//pairs_iid_pred.txt",col.names=T,row.names=F,quote=F,sep="\t")
```

### Part 2: Load in the IMEx dataset

Now I load the IMEx dataset.


```r
system("Rscript ./scripts/IMEx_ds_generator.R --save")
```

```r
imex_full <- read.delim("./results/imex_full.txt", header=T, sep="\t",colClasses="character")
```

### Part 3: IMEx-predicted IID comparison

I compare both sets and save the results for further use. 


```r
imex_sel <- unique(select(imex_full,pair_id_clean,pair_id_clean,id_a_clean,id_b_clean,taxid_a,taxid_b,pubid))
imex_sel$imex <- "yes"
imex_pairs <- unique(select(imex_sel, pair_id=pair_id_clean,imex))

comp <- unique(merge(iid_pred_sel,imex_pairs,by="pair_id",all=T))

comp <- mutate(comp, db_pair =
                       ifelse(iid == "yes" & is.na(imex), "iid",
                       ifelse(is.na(iid) & imex == "yes", "imex",
                       ifelse(iid == "yes" & imex == "yes", "iid & imex",
                       "check"))))

comp$db_pair <- as.factor(comp$db_pair)

table(comp$db_pair,useNA="ifany")
```

```

       iid iid & imex       imex 
    683549      24368     420373 
```

I perform an additional check to see if there are interactions that IID obtained from IMEx-complying databases that were not properly mapped with my approach. 


```r
faulty_maps_imex <- unique(subset(comp, (grepl("intact",comp$dbs) | grepl("mint",comp$dbs) | grepl("dip",comp$dbs) | grepl("innatedb",comp$dbs)) & db_pair=="iid"))
```

There are 3893 interacting pairs that should have been mapped to IMEx and were not, which represents 0.35% of the total amount of pairs compared. I need to investigate further why this happens. For now, I will ignore it and generate a provisional file for further comparison with other approaches. 


```r
comp_simple <- unique(select(comp, pair_id,db_pair))

write.table(comp_simple,"./results/pairs_iid_vs_imex.txt",col.names=T,row.names=F,quote=F,sep="\t")
```
