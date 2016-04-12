---
title: Reactome comparison to IntAct
date: 2016-04-12
author: Pablo Porras
---

Estimating the size of the uncurated interactome
========================================================

### Synopsis

After producing tidy datasets comparing different resources to IMEx data, we put together the data and compare the overlap. 

### Part 1: Load datasets

#### IMEx dataset

I select only purely human interactions here (interactions where both proteins are human). 


```r
system("Rscript ./scripts/IMEx_ds_generator.R --save")
```

```r
imex_full <- read.delim("./results/imex_full.txt", header=T, sep="\t",colClasses="character")
imex_human <- unique(subset(imex_full,taxid_a=="9606" & taxid_b=="9606"))
library(dplyr)
imex_pairs <- unique(select(imex_human,pair_id=pair_id_clean))
imex_pairs$imex <- 1
```

#### Reactome data


```r
reactome_pairs <- read.csv("../reactome_interactions/results/pairs_reactome.txt",header=T,sep="\t",colClasses=c("character","numeric"))
```

#### Text-mining EPMC data


```r
tm_pairs <- read.csv("../epmc_text_mining/results/pairs_tm.txt",header=T,sep="\t",colClasses=c("character","numeric"))
```

#### IID predictions data


```r
iid_pred_pairs <- read.csv("../iid_predictions/results/pairs_iid_pred.txt",header=T,sep="\t",colClasses=c("character","numeric"))
```

### Part 2: Generating comparison dataset


```r
all_df <- list(imex_pairs,reactome_pairs,tm_pairs,iid_pred_pairs)

comp_table <- Reduce(function(...) merge(..., all=TRUE), all_df)

# I clean and replace all NAs if present.

comp_table_final <- comp_table
comp_table_final[is.na(omp_table_final <- comp_table)] <- 0
write.table(comp_table_final,"./results/comp_table_final.txt",col.names=T,row.names=F,sep="\t",quote=F)
```

The comparison set gives a total number of 1715496 potentially interacting pairs, of which 1599235 (93.22%) are not curated in IMEx. 

I finally produce a plot with the summary of the overlap between the different datasets evaluated. 

##### Figure 1: Comparison between different protein association datasets
![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

I also produce a sub-plot showing exclusively the intersections, so a more informative scale is used. 

##### Figure 2: Intersection between different protein association datasets
![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)


********************************************************************************************

### Appendix

#### Code for the figures

##### Figure 1: Intersection between different protein association datasets

```r
library(UpSetR)
upset(comp_table_final, 
      nsets = 4, 
      point.size = 6, 
      name.size = 12, 
      line.size = 2, 
      mainbar.y.label = "Common protein pairs", 
      sets.x.label = "Nr of protein pairs in dataset", 
      order.by="freq",
      decreasing=FALSE,
      queries = list(
              list(query = intersects, params = list("reactome","iid_pred","tm_epmc"), color = "blue", active = T),
              list(query = intersects, params = list("iid_pred","tm_epmc"), color= "cornflowerblue",active = T),
              list(query = intersects, params = list("reactome","iid_pred"), color= "cornflowerblue",active = T),
              list(query = intersects, params = list("reactome","tm_epmc"), color= "cornflowerblue",active = T),
              list(query = intersects, params = list("imex","iid_pred","reactome","tm_epmc"), color= "darkorange2",active = T),
              list(query = intersects, params = list("imex","reactome","tm_epmc"), color= "orange1",active = T),
              list(query = intersects, params = list("imex","iid_pred","tm_epmc"), color= "orange1",active = T),
              list(query = intersects, params = list("imex","iid_pred","reactome"), color= "orange1",active = T),
              list(query = intersects, params = list("imex","iid_pred"), color= "lightgoldenrod1",active = T),
              list(query = intersects, params = list("imex","reactome"), color= "lightgoldenrod1",active = T),
              list(query = intersects, params = list("imex","tm_epmc"), color= "lightgoldenrod1",active = T),
              list(query = intersects, params = list("imex"), color= "gray70",active = T),
              list(query = intersects, params = list("reactome"), color= "gray70",active = T),
              list(query = intersects, params = list("iid_pred"), color= "gray70",active = T),
              list(query = intersects, params = list("tm_epmc"), color= "gray70",active = T)))
```
##### Figure 2: Intersection between different protein association datasets (intersections only)

```r
upset(comp_table_final, 
      nsets = 4, 
      point.size = 6, 
      name.size = 12, 
      line.size = 2, 
      mainbar.y.label = "Common protein pairs", 
      sets.x.label = "Nr of protein pairs in dataset",
      order.by="freq",
      decreasing=FALSE,
      queries = list(
              list(query = intersects, params = list("reactome","iid_pred","tm_epmc"), color = "blue", active = T),
              list(query = intersects, params = list("iid_pred","tm_epmc"), color= "cornflowerblue",active = T),
              list(query = intersects, params = list("reactome","iid_pred"), color= "cornflowerblue",active = T),
              list(query = intersects, params = list("reactome","tm_epmc"), color= "cornflowerblue",active = T),
              list(query = intersects, params = list("imex","iid_pred","reactome","tm_epmc"), color= "darkorange2",active = T),
              list(query = intersects, params = list("imex","reactome","tm_epmc"), color= "orange1",active = T),
              list(query = intersects, params = list("imex","iid_pred","tm_epmc"), color= "orange1",active = T),
              list(query = intersects, params = list("imex","iid_pred","reactome"), color= "orange1",active = T),
              list(query = intersects, params = list("imex","iid_pred"), color= "lightgoldenrod1",active = T),
              list(query = intersects, params = list("imex","reactome"), color= "lightgoldenrod1",active = T),
              list(query = intersects, params = list("imex","tm_epmc"), color= "lightgoldenrod1",active = T)),
      intersections = list(
              list("reactome","iid_pred","tm_epmc"),
              list("iid_pred","tm_epmc"),
              list("reactome","iid_pred"),
              list("reactome","tm_epmc"),
              list("imex","iid_pred","reactome","tm_epmc"),
              list("imex","reactome","tm_epmc"),
              list("imex","iid_pred","tm_epmc"),
              list("imex","iid_pred","reactome"),
              list("imex","iid_pred"),
              list("imex","reactome"),
              list("imex","tm_epmc")))
```
