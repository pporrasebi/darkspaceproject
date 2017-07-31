# Comparison of IMEx with other protein association datasets
Pablo Porras  
2017-07-28  

Estimating the size of the uncurated interactome
========================================================

### Synopsis

After producing tidy datasets comparing different resources to IMEx data, we put together the data and compare the overlap. 

IMPORTANT: This set of scripts assume the different referenced datasets have been freshly updated. If you need to update results of this set, please re-run the corresponding source as well. 

#### Required libraries

```r
# install packages if some are not already installed
packages = c("plyr","dplyr","data.table","UpSetR", "splitstackshape", "ggplot2", "shiny", "htmlwidgets")
if(mean(packages %in% installed.packages()) != 1){
        install.packages(packages[!packages %in% installed.packages()])
}
suppressPackageStartupMessages({
library(plyr)
library(dplyr)
library(data.table)
library(UpSetR)
library(splitstackshape)
library(ggplot2)
library(htmlwidgets)
})
```
### Part 1: Load datasets

#### Interaction datasets

##### IMEx dataset

I select only purely human interactions here (interactions where both proteins are human). 


```r
imex_full <- fread("../IMEx/results/imex_full.txt", header=T, sep="\t",colClasses="character",data.table = T)
imex_human <- unique(subset(imex_full,taxid_a=="9606" & taxid_b=="9606"))
imex_human$imex <- 1
imex_human_sel <- unique(select(imex_human,pair_id=pair_id_clean,pmid=pubid,imex))
imex_pairs <- unique(select(imex_human,pair_id=pair_id_clean,imex))
imex_pmids <- unique(select(imex_human,pmid=pubid,imex))
```

The dataset contains 158689 protein interactions recorded in 8988 publications.

##### BioGRID data


```r
setwd("../BioGRID/results/")
BioGRID_pairs_pmids <- fread("pairs_pmids_biogrid.txt",header=T,sep="\t",colClasses=c("character","character","numeric"),data.table = T)
setwd("../../dsp_comparison/")

BioGRID_pairs_pmids_sel <- unique(select(BioGRID_pairs_pmids, pair_id=pair_id_clean, pmid=pubid, BioGRID = biogrid))
BioGRID_pairs <- unique(select(BioGRID_pairs_pmids, pair_id=pair_id_clean,BioGRID = biogrid))
BioGRID_pmids <- unique(select(BioGRID_pairs_pmids, pmid=pubid,BioGRID = biogrid))
```

The BioGRID dataset contains 199566 protein associations recorded in 24241 publications.

##### GO IPI data (EBI_GOA_nonIntAct)


```r
setwd("../GO_IPI/results/")
GO_IPI_pairs_pmids <- fread("pairs_pmids_EBI_GOA_nonIntAct.txt",header=T,sep="\t",colClasses=c("character","character","numeric"),data.table = T)
setwd("../../dsp_comparison/")

GO_IPI_pairs_pmids_sel <- unique(select(GO_IPI_pairs_pmids, pair_id=pair_id_clean, pmid=pubid, GO_IPI = EBI_GOA_nonIntAct))
GO_IPI_pairs <- unique(select(GO_IPI_pairs_pmids, pair_id=pair_id_clean,GO_IPI = EBI_GOA_nonIntAct))
GO_IPI_pmids <- unique(select(GO_IPI_pairs_pmids, pmid=pubid,GO_IPI = EBI_GOA_nonIntAct))
```

The GO IPI dataset contains 9095 protein associations recorded in 5032 publications.

#### Pathway datasets

##### Reactome data


```r
setwd("../reactome_interactions/results/")
system("gunzip -k pairs_pmid_reactome.txt.gz")
reactome_pairs_pmids <- fread("pairs_pmid_reactome.txt",header=T,sep="\t",colClasses=c("character","character","numeric"),data.table = T)
system("rm pairs_pmid_reactome.txt")
setwd("../../dsp_comparison/")

reactome_pairs <- unique(select(reactome_pairs_pmids,pair_id,reactome))
reactome_pmids <- unique(select(reactome_pairs_pmids,pmid,reactome))
```

The reactome dataset contains 6419 protein associations recorded in 3046 publications. 
##### OmniPath interaction data 


```r
setwd("../OmniPath/results/")
OmniPath_interaction_pairs_pmids <- fread("pairs_pmids_OmniPath_interactions_minimal.txt",header=T,sep="\t",colClasses=c("character","character","character","character","character","numeric"),data.table = T)
setwd("../../dsp_comparison/")

OmniPath_interaction_pairs_pmids_sel <- unique(select(OmniPath_interaction_pairs_pmids, pair_id=pair_id_clean, pmid=pubid, OmniPath_interactions))
OmniPath_interaction_pairs <- unique(select(OmniPath_interaction_pairs_pmids, pair_id=pair_id_clean,OmniPath_interactions))
OmniPath_interaction_pmids <- unique(select(OmniPath_interaction_pairs_pmids, pmid=pubid,OmniPath_interactions))
```

The OmniPath-interactions dataset contains 5917 protein associations recorded in 6596 publications.

##### OmniPath ptm (post-translational modifications) data 


```r
setwd("../OmniPath/results/")
OmniPath_ptm_pairs_pmids <- fread("pairs_pmids_OmniPath_ptm_interactions_minimal.txt",header=T,sep="\t",colClasses=c("character","character","character","character","character","numeric"),data.table = T)
setwd("../../dsp_comparison/")

OmniPath_ptm_pairs_pmids_sel <- unique(select(OmniPath_ptm_pairs_pmids, pair_id=pair_id_clean, pmid=pubid, OmniPath_ptm))
OmniPath_ptm_pairs <- unique(select(OmniPath_ptm_pairs_pmids, pair_id=pair_id_clean,OmniPath_ptm))
OmniPath_ptm_pmids <- unique(select(OmniPath_ptm_pairs_pmids, pmid=pubid,OmniPath_ptm))
```

The OmniPath-ptms dataset contains 5596 protein associations recorded in 2609 publications.

##### Pathway-inferred STRING data

```r
string_pi <- fread("../STRING/results/pairs_STRING_pathway_inference.txt",header=T,colClasses=c("character","character","character","character","numeric","character","numeric"),data.table = T)
string_pi_pairs <- unique(string_pi[,.(STRING_pi_score_ave=mean(STRING_score),STRING_pi=STRING_pathway_inference),by=pair_id_clean])
string_pi_pairs <- string_pi_pairs[,.(pair_id=pair_id_clean,STRING_pi_score_ave,STRING_pi)]
```

The STRING pathway-inference dataset contains 191335 protein associations, no PMIDs provided in this case.

#### Text-mining datasets

##### Text-mining EPMC data


```r
setwd("../epmc_text_mining/results/")
system("gunzip -k pairs_pmids_tm.txt.gz")
tm_pairs_pmids <- fread("pairs_pmids_tm.txt",header=T,sep="\t",colClasses=c("character","character","numeric","character","numeric","numeric"),data.table = T)
system("rm pairs_pmids_tm.txt")
setwd("../../dsp_comparison/")

tm_pairs_pmids_sel <- unique(tm_pairs_pmids[,.(pair_id,pmid,tm_epmc=tm,tm_pr_times_found,tm_dm_times_found)])
tm_pairs <- unique(tm_pairs_pmids[,.(tm_epmc=tm,tm_pr_total_times_found=sum(tm_pr_times_found)),by=pair_id])
tm_pmids <- unique(tm_pairs_pmids[,.(pmid,tm_epmc=tm,tm_dm_times_found)])
```

The text-mining EPMC dataset contains 97705 protein associations recorded in 34048 publications.

##### Text-mining EVEX data

EVEX provides a confidence score, thehigher the score, the more likely the pair identified in a given publication is to be true. Negative values have no special meaning. 

```r
setwd("../EVEX/results/")
EVEX_pairs_pmids <- fread("pairs_pmids_EVEX_minimal.txt",header=T,sep="\t",colClasses=c("character","character","character","character","character","numeric","numeric"),data.table = T)
setwd("../../dsp_comparison/")

EVEX_pairs_pmids_sel <- unique(EVEX_pairs_pmids[, .(pair_id=pair_id_clean, pmid=pubid,EVEX,evex_score)])
EVEX_pairs <- unique(EVEX_pairs_pmids[, .(EVEX,evex_score_ave=mean(evex_score)),by=pair_id_clean])
EVEX_pairs <- EVEX_pairs[,.(pair_id=pair_id_clean,EVEX,evex_score_ave)]
EVEX_pmids <- unique(EVEX_pairs_pmids[,.(EVEX,evex_score_ave=mean(evex_score)),by=pubid])
EVEX_pmids <- EVEX_pmids[,.(pmid=pubid,EVEX,evex_score_ave)]
```

The text-mining EVEX dataset contains 540400 protein associations recorded in 54320 publications.

##### Text-mining STRING data

```r
string_tm <- fread("../STRING/results/pairs_STRING_textmining.txt",header=T,colClasses=c("character","character","character","character","numeric","character","numeric"),data.table = T)
string_tm_pairs <- unique(string_tm[,.(STRING_tm_score_ave=mean(STRING_score),STRING_textmining),by=pair_id_clean])
string_tm_pairs <- string_tm_pairs[,.(pair_id=pair_id_clean,STRING_tm_score_ave,STRING_textmining)]
```

The STRING text-mining dataset contains 167668 protein associations, no PMIDs provided in this case.

#### Prediction datasets

##### IID predictions data


```r
iid_pred_pairs <- fread("../iid_predictions/results/pairs_iid_pred.txt",header=T,sep="\t",colClasses=c("character","numeric"),data.table=T)
```

The IID-predictions dataset contains 706228 protein associations.

##### STRING predictions data

```r
string_phylo <- fread("../STRING/results/pairs_STRING_phylo_predictions.txt",header=T,sep="\t",colClasses=c("character","character","character","character","numeric","character","numeric"),data.table = T)
string_phylo_pairs <- unique(string_phylo[,.(STRING_phylo_score_ave=mean(STRING_score),STRING_phylo=STRING_phylo_predictions),by=pair_id_clean])
string_phylo_pairs <- string_phylo_pairs[,.(pair_id=pair_id_clean,STRING_phylo_score_ave,STRING_phylo)]
```

The STRING phylogenetic predictions dataset contains 242503 protein associations.

### Part 2: Generating comparison dataset at the pair level


```r
# Code below generates (to allow any number and any column names) and evaluates this:
# paste0("all_df <- list(",paste0(grep("_pairs$", ls(), value = T), collapse= ","),")")
# eval(parse(text=paste0("all_df <- list(",paste0(grep("_pairs$", ls(), value = T), collapse= ","),")")))
# results in errors upstream

all_df <- list(imex_pairs,reactome_pairs,tm_pairs,iid_pred_pairs, EVEX_pairs, BioGRID_pairs, GO_IPI_pairs, OmniPath_interaction_pairs, OmniPath_ptm_pairs,string_phylo_pairs,string_pi_pairs,string_tm_pairs)

paircomp_table <- Reduce(function(...) merge(..., all=TRUE), all_df)

# I clean and replace all NAs if present.

paircomp_table_final <- paircomp_table
paircomp_table_final[is.na(paircomp_table_final <- paircomp_table)] <- 0
paircomp_table_final = unique(paircomp_table_final)
fwrite(paircomp_table_final,"./results/paircomp_table_final.txt",col.names=T,row.names=F,sep="\t",quote=F)
system("gzip ./results/paircomp_table_final.txt --force")
unlink("./results/paircomp_table_final.txt")
```

The comparison set gives a total number of 1937544 potentially interacting pairs, of which 1778855 (91.81%) are not curated in IMEx. 

### Part 3: Generating comparison dataset at the publication level


```r
allpub_df <- list(imex_pmids,reactome_pmids,tm_pmids, EVEX_pmids, BioGRID_pmids, GO_IPI_pmids, OmniPath_interaction_pmids, OmniPath_ptm_pmids)

pubcomp_table <- Reduce(function(...) merge(..., all=TRUE), allpub_df)

# I clean and replace all NAs if present.

pubcomp_table_final <- pubcomp_table
pubcomp_table_final[is.na(pubcomp_table_final <- pubcomp_table)] <- 0
fwrite(pubcomp_table_final,"./results/pubcomp_table_final.txt",col.names=T,row.names=F,sep="\t",quote=F)
system("gzip ./results/pubcomp_table_final.txt --force")
```

The comparison set gives a total number of 114374 publications, of which 105386 (92.14%) are not curated in IMEx. 

### Part 4: Generating comparison dataset at the pair level taking the publication into account


```r
allpubpair_df <- list(imex_human_sel,reactome_pairs_pmids,tm_pairs_pmids_sel, EVEX_pairs_pmids_sel, BioGRID_pairs_pmids_sel, GO_IPI_pairs_pmids_sel, OmniPath_interaction_pairs_pmids_sel, OmniPath_ptm_pairs_pmids_sel)

prepubpaircomp_table_1 <- Reduce(function(...) merge(..., by=c("pair_id","pmid"), all=TRUE), allpubpair_df)

paironly_df <- list(iid_pred_pairs,string_pi_pairs,string_tm_pairs)

prepubpaircomp_table_2 <- Reduce(function(...) merge(..., by=c("pair_id"), all=TRUE), paironly_df)

pubpaircomp_table <- unique(merge(prepubpaircomp_table_1,prepubpaircomp_table_2,by=c("pair_id"),all=T))

# I clean and replace all NAs if present.

pubpaircomp_table_form <- pubpaircomp_table
pubpaircomp_table_form[is.na(pubpaircomp_table_form <- pubpaircomp_table)] <- 0
```

The comparison set gives a total number of 2093151 protein association pairs, of which 1912240 (91.36%) are not curated in IMEx. In all these pairs the publication from which they were derived was also matched, so the overlaps and numbers differ from my previous comparisons. 

### Part 5: Checking how many of the pair/publication combos involved uncurated proteins

I use a list of proteins non-curated in IMEx produced by Vitalii Kleschevnikov (IntAct group). He produced several versions of this list, using the different versions of the UniProtKB as reference. 

##### Pre-formatting pair/puplication combo comparisons

I need to preformat the 'pubpaircomp_table_final' data frame to compare it with VK lists.


```r
pubpaircomp_table_check_pt1 <- pubpaircomp_table_form
pubpaircomp_table_check_pt1$prots <- pubpaircomp_table_form$pair_id
pubpaircomp_table_check_pt2 <- cSplit(pubpaircomp_table_check_pt1, "prots", sep = "_", direction = "long")
```

##### Upload IMEx non-curated protein lists

I do not use the list that considers isoforms, since they were not considered in any of the datasets used in the comparison. 

```r
noimex_spnoisof <- unique(fread("./imex_non_curated/Swissprot_without_isoforms_missing_in_IntAct.txt",header=F,colClasses = "character",data.table = F))

noimex_upnoisof <- unique(fread("./imex_non_curated/UniprotKB_without_isoforms_missing_in_IntAct.txt",header=F,colClasses = "character",data.table = F))
```

##### Check how many of the proteins are missing from IMEx

I will identify those pairs that have proteins that have not been 

```r
pubpaircomp_table_check_pt3 <- mutate(pubpaircomp_table_check_pt2,
                                  noncur_prot = 
                                          ifelse(
                                                  prots %in% noimex_upnoisof$V1,
                                                  1,
                                                  0))

pubpaircomp_noncur_pairs <- unique(subset(pubpaircomp_table_check_pt3,noncur_prot==1,select= c("pair_id","noncur_prot")))

pubpaircomp_table_check_pt4 <- unique(merge(pubpaircomp_table_form,pubpaircomp_noncur_pairs,by="pair_id",all.x = T,all.y = F))
pubpaircomp_table_final <- pubpaircomp_table_check_pt4
pubpaircomp_table_final[is.na(pubpaircomp_table_final <- pubpaircomp_table_check_pt4)] <- 0
        
write.table(pubpaircomp_table_final,"./results/pubpaircomp_table_final.txt",col.names=T,row.names=F,sep="\t",quote=F)

setwd("./results")
system("tar -czvf pubpaircomp_table_final.txt.tar.gz pubpaircomp_table_final.txt && rm pubpaircomp_table_final.txt")
setwd("../")
```

### Part 6: Test sets for Reactome and text-mined datasets

I generate a couple of random samples from the subset of the text-mined and Reactome datasets that is not represented in IMEx. These samples are to be checked by curators to estimate what is the percentage of true/false positives in these datasets. I decide to just provide the PMIDs and request curators to check if the PMID contains interactors and identify the proteins if that's so. 


```r
pubpairs_not_imex <- pubpaircomp_table_final[pubpaircomp_table_final$imex==0 & pubpaircomp_table_final$pair_id != "",]
pubpairs_not_imex_sel <- unique(select(pubpairs_not_imex,pmid,pair_id))

publs_not_imex <- unique(select(pubpairs_not_imex,pmid,imex,reactome,tm_epmc,iid_pred,noncur_prot))

set.seed(88)

reactpubl_not_imex <- publs_not_imex[publs_not_imex$reactome==1,]
reactpubl_sample <- reactpubl_not_imex[sample(1:nrow(reactpubl_not_imex),100),]
react_sample <- unique(merge(reactpubl_sample,pubpairs_not_imex_sel,by="pmid"))

tmpubl_not_imex <- publs_not_imex[publs_not_imex$tm_epmc==1,]
tmpubl_sample <- tmpubl_not_imex[sample(1:nrow(tmpubl_not_imex),500),]
tm_sample <- unique(merge(tmpubl_sample,pubpairs_not_imex_sel,by="pmid"))

write.table(reactpubl_sample,"./results/reactpubl_sample.txt",col.names = T,row.names = F,sep="\t",quote=F)
write.table(tmpubl_sample,"./results/tmpubl_sample.txt",col.names = T,row.names = F,sep="\t",quote=F)
```

Low-hanging fruit and Reactome/TM evaluation lists can be checked and accessed at https://docs.google.com/spreadsheets/d/1tL1HtVD3-BxHxKuXbIYhcFjmCptGVEOD5aFJCZw6CZk/edit?usp=sharing. 

********************************************************************************************
