---
title: Reactome dataset generation
date: 2016-04-13
author: Pablo Porras
---

Generation of the Reactome dataset: Technical report
========================================================

### Synopsis

The goal of this analysis is to use Reactome to try to find potentially non-curated protein interactions for which there is experimental evidence. I define a set of associated pairs as found in Reactome. Those proteins that are found to be present in the same reaction are combined according to their role.  

### Part 1: Loading Reactome data

The files used for this analysis were supplied by Antonio Fabregat at Reactome around end of June, 2015. I cannot use the interactions as exported by Reactome since they do not give a PMID for each inferred interacting pair and are then not usable. As a pre-load operation, I clean up the UniProtKB accessions given by Reactome to avoid using isoform identifiers. 


```r
system("perl ./scripts/reactome_prot_cleaner.pl ./source_files/rr.proteins.txt ./processed_files/rr.proteins.clean.txt")
rr_prots <- read.delim("./processed_files/rr.proteins.clean.txt", header = F, sep = "\t", colClasses = "character")
system("rm ./processed_files/rr.proteins.clean.txt")
colnames(rr_prots) <- c("upac","upac_clean","rr_id", "type")
library(dplyr)

rr_prots <- unique(select(rr_prots,upac=upac_clean,rr_id,type))

rr_pmids <- read.delim("./source_files/rr.pubmed.txt", header = F, sep = "\t", colClasses = "character")
colnames(rr_pmids) <- c("rr_id", "pmid")

rr_upacs <- unique(data.frame(rr_prots$upac))
```

```r
write.table(rr_upacs, "./processed_files/rr.upacs.txt", quote = F, row.names = F, col.names = F)
```

### Part 2: Select only human proteins

I need to query uniprot via a script as defined on their webpage (see [here](http://www.uniprot.org/help/programmatic_access#batch_retrieval_of_entries)), in order to get the taxIDs of the different organisms represented in the dataset. I will restrict the search to human interactions, so it is required to limit the search to those. 


```r
system("perl ./scripts/up_batch_retr.pl ./processed_files/rr.upacs.txt > ./processed_files/rr.up_taxid.txt")     # This takes a LONG time, do not run unless 100% needed. 
```

The previous script takes a long time and cannot actually deal with isoforms, so I manually perform the query on the UniProtKB website using the retrieve/ID mapping tool. This allows for isoform mapping, but I also need to clean the results up so each UniProtKB accession is listed in a single line. 


```r
system("perl ./scripts/multiplier.pl ./processed_files/rr.up_taxid_man.txt ./processed_files/rr.up_taxid_clean.txt")
```

```r
rr_up_taxid <- read.delim("./processed_files/rr.up_taxid_clean.txt", header = T, sep = "\t", colClasses = "character")
```

And I select human proteins only, which I use to select only those reactions involving human proteins from the original list. 


```r
rr_up_human <- data.frame(rr_up_taxid[rr_up_taxid$Organism.ID == 9606,1])
colnames(rr_up_human) <- c("upac")

rr_prots_human <- unique(merge(rr_up_human, rr_prots, by = "upac", all.x=T, all.y=F))

rr_reacts_human <- unique(data.frame(rr_prots_human$rr_id))
colnames(rr_reacts_human) <- "rr_id"

rr_pmids_human <- unique(merge(rr_reacts_human, rr_pmids, by="rr_id", all=F))

react_pmids_human <- unique(subset(rr_pmids_human, !is.na(rr_pmids_human$pmid), select = "pmid"))
react_pmids_human$reactome <- "yes"
```

### Part 3: Generate Reactome association pairs

First I classify the human proteins in reactome using their types, in order to construct pairs to compare to intact interactions. I consider only the role types 'INPUT', 'OUTPUT' and 'CATALYST'. 



```r
rr_p_h_input <- data.frame(rr_prots_human[rr_prots_human$type == "INPUT",1:2])
rr_p_h_output <- data.frame(rr_prots_human[rr_prots_human$type == "OUTPUT",1:2])
rr_p_h_catalyst <- data.frame(rr_prots_human[rr_prots_human$type == "CATALYST",1:2])
```

I generate the pairs of potential interactions, following these rules:

  *in-out pairs: proteins in a pair are annotated as 'INPUT' and 'OUTPUT' of a reaction, respectively.   
  *in-cat pairs: proteins annotated as 'INPUT' and 'CATALYST' in a reaction.  
  *out-cat pairs: proteins annotated as 'OUTPUT' and 'CATALYST' in a reaction.  


```r
rr_pair_inout <- unique(merge(rr_p_h_input, rr_p_h_output, by="rr_id", all=F))
rr_pair_incat <- unique(merge(rr_p_h_input, rr_p_h_catalyst, by="rr_id", all=F))
rr_pair_outcat <- unique(merge(rr_p_h_output, rr_p_h_catalyst, by="rr_id", all=F))

rr_pairs_1 <- unique(rbind(rr_pair_inout, rr_pair_incat))
rr_pairs <- unique(rbind(rr_pairs_1, rr_pair_outcat))
```

I finally add the PMID information to each pair of associations and I check if there are reactions for which no PMID is available. 


```r
rr_pairs_pmid_all <- unique(merge(rr_pairs, rr_pmids_human, by="rr_id",all=F))
```

### Part 4: Final formatting steps for generation of the Reactome dataset

I generate a unique pair identifier per line and save the Reactome pairs associated with the PMIDs as a text file for further comparison with other datasets.


```r
library(dplyr)

rr_pairs_pmid_all$pair_id <- apply(rr_pairs_pmid_all[,2:3], 1,function(i){
  paste(sort(i),collapse = "_")
})
reactome_pairs_pmid <- unique(select(rr_pairs_pmid_all, pair_id, pmid))
reactome_pairs_pmid$reactome <- 1
write.table(reactome_pairs_pmid,"./results/pairs_pmid_reactome.txt",col.names=T,row.names=F,quote=F,sep="\t")
system("gzip ./results/pairs_pmid_reactome.txt")
```

After collapsing redundant pairs, I end up with a list of 5180082 unique protein pair-publication associations extracted from Reactome in human. 
