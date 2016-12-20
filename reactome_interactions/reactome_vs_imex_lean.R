## ----set-options, echo=FALSE---------------------------------------------
options(width = 80)

## ----eval=T,message=FALSE,warning=FALSE----------------------------------
system("perl ./scripts/reactome_prot_cleaner.pl ./source_files/rr.proteins.txt ./processed_files/rr.proteins.clean.txt")
rr_prots <- read.delim("./processed_files/rr.proteins.clean.txt", header = F, sep = "\t", colClasses = "character")
system("rm ./processed_files/rr.proteins.clean.txt")
colnames(rr_prots) <- c("upac","upac_clean","type", "rr_id")
library(dplyr)

rr_prots <- unique(select(rr_prots,upac=upac_clean,rr_id,type))

rr_pmids <- read.delim("./source_files/rr.pubmed.txt", header = F, sep = "\t", colClasses = "character")
colnames(rr_pmids) <- c("rr_id", "pmid")

rr_upacs <- unique(data.frame(rr_prots$upac))

## ----eval=FALSE----------------------------------------------------------
write.table(rr_upacs, "./processed_files/rr.upacs.txt", quote = F, row.names = F, col.names = F)

## ----eval=FALSE----------------------------------------------------------
## system("perl ./scripts/up_batch_retr.pl ./processed_files/rr.upacs.txt > ./processed_files/rr.up_taxid.txt")     # This takes a LONG time, do not run unless 100% needed.

## ----eval=FALSE----------------------------------------------------------
## system("perl ./scripts/multiplier.pl ./processed_files/rr.up_taxid_man.txt ./processed_files/rr.up_taxid_clean.txt")

## ------------------------------------------------------------------------
rr_up_taxid <- read.delim("./processed_files/rr.up_taxid_clean.txt", header = T, sep = "\t", colClasses = "character")

## ----eval=TRUE-----------------------------------------------------------
rr_up_human <- data.frame(rr_up_taxid[rr_up_taxid$Organism.ID == 9606,1])
colnames(rr_up_human) <- c("upac")

rr_prots_human <- unique(merge(rr_up_human, rr_prots, by = "upac", all.x=T, all.y=F))

rr_reacts_human <- unique(data.frame(rr_prots_human$rr_id))
colnames(rr_reacts_human) <- "rr_id"

rr_pmids_human <- unique(merge(rr_reacts_human, rr_pmids, by="rr_id", all=F))

react_pmids_human <- unique(subset(rr_pmids_human, !is.na(rr_pmids_human$pmid), select = "pmid"))
react_pmids_human$reactome <- "yes"

## ----eval=TRUE-----------------------------------------------------------
rr_p_h_input <- data.frame(rr_prots_human[rr_prots_human$type == "input",c(1,3)])
rr_p_h_output <- data.frame(rr_prots_human[rr_prots_human$type == "output",c(1,3)])
rr_p_h_catalyst <- data.frame(rr_prots_human[rr_prots_human$type == "catalystActivity",c(1,3)])

## ----eval=TRUE-----------------------------------------------------------
rr_pair_inout <- unique(merge(rr_p_h_input, rr_p_h_output, by="rr_id", all=F))
rr_pair_incat <- unique(merge(rr_p_h_input, rr_p_h_catalyst, by="rr_id", all=F))
rr_pair_outcat <- unique(merge(rr_p_h_output, rr_p_h_catalyst, by="rr_id", all=F))

rr_pairs_1 <- unique(rbind(rr_pair_inout, rr_pair_incat))
rr_pairs <- unique(rbind(rr_pairs_1, rr_pair_outcat))

## ------------------------------------------------------------------------
rr_pairs_pmid_all <- unique(merge(rr_pairs, rr_pmids_human, by="rr_id",all=F))

## ----message=FALSE,warning=FALSE-----------------------------------------
library(dplyr)

rr_pairs_pmid_all$pair_id <- apply(rr_pairs_pmid_all[,2:3], 1,function(i){
  paste(sort(i),collapse = "_")
})
reactome_pairs_pmid <- unique(select(rr_pairs_pmid_all, pair_id, pmid))
reactome_pairs_pmid$reactome <- 1
write.table(reactome_pairs_pmid,"./results/pairs_pmid_reactome.txt",col.names=T,row.names=F,quote=F,sep="\t")
system("gzip ./results/pairs_pmid_reactome.txt")

