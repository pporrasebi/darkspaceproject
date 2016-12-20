## ----set-options, echo=FALSE---------------------------------------------
options(width = 80)

## ----libraries,message=FALSE,warning=FALSE-------------------------------
library(dplyr)
library(data.table)
library(UpSetR)
library(splitstackshape)

## ----warning=FALSE,message=FALSE-----------------------------------------
imex_full <- fread("../IMEx/results/imex_full.txt", header=T, sep="\t",colClasses="character",data.table = F)
imex_human <- unique(subset(imex_full,taxid_a=="9606" & taxid_b=="9606"))
imex_human$imex <- 1
imex_human_sel <- unique(select(imex_human,pair_id=pair_id_clean,pmid=pubid,imex))
imex_pairs <- unique(select(imex_human,pair_id=pair_id_clean,imex))
imex_pmids <- unique(select(imex_human,pmid=pubid,imex))

## ------------------------------------------------------------------------
setwd("../reactome_interactions/results/")
system("gunzip -k pairs_pmid_reactome.txt.gz")
reactome_pairs_pmids <- fread("pairs_pmid_reactome.txt",header=T,sep="\t",colClasses=c("character","character","numeric"),data.table = F)
system("rm pairs_pmid_reactome.txt")
setwd("../../dsp_comparison/")

reactome_pairs <- unique(select(reactome_pairs_pmids,pair_id,reactome))
reactome_pmids <- unique(select(reactome_pairs_pmids,pmid,reactome))

## ------------------------------------------------------------------------
setwd("../epmc_text_mining/results/")
system("gunzip -k pairs_pmids_tm.txt.gz")
tm_pairs_pmids <- fread("pairs_pmids_tm.txt",header=T,sep="\t",colClasses=c("character","character","numeric"),data.table = F)
system("rm pairs_pmids_tm.txt")
setwd("../../dsp_comparison/")

tm_pairs_pmids_sel <- unique(select(tm_pairs_pmids,pair_id,pmid,tm_epmc=tm))
tm_pairs <- unique(select(tm_pairs_pmids,pair_id,tm_epmc=tm))
tm_pmids <- unique(select(tm_pairs_pmids,pmid,tm_epmc=tm))

## ------------------------------------------------------------------------
iid_pred_pairs <- fread("../iid_predictions/results/pairs_iid_pred.txt",header=T,sep="\t",colClasses=c("character","numeric"),data.table=F)

## ------------------------------------------------------------------------
all_df <- list(imex_pairs,reactome_pairs,tm_pairs,iid_pred_pairs)

comp_table <- Reduce(function(...) merge(..., all=TRUE), all_df)

# I clean and replace all NAs if present.

comp_table_final <- comp_table
comp_table_final[is.na(comp_table_final <- comp_table)] <- 0
write.table(comp_table_final,"./results/comp_table_final.txt",col.names=T,row.names=F,sep="\t",quote=F)
system("gzip ./results/comp_table_final.txt")

## ----echo=FALSE, fig.width=10,fig.height=6-------------------------------
upset(comp_table_final, 
      nsets = 4, 
      point.size = 6, 
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
              list(query = intersects, params = list("imex","iid_pred"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex","reactome"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex","tm_epmc"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex"), color= "gray70",active = T),
              list(query = intersects, params = list("reactome"), color= "gray70",active = T),
              list(query = intersects, params = list("iid_pred"), color= "gray70",active = T),
              list(query = intersects, params = list("tm_epmc"), color= "gray70",active = T)))


## ----echo=FALSE,fig.width=10,fig.height=6--------------------------------
upset(comp_table_final, 
      nsets = 4, 
      point.size = 6, 
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
              list(query = intersects, params = list("imex","iid_pred"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex","reactome"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex","tm_epmc"), color= "lightgoldenrod2",active = T)),
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

## ------------------------------------------------------------------------
allpub_df <- list(imex_pmids,reactome_pmids,tm_pmids)

pubcomp_table <- Reduce(function(...) merge(..., all=TRUE), allpub_df)

# I clean and replace all NAs if present.

pubcomp_table_final <- pubcomp_table
pubcomp_table_final[is.na(pubcomp_table_final <- pubcomp_table)] <- 0
write.table(pubcomp_table_final,"./results/pubcomp_table_final.txt",col.names=T,row.names=F,sep="\t",quote=F)
system("gzip ./results/pubcomp_table_final.txt")

## ----echo=FALSE, fig.width=10,fig.height=6-------------------------------
upset(pubcomp_table_final, 
      nsets = 3, 
      point.size = 6, 
      line.size = 2, 
      mainbar.y.label = "Common publications", 
      sets.x.label = "Nr of publications in dataset", 
      order.by="freq",
      decreasing=FALSE,
      queries = list(
              list(query = intersects, params = list("reactome","tm_epmc"), color = "blue", active = T),
              list(query = intersects, params = list("imex","reactome","tm_epmc"), color= "darkorange2",active = T),
              list(query = intersects, params = list("imex","reactome"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex","tm_epmc"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex"), color= "gray70",active = T),
              list(query = intersects, params = list("reactome"), color= "gray70",active = T),
              list(query = intersects, params = list("tm_epmc"), color= "gray70",active = T)))

## ----echo=FALSE, fig.width=10,fig.height=6-------------------------------
upset(pubcomp_table_final, 
      nsets = 3, 
      point.size = 6, 
      line.size = 2, 
      mainbar.y.label = "Common protein pairs", 
      sets.x.label = "Nr of protein pairs in dataset",
      order.by="freq",
      decreasing=FALSE,
      queries = list(
              list(query = intersects, params = list("reactome","tm_epmc"), color = "blue", active = T),
              list(query = intersects, params = list("imex","reactome","tm_epmc"), color= "darkorange2",active = T),
              list(query = intersects, params = list("imex","reactome"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex","tm_epmc"), color= "lightgoldenrod2",active = T)),
      intersections = list(
              list("reactome","tm_epmc"),
              list("imex","reactome","tm_epmc"),
              list("imex","reactome"),
              list("imex","tm_epmc")))

## ------------------------------------------------------------------------
prepubpaircomp_table_1 <- unique(merge(imex_human_sel,reactome_pairs_pmids,by=c("pair_id","pmid"),all=T))
prepubpaircomp_table_2 <- unique(merge(prepubpaircomp_table_1,tm_pairs_pmids_sel,by=c("pair_id","pmid"),all=T))
pubpaircomp_table <- unique(merge(prepubpaircomp_table_2,iid_pred_pairs,by=c("pair_id"),all=T))

# I clean and replace all NAs if present.

pubpaircomp_table_form <- pubpaircomp_table
pubpaircomp_table_form[is.na(pubpaircomp_table_form <- pubpaircomp_table)] <- 0

## ----echo=FALSE, fig.width=10,fig.height=6-------------------------------
upset(pubpaircomp_table_form, 
      nsets = 4, 
      point.size = 6, 
      line.size = 2, 
      mainbar.y.label = "Common pair/publication combos", 
      sets.x.label = "Nr of pair/publication combos in dataset", 
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
              list(query = intersects, params = list("imex","iid_pred"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex","reactome"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex","tm_epmc"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex"), color= "gray70",active = T),
              list(query = intersects, params = list("reactome"), color= "gray70",active = T),
              list(query = intersects, params = list("iid_pred"), color= "gray70",active = T),
              list(query = intersects, params = list("tm_epmc"), color= "gray70",active = T)))

## ----echo=FALSE,fig.width=10,fig.height=6--------------------------------
upset(pubpaircomp_table_form, 
      nsets = 4, 
      point.size = 6, 
      line.size = 2, 
      mainbar.y.label = "Common pair/publication combos", 
      sets.x.label = "Nr of pair/publication combos in dataset",
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
              list(query = intersects, params = list("imex","iid_pred"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex","reactome"), color= "lightgoldenrod2",active = T),
              list(query = intersects, params = list("imex","tm_epmc"), color= "lightgoldenrod2",active = T)),
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

## ----pp_combo_preformat--------------------------------------------------
pubpaircomp_table_check_pt1 <- pubpaircomp_table_form
pubpaircomp_table_check_pt1$prots <- pubpaircomp_table_form$pair_id
pubpaircomp_table_check_pt2 <- cSplit(pubpaircomp_table_check_pt1, "prots", sep = "_", direction = "long")

## ----noncur_upload-------------------------------------------------------
noimex_spnoisof <- unique(fread("./imex_non_curated/Swissprot_without_isoforms_missing_in_IntAct.txt",header=F,colClasses = "character",data.table = F))

noimex_upnoisof <- unique(fread("./imex_non_curated/UniprotKB_without_isoforms_missing_in_IntAct.txt",header=F,colClasses = "character",data.table = F))

## ----missing_check-------------------------------------------------------
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

## ----final_lhf-----------------------------------------------------------
lhf <- unique(subset(pubpaircomp_table_final,
                     imex == 0 &
                             reactome == 1 &
                             tm_epmc == 1 &
                             iid_pred == 1))

lhf <- arrange(lhf,desc(noncur_prot),desc(pmid))

lhf_pairs <- unique(select(lhf,pair_id))
lhf_pmids <- unique(select(lhf,pmid))

