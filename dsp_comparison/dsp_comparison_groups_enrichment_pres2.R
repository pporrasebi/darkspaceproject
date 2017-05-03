## ----setup, include=FALSE------------------------------------------------
opts_chunk$set(cache=F, message=FALSE,warning=FALSE, echo=FALSE)

## ----libraries-----------------------------------------------------------
# install packages if some are not already installed
# packages = c("plyr","dplyr","data.table","UpSetR", "splitstackshape", "ggplot2", "shiny", "VennDiagram")
# if(mean(packages %in% installed.packages()) != 1){
#        install.packages(packages[!packages %in% installed.packages()])
#}
suppressPackageStartupMessages({
library(plyr)
library(dplyr)
library(data.table)
library(UpSetR)
library(splitstackshape)
library(ggplot2)
library(shiny)
library(VennDiagram)
})

## ----results='hide'------------------------------------------------------
imex_full <- fread("../IMEx/results/imex_full.txt", header=T, sep="\t",colClasses="character",data.table = T)
imex_human = imex_full[taxid_a == "9606" & taxid_b == "9606",]
imex_human$imex <- 1
imex_human_sel <- unique(select(imex_human,pair_id=pair_id_clean,pmid=pubid,imex))
imex_pairs <- unique(select(imex_human,pair_id=pair_id_clean,imex))
imex_pmids <- unique(select(imex_human,pmid=pubid,imex))

## ---- results='hide', warning=FALSE,message=FALSE, echo=FALSE, cache = F----
setwd("../reactome_interactions/results/")
system("gunzip -k pairs_pmid_reactome.txt.gz")
reactome_pairs_pmids <- fread("pairs_pmid_reactome.txt",header=T,sep="\t",colClasses=c("character","character","numeric"),data.table = T)
system("rm pairs_pmid_reactome.txt")
setwd("../../dsp_comparison/")

reactome_pairs <- unique(select(reactome_pairs_pmids,pair_id,reactome))
reactome_pmids <- unique(select(reactome_pairs_pmids,pmid,reactome))

## ---- fig.height=8, fig.width=8, results='hide'--------------------------
setwd("../epmc_text_mining/results/")
system("gunzip -k pairs_pmids_tm.txt.gz")
tm_pairs_pmids <- fread("pairs_pmids_tm.txt",header=T,sep="\t",colClasses=c("character", "character","numeric","character"),data.table = T)
system("rm pairs_pmids_tm.txt")
setwd("../../dsp_comparison/")

tm_pairs_pmids_sel <- unique(select(tm_pairs_pmids,pair_id,pmid,tm_epmc=tm))
tm_pairs <- unique(select(tm_pairs_pmids,pair_id,tm_epmc=tm))
tm_pmids <- unique(select(tm_pairs_pmids,pmid,tm_epmc=tm))

## ----results='hide'------------------------------------------------------
iid_pred_pairs <- fread("../iid_predictions/results/pairs_iid_pred.txt",header=T,sep="\t",colClasses=c("character","numeric"),data.table=T)

## ----results='hide'------------------------------------------------------
setwd("../EVEX/results/")
EVEX_pairs_pmids <- fread("pairs_pmids_EVEX_shallow.txt",header=T,sep="\t",colClasses=c("character","character","character","character","character","numeric"),data.table = T)
setwd("../../dsp_comparison/")

EVEX_pairs_pmids_sel <- unique(select(EVEX_pairs_pmids, pair_id=pair_id_clean, pmid=pubid,EVEX))
EVEX_pairs <- unique(select(EVEX_pairs_pmids, pair_id=pair_id_clean,EVEX))
EVEX_pmids <- unique(select(EVEX_pairs_pmids, pmid=pubid,EVEX))

## ----results='hide'------------------------------------------------------
setwd("../BioGRID/results/")
BioGRID_pairs_pmids <- fread("pairs_pmids_biogrid.txt",header=T,sep="\t",colClasses=c("character","character","numeric"),data.table = T)
setwd("../../dsp_comparison/")

BioGRID_pairs_pmids_sel <- unique(select(BioGRID_pairs_pmids, pair_id=pair_id_clean, pmid=pubid, BioGRID = biogrid))
BioGRID_pairs <- unique(select(BioGRID_pairs_pmids, pair_id=pair_id_clean,BioGRID = biogrid))
BioGRID_pmids <- unique(select(BioGRID_pairs_pmids, pmid=pubid,BioGRID = biogrid))

## ----results='hide'------------------------------------------------------
setwd("../GO_IPI/results/")
GO_IPI_pairs_pmids <- fread("pairs_pmids_EBI_GOA_nonIntAct.txt",header=T,sep="\t",colClasses=c("character","character","numeric"),data.table = T)
setwd("../../dsp_comparison/")

GO_IPI_pairs_pmids_sel <- unique(select(GO_IPI_pairs_pmids, pair_id=pair_id_clean, pmid=pubid, GO_IPI = EBI_GOA_nonIntAct))
GO_IPI_pairs <- unique(select(GO_IPI_pairs_pmids, pair_id=pair_id_clean,GO_IPI = EBI_GOA_nonIntAct))
GO_IPI_pmids <- unique(select(GO_IPI_pairs_pmids, pmid=pubid,GO_IPI = EBI_GOA_nonIntAct))

## ----results='hide'------------------------------------------------------
setwd("../OmniPath/results/")
OmniPath_interaction_pairs_pmids <- fread("pairs_pmids_OmniPath_interactions_minimal.txt",header=T,sep="\t",colClasses=c("character","character","character","character","character","numeric"),data.table = T)
setwd("../../dsp_comparison/")

setwd("../OmniPath/results/")
OmniPath_ptm_pairs_pmids <- fread("pairs_pmids_OmniPath_ptm_interactions_minimal.txt",header=T,sep="\t",colClasses=c("character","character","character","character","character","numeric"),data.table = T)
setwd("../../dsp_comparison/")

OmniPath_interaction_pairs_pmids_sel <- unique(select(OmniPath_interaction_pairs_pmids, pair_id=pair_id_clean, pmid=pubid, OmniPath_interactions))
OmniPath_interaction_pairs <- unique(select(OmniPath_interaction_pairs_pmids, pair_id=pair_id_clean,OmniPath_interactions))
OmniPath_interaction_pmids <- unique(select(OmniPath_interaction_pairs_pmids, pmid=pubid,OmniPath_interactions))

OmniPath_ptm_pairs_pmids_sel <- unique(select(OmniPath_ptm_pairs_pmids, pair_id=pair_id_clean, pmid=pubid, OmniPath_ptm))
OmniPath_ptm_pairs <- unique(select(OmniPath_ptm_pairs_pmids, pair_id=pair_id_clean,OmniPath_ptm))
OmniPath_ptm_pmids <- unique(select(OmniPath_ptm_pairs_pmids, pmid=pubid,OmniPath_ptm))

## ----omnipath_merge, results='hide'--------------------------------------
setnames(OmniPath_interaction_pairs_pmids, colnames(OmniPath_interaction_pairs_pmids)[6], "OmniPath")
setnames(OmniPath_ptm_pairs_pmids, colnames(OmniPath_ptm_pairs_pmids)[6], "OmniPath")
OmniPath_pairs_pmids = rbind(OmniPath_interaction_pairs_pmids, OmniPath_ptm_pairs_pmids)

OmniPath_pairs_pmids_sel <- unique(select(OmniPath_pairs_pmids, pair_id=pair_id_clean, pmid=pubid, OmniPath))
OmniPath_pairs <- unique(select(OmniPath_pairs_pmids, pair_id=pair_id_clean,OmniPath))
OmniPath_pmids <- unique(select(OmniPath_pairs_pmids, pmid=pubid,OmniPath))

## ----STRING, results='hide'----------------------------------------------
# read STRING_textmining table
setwd("../STRING/results/")
STRING_textmining <- fread("pairs_STRING_textmining.txt",header=T,sep="\t",colClasses=c("character","character","character","character","numeric","character","numeric"),data.table = T)
STRING_textmining_pairs <- unique(select(STRING_textmining, pair_id=pair_id_clean,STRING_textmining))

# read STRING_pathway_inference table
STRING_pathway_inference <- fread("pairs_STRING_pathway_inference.txt",header=T,sep="\t",colClasses=c("character","character","character","character","numeric","character","numeric"),data.table = T)
STRING_pathway_inference_pairs <- unique(select(STRING_pathway_inference, pair_id=pair_id_clean,STRING_pathway_inference))

# read STRING_phylo_predictions table
STRING_phylo_predictions <- fread("pairs_STRING_phylo_predictions.txt",header=T,sep="\t",colClasses=c("character","character","character","character","numeric","character","numeric"),data.table = T)
STRING_phylo_predictions_pairs <- unique(select(STRING_phylo_predictions, pair_id=pair_id_clean,STRING_phylo_predictions))

setwd("../../dsp_comparison/")


## ------------------------------------------------------------------------
# Code below generates (to allow any number and any column names) and evaluates this:
# paste0("all_df <- list(",paste0(grep("_pairs$", ls(), value = T), collapse= ","),")")
# eval(parse(text=paste0("all_df <- list(",paste0(grep("_pairs$", ls(), value = T), collapse= ","),")")))
# results in errors upstream

all_df <- list(imex_pairs,reactome_pairs,tm_pairs,iid_pred_pairs, EVEX_pairs, BioGRID_pairs, GO_IPI_pairs, OmniPath_pairs, STRING_textmining_pairs, STRING_pathway_inference_pairs, STRING_phylo_predictions_pairs)

comp_table <- Reduce(function(...) merge(..., all=TRUE), all_df)

# I clean and replace all NAs if present.

comp_table_final <- comp_table
comp_table_final[is.na(comp_table_final <- comp_table)] <- 0
comp_table_final = unique(comp_table_final)
fwrite(comp_table_final,"./results/comp_table_final.txt",col.names=T,row.names=F,sep="\t",quote=F)
system("gzip ./results/comp_table_final.txt --force")
unlink("./results/comp_table_final.txt")

## ----visualize_overlaps, fig.height=9, fig.width=12, results='hide'------
groups = colnames(comp_table_final)[2:ncol(comp_table_final)]
# generate expression like this "comp_table_final[, overlaps := paste0(imex,reactome,tm_epmc,iid_pred,EVEX,BioGRID,GO_IPI,OmniPath)]" and evaluate it - this way the code is independent of column names
eval(parse(text=paste0("comp_table_final[, overlaps := paste0(",paste0(groups, collapse= ","),")]")))
eval(parse(text=paste0("comp_table_final[, N_overlaps := sum(",paste0(groups, collapse= ","),"), by = pair_id]")))

comp_table_final[, N_per_overlaps := .N, by = overlaps]
comp_table_final[, N_per_N_overlaps := .N, by = N_overlaps]

N_overlaps = unique(comp_table_final[,.(N_overlaps, N_per_N_overlaps)])
N_overlaps = N_overlaps[order(N_per_N_overlaps, decreasing = T),]
qplot(label = N_overlaps$N_per_N_overlaps, y = N_overlaps$N_per_N_overlaps, x = N_overlaps$N_overlaps,  geom = "text") + scale_y_log10(labels = 10^c(1:7), breaks = 10^c(1:7)) +
        ggtitle("The number of interacting pairs and the number of resources \n the pairs are shared between") + ylab("the number of interacting pairs, log10 scale")+ xlab("the number of resources")+scale_x_continuous(labels = N_overlaps$N_overlaps, breaks = N_overlaps$N_overlaps)+theme(text = element_text(size =18), title = element_text(size =24), legend.position="none")

overlaps = unique(comp_table_final[,.(overlaps, N_per_overlaps)])
overlaps = overlaps[order(N_per_overlaps, decreasing = T),]

## ---- eval=T, results='hide'---------------------------------------------
allpub_df <- list(imex_pmids,reactome_pmids,tm_pmids, EVEX_pmids, BioGRID_pmids, GO_IPI_pmids, OmniPath_pmids)

pubcomp_table <- Reduce(function(...) merge(..., all=TRUE), allpub_df)

# I clean and replace all NAs if present.

pubcomp_table_final <- pubcomp_table
pubcomp_table_final[is.na(pubcomp_table_final <- pubcomp_table)] <- 0
fwrite(pubcomp_table_final,"./results/pubcomp_table_final.txt",col.names=T,row.names=F,sep="\t",quote=F)
system("gzip ./results/pubcomp_table_final.txt --force")

## ----visualize_overlaps_publications, fig.height=9, fig.width=12, results='hide'----
pubgroups = colnames(pubcomp_table_final)[2:ncol(pubcomp_table_final)]
pubcomp_table_final[, puboverlaps := paste0(imex, reactome, tm_epmc, EVEX, BioGRID, GO_IPI, OmniPath)]
pubcomp_table_final[, N_puboverlaps := sum(imex, reactome, tm_epmc, EVEX, BioGRID, GO_IPI, OmniPath), by = pmid]
pubcomp_table_final[, N_per_puboverlaps := .N, by = puboverlaps]
pubcomp_table_final[, N_per_N_puboverlaps := .N, by = N_puboverlaps]

N_puboverlaps = unique(pubcomp_table_final[,.(N_puboverlaps, N_per_N_puboverlaps)])
N_puboverlaps = N_puboverlaps[order(N_per_N_puboverlaps, decreasing = T),]
qplot(label = N_puboverlaps$N_per_N_puboverlaps, y = N_puboverlaps$N_per_N_puboverlaps, x = N_puboverlaps$N_puboverlaps,  geom = "text") + scale_y_log10(labels = 10^c(1:7), breaks = 10^c(1:7)) +
        ggtitle("The number of publications and the number of resources \n the publications are shared between") + ylab("the number of publications, log10 scale")+ xlab("the number of resources")+scale_x_continuous(labels = N_puboverlaps$N_puboverlaps, breaks = N_puboverlaps$N_puboverlaps)+theme(text = element_text(size =15), title = element_text(size =24))

## ---- eval=T, results='hide'---------------------------------------------
prepubpaircomp_table_1 <- unique(merge(imex_human_sel,reactome_pairs_pmids,by=c("pair_id","pmid"),all=T))
prepubpaircomp_table_1 <- unique(merge(prepubpaircomp_table_1,tm_pairs_pmids_sel,by=c("pair_id","pmid"),all=T))
prepubpaircomp_table_1 <- unique(merge(prepubpaircomp_table_1,EVEX_pairs_pmids_sel,by=c("pair_id","pmid"),all=T))
prepubpaircomp_table_1 <- unique(merge(prepubpaircomp_table_1,BioGRID_pairs_pmids_sel,by=c("pair_id","pmid"),all=T))
prepubpaircomp_table_1 <- unique(merge(prepubpaircomp_table_1,GO_IPI_pairs_pmids_sel,by=c("pair_id","pmid"),all=T))
prepubpaircomp_table_1 <- unique(merge(prepubpaircomp_table_1,OmniPath_pairs_pmids_sel,by=c("pair_id","pmid"),all=T))

pubpaircomp_table <- unique(merge(prepubpaircomp_table_1,iid_pred_pairs,by=c("pair_id"),all=T))

# I clean and replace all NAs if present.

pubpaircomp_table_form <- pubpaircomp_table
pubpaircomp_table_form[is.na(pubpaircomp_table_form)] <- 0

## ----Pathways_1, fig.width = 16, fig.height= 6---------------------------

grid.newpage()
plotname = paste0("Comparison of pathway resources: Reactome, Omnipath, STRING pathway-based inference")
# setup grid
pushViewport(viewport(layout=grid.layout(nrow = 3, ncol=5, widths = unit(c(1/14,4/14,4/14,4/14,1/14), "npc"), 
                                         heights = unit(c(1/7, 1/7, 5/7), "npc"))))
# plot names
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 1))
x =grid.text(plotname, x = unit(1.6, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

# individual labels
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 2))
x =grid.text("interacting pairs", x = unit(0.5, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

pushViewport(viewport(layout.pos.col=3, layout.pos.row = 2))
x =grid.text("publications", x = unit(0.5, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

pushViewport(viewport(layout.pos.col=4, layout.pos.row = 2))
x =grid.text("interacting pairs and publications", x = unit(0.5, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

# venn diagrams
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 3))

Pathways_1 = unique(comp_table_final[,.(pair_id, reactome, OmniPath, STRING_pathway_inference)])
venn_Pathways_1 = draw.triple.venn(area1 = Pathways_1[reactome == 1,.N], 
                          area2 = Pathways_1[OmniPath == 1,.N], 
                          area3 = Pathways_1[STRING_pathway_inference == 1,.N], 
                          n12 = Pathways_1[reactome == 1 & OmniPath == 1,.N], 
                          n23 = Pathways_1[OmniPath == 1 & STRING_pathway_inference == 1,.N],
                          n13 = Pathways_1[reactome == 1 & STRING_pathway_inference == 1,.N],
                          n123 = Pathways_1[reactome == 1 & OmniPath == 1 & STRING_pathway_inference == 1,.N], 
                          category = c("Reactome", 
                                       "OmniPath", 
                                       "STRING_pathway_inference"), 
                          lty = rep("blank", 3), 
                          fill = c("blue", "red", "green"), 
                          alpha = rep(0.5, 3), cat.pos = c(350, 25, 160), 
                          cat.dist = c(0.08,0.035, 0.08), 
                          cat.cex = c(1.4,1.4, 1.4), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 1.5
)
popViewport()

pushViewport(viewport(layout.pos.col=3, layout.pos.row = 3))
Pathways_2 = unique(pubcomp_table_final[,.(pmid, reactome, OmniPath)])
venn_Pathways_2 = draw.pairwise.venn(area1 = Pathways_2[reactome == 1,.N], 
                          area2 = Pathways_2[OmniPath == 1,.N], 
                          cross.area = Pathways_2[reactome == 1 & OmniPath == 1,.N], 
                          category = c("Reactome", 
                                       "OmniPath"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(350, 25), 
                          cat.dist = c(0.08,0.035), 
                          cat.cex = c(1.4,1.4), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 1.5
)
popViewport()

pushViewport(viewport(layout.pos.col=4, layout.pos.row = 3))
Pathways_3 = unique(pubpaircomp_table_form[,.(pair_id, pmid, reactome, OmniPath)])
venn_Pathways_3 = draw.pairwise.venn(area1 = Pathways_3[reactome == 1,.N], 
                          area2 = Pathways_3[OmniPath == 1,.N], 
                          cross.area = Pathways_3[reactome == 1 & OmniPath == 1,.N], 
                          category = c("Reactome", 
                                       "OmniPath"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(350, 25), 
                          cat.dist = c(0.08,0.035), 
                          cat.cex = c(1.4,1.4), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 1.5
)
popViewport()


## ----tm_1, fig.width = 16, fig.height= 6---------------------------------
{ 
grid.newpage()
plotname = paste0("Comparison of text-mining (NLP?) resources: EVEX and custom EPMC text-mining datasets")
# setup grid
pushViewport(viewport(layout=grid.layout(nrow = 3, ncol=5, widths = unit(c(1/14,4/14,4/14,4/14,1/14), "npc"), 
                                         heights = unit(c(1/7, 1/7, 5/7), "npc"))))
# plot names
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 1))
x =grid.text(plotname, x = unit(1.3, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

# individual labels
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 2))
x =grid.text("interacting pairs", x = unit(0.5, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

pushViewport(viewport(layout.pos.col=3, layout.pos.row = 2))
x =grid.text("publications", x = unit(0.5, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

pushViewport(viewport(layout.pos.col=4, layout.pos.row = 2))
x =grid.text("interacting pairs and publications", x = unit(0.5, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

# venn diagrams
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 3))
textmining_1 = unique(comp_table_final[,.(pair_id, tm_epmc, EVEX, STRING_textmining)])
venn_Interactions_1 = draw.triple.venn(area1 = textmining_1[tm_epmc == 1,.N], 
                          area2 = textmining_1[EVEX == 1,.N], 
                          area3 = textmining_1[STRING_textmining == 1,.N], 
                          n12 = textmining_1[tm_epmc == 1 & EVEX == 1,.N], 
                          n23 = textmining_1[EVEX == 1 & STRING_textmining == 1,.N],
                          n13 = textmining_1[tm_epmc == 1 & STRING_textmining == 1,.N],
                          n123 = textmining_1[tm_epmc == 1 & EVEX == 1 & STRING_textmining == 1,.N], 
                          category = c("EPMC text-mining", 
                                       "EVEX", 
                                       "STRING_textmining"), 
                          lty = rep("blank", 3), 
                          fill = c("blue", "red", "green"), 
                          alpha = rep(0.5, 3), cat.pos = c(350, 25, 160), 
                          cat.dist = c(0.08,0.035, 0.08), 
                          cat.cex = c(1.4,1.4, 1.4), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 1.5
)
popViewport()

pushViewport(viewport(layout.pos.col=3, layout.pos.row = 3))
textmining_2 = unique(pubcomp_table_final[,.(pmid, tm_epmc, EVEX)])
venn_Pathways_2 = draw.pairwise.venn(area1 = textmining_2[tm_epmc == 1,.N], 
                          area2 = textmining_2[EVEX == 1,.N], 
                          cross.area = textmining_2[tm_epmc == 1 & EVEX == 1,.N], 
                          category = c("EPMC text-mining", 
                                       "EVEX"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(350, 25), 
                          cat.dist = c(0.08,0.035), 
                          cat.cex = c(1.4,1.4), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 1.5
)
popViewport()

pushViewport(viewport(layout.pos.col=4, layout.pos.row = 3))
textmining_3 = unique(pubpaircomp_table_form[,.(pair_id, pmid, tm_epmc, EVEX)])
venn_Pathways_3 = draw.pairwise.venn(area1 = textmining_3[tm_epmc == 1,.N], 
                          area2 = textmining_3[EVEX == 1,.N], 
                          cross.area = textmining_3[tm_epmc == 1 & EVEX == 1,.N], 
                          category = c("EPMC text-mining", 
                                       "EVEX"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(350, 25), 
                          cat.dist = c(0.08,0.035), 
                          cat.cex = c(1.4,1.4), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 1.5
)
popViewport()
}

## ----Predictions_1,  fig.width = 16, fig.height= 6-----------------------
{ 
grid.newpage()
plotname = paste0("Comparison of predictions resources: IID and STRING phylogeny-/orthology-based predictions")
# setup grid
pushViewport(viewport(layout=grid.layout(nrow = 3, ncol=3, widths = unit(c(3/14,8/14,3/14), "npc"), 
                                         heights = unit(c(1/7, 1/7, 5/7), "npc"))))
# plot names
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 1))
x =grid.text(plotname, x = unit(0.3, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

# individual labels
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 2))
x =grid.text("interacting pairs", x = unit(0.5, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

# venn diagrams
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 3))
Predictions_1 = unique(comp_table_final[,.(pair_id, iid_pred, STRING_phylo_predictions)])
venn_Pathways_2 = draw.pairwise.venn(area1 = Predictions_1[iid_pred == 1,.N], 
                          area2 = Predictions_1[STRING_phylo_predictions == 1,.N], 
                          cross.area = Predictions_1[iid_pred == 1 & STRING_phylo_predictions == 1,.N], 
                          category = c("IID predictions", 
                                       "STRING predictions"), 
                          lty = rep("blank", 2), 
                          fill = c("blue", "red"), 
                          alpha = rep(0.5, 2), cat.pos = c(350, 25), 
                          cat.dist = c(0.08,0.035), 
                          cat.cex = c(1.4,1.4), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 1.5
)
popViewport()
}

## ----Interactions_1, fig.width = 16, fig.height= 6-----------------------

{
grid.newpage()
plotname = paste0("Comparison of interaction resources: IMEx, Biogrid, GO IPI")
# setup grid
pushViewport(viewport(layout=grid.layout(nrow = 3, ncol=5, widths = unit(c(1/14,4/14,4/14,4/14,1/14), "npc"), 
                                         heights = unit(c(1/7, 1/7, 5/7), "npc"))))
# plot names
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 1))
x =grid.text(plotname, x = unit(0.7, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

# individual labels
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 2))
x =grid.text("interacting pairs", x = unit(0.5, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

pushViewport(viewport(layout.pos.col=3, layout.pos.row = 2))
x =grid.text("publications", x = unit(0.5, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

pushViewport(viewport(layout.pos.col=4, layout.pos.row = 2))
x =grid.text("interacting pairs and publications", x = unit(0.5, "npc"),y= unit(0.5, "npc"), gp=gpar(fontsize=20))
popViewport()

# venn diagrams
pushViewport(viewport(layout.pos.col=2, layout.pos.row = 3))
Interactions_1 = unique(comp_table_final[,.(pair_id, imex, BioGRID, GO_IPI)])
N_Interactions_1 = nrow(Interactions_1)
venn_Interactions_1 = draw.triple.venn(area1 = Interactions_1[imex == 1,.N], 
                          area2 = Interactions_1[BioGRID == 1,.N], 
                          area3 = Interactions_1[GO_IPI == 1,.N], 
                          n12 = Interactions_1[imex == 1 & BioGRID == 1,.N], 
                          n23 = Interactions_1[BioGRID == 1 & GO_IPI == 1,.N],
                          n13 = Interactions_1[imex == 1 & GO_IPI == 1,.N],
                          n123 = Interactions_1[imex == 1 & BioGRID == 1 & GO_IPI == 1,.N], 
                          category = c("IMEx", 
                                       "BioGRID", 
                                       "GO_IPI"), 
                          lty = rep("blank", 3), 
                          fill = c("blue", "red", "green"), 
                          alpha = rep(0.5, 3), cat.pos = c(350, 25, 160), 
                          cat.dist = c(0.08,0.035, 0.08), 
                          cat.cex = c(1.4,1.4, 1.4), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 1.5
)
popViewport()

pushViewport(viewport(layout.pos.col=3, layout.pos.row = 3))
Interactions_2 = unique(pubcomp_table_final[,.(pmid, imex, BioGRID, GO_IPI)])
N_Interactions_2 = nrow(Interactions_2)
venn_Interactions_2 = draw.triple.venn(area1 = Interactions_2[imex == 1,.N], 
                          area2 = Interactions_2[BioGRID == 1,.N], 
                          area3 = Interactions_2[GO_IPI == 1,.N], 
                          n12 = Interactions_2[imex == 1 & BioGRID == 1,.N], 
                          n23 = Interactions_2[BioGRID == 1 & GO_IPI == 1,.N],
                          n13 = Interactions_2[imex == 1 & GO_IPI == 1,.N],
                          n123 = Interactions_2[imex == 1 & BioGRID == 1 & GO_IPI == 1,.N], 
                          category = c("IMEx", 
                                       "BioGRID", 
                                       "GO_IPI"), 
                          lty = rep("blank", 3), 
                          fill = c("blue", "red", "green"), 
                          alpha = rep(0.5, 3), cat.pos = c(350, 25, 160), 
                          cat.dist = c(0.08,0.035, 0.08), 
                          cat.cex = c(1.4,1.4, 1.4), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 1.5
)
popViewport()

pushViewport(viewport(layout.pos.col=4, layout.pos.row = 3))
Interactions_3 = unique(pubpaircomp_table_form[,.(pair_id, pmid, imex, BioGRID, GO_IPI)])
N_Interactions_3 = nrow(Interactions_3)
venn_Interactions_3 = draw.triple.venn(area1 = Interactions_3[imex == 1,.N], 
                          area2 = Interactions_3[BioGRID == 1,.N], 
                          area3 = Interactions_3[GO_IPI == 1,.N], 
                          n12 = Interactions_3[imex == 1 & BioGRID == 1,.N], 
                          n23 = Interactions_3[BioGRID == 1 & GO_IPI == 1,.N],
                          n13 = Interactions_3[imex == 1 & GO_IPI == 1,.N],
                          n123 = Interactions_3[imex == 1 & BioGRID == 1 & GO_IPI == 1,.N], 
                          category = c("IMEx", 
                                       "BioGRID", 
                                       "GO_IPI"), 
                          lty = rep("blank", 3), 
                          fill = c("blue", "red", "green"), 
                          alpha = rep(0.5, 3), cat.pos = c(350, 25, 160), 
                          cat.dist = c(0.08,0.035, 0.08), 
                          cat.cex = c(1.4,1.4, 1.4), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 1.5
)
popViewport()
}

## ----GeneID_to_publication_IDs, results='hide'---------------------------
## getting the publication ID for each gene
{
geneID2pubmed_url = "ftp://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz"
geneID2pubmed_filename = paste0("./geneID2pubmed_release_", format(Sys.Date(), "%m-%Y.gz"))
geneID2pubmed_filename_txt = substr(geneID2pubmed_filename, 1, nchar(geneID2pubmed_filename)-3)
if(!file.exists(geneID2pubmed_filename_txt)) {
  downloader::download(geneID2pubmed_url, geneID2pubmed_filename)
  R.utils::gunzip(geneID2pubmed_filename)
  gitignore = c(readLines("../.gitignore"), paste0("/dsp_comparison/geneID2pubmed_release_", format(Sys.Date(), "%m-%Y")))
  write(gitignore,"../.gitignore")
}
}

# getting all pmids which have genes associated with them
geneID2pubmed = fread(geneID2pubmed_filename_txt, colClasses = c("character","character","character"))
pubmedIDs = geneID2pubmed[`#tax_id` == "9606",unique(PubMed_ID)]

## ------------------------------------------------------------------------
# enrichment of resources in imex-curated publications
# all gene-assiciated pmids
total_Npmids = length(pubmedIDs)
imex_Npmids = length(pubcomp_table_final[imex == 1,unique(pmid)])
# counting pmids in each of the resources
resources = colnames(pubcomp_table_final)[2:8]
resources_Npmids = sapply(resources, function(x) eval(parse(text=paste0("length(pubcomp_table_final[",x," == 1,unique(pmid)])"))))
# how many pmids overlap with IMEx
resources_Npmids_in_IMEx = sapply(resources, function(x) eval(parse(text=paste0("length(pubcomp_table_final[imex == 1 & ",x," == 1,unique(pmid)])"))))
resources_Npmids2 = rbind(resources_Npmids, resources_Npmids_in_IMEx)
# hypergeometric test
apply(resources_Npmids2,2, function(x){
                phyper(q = x[2], m = imex_Npmids, n = total_Npmids - imex_Npmids, k = x[1], lower.tail = F, log.p = FALSE)
})

## ------------------------------------------------------------------------
# all pmids in our resources
total_Npmids_our = length(pubcomp_table_final[,unique(pmid)])
# hypergeometric test
apply(resources_Npmids2,2, function(x){
                phyper(q = x[2], m = imex_Npmids, n = total_Npmids_our - imex_Npmids, k = x[1], lower.tail = F, log.p = FALSE)
})

## ------------------------------------------------------------------------
BioGRID_Npmids = length(pubcomp_table_final[BioGRID == 1,unique(pmid)])
resources_Npmids_in_BioGRID = sapply(resources, function(x) eval(parse(text=paste0("length(pubcomp_table_final[BioGRID == 1 & ",x," == 1,unique(pmid)])"))))
resources_Npmids2 = rbind(resources_Npmids, resources_Npmids_in_BioGRID)
apply(resources_Npmids2,2, function(x){
                phyper(q = x[2], m = BioGRID_Npmids, n = total_Npmids - BioGRID_Npmids, k = x[1], lower.tail = F, log.p = FALSE)
})

## ---- fig.width = 12, fig.height= 6--------------------------------------
# enrichment of resources in imex-curated publications
# all gene-assiciated pmids
total_NpmidsIMEx = length(pubmedIDs)
imex_NpmidsIMEx = length(pubcomp_table_final[imex == 1,unique(pmid)])
# counting pmids in each of the resources
resourcesIMEx = colnames(pubcomp_table_final)[2:8]
resources_NpmidsIMEx = sapply(resourcesIMEx, function(x) eval(parse(text=paste0("length(pubcomp_table_final[",x," == 1,unique(pmid)])"))))
# how many pmids overlap with IMEx
resources_Npmids_in_IMEx = sapply(resourcesIMEx, function(x) eval(parse(text=paste0("length(pubcomp_table_final[imex == 1 & ",x," == 1,unique(pmid)])"))))
resources_Npmids2IMEx = rbind(resources_NpmidsIMEx, resources_Npmids_in_IMEx)
resources_Npmids2IMEx = cbind(resources_Npmids2IMEx, gene2pubmed = c(total_NpmidsIMEx, imex_NpmidsIMEx))
resources_Npmids2IMEx = rbind(resources_Npmids2IMEx, fraction = c(resources_Npmids2IMEx[2,]/resources_Npmids2IMEx[1,]))

# enrichment of resources in BioGRID-curated publications
# all gene-assiciated pmids
total_NpmidsBioGRID = length(pubmedIDs)
BioGRID_NpmidsBioGRID = length(pubcomp_table_final[BioGRID == 1,unique(pmid)])
# counting pmids in each of the resources
resourcesBioGRID = colnames(pubcomp_table_final)[2:8]
resources_NpmidsBioGRID = sapply(resourcesBioGRID, function(x) eval(parse(text=paste0("length(pubcomp_table_final[",x," == 1,unique(pmid)])"))))
# how many pmids overlap with BioGRID
resources_Npmids_in_BioGRID = sapply(resourcesBioGRID, function(x) eval(parse(text=paste0("length(pubcomp_table_final[BioGRID == 1 & ",x," == 1,unique(pmid)])"))))
resources_Npmids2BioGRID = rbind(resources_NpmidsBioGRID, resources_Npmids_in_BioGRID)
resources_Npmids2BioGRID = cbind(resources_Npmids2BioGRID, gene2pubmed = c(total_NpmidsBioGRID, BioGRID_NpmidsBioGRID))
resources_Npmids2BioGRID = rbind(resources_Npmids2BioGRID, fraction = c(resources_Npmids2BioGRID[2,]/resources_Npmids2BioGRID[1,]))

resources_Npmids2IMEx = as.data.table(t(resources_Npmids2IMEx),keep.rownames = "resource")
resources_Npmids2IMEx[, imex_biogrid := "IMEx"]
colnames(resources_Npmids2IMEx)[2:3] = c("resources_Npmids", "resources_Npmids_in_resource")
resources_Npmids2BioGRID = as.data.table(t(resources_Npmids2BioGRID),keep.rownames = "resource")
resources_Npmids2BioGRID[, imex_biogrid := "BioGRID"]
colnames(resources_Npmids2BioGRID)[2:3] = c("resources_Npmids", "resources_Npmids_in_resource")

resources_Npmids2 = rbind(resources_Npmids2IMEx, resources_Npmids2BioGRID)
resources_Npmids2[, label_col := paste0(signif(as.numeric(resources_Npmids_in_resource),3)," /\n",signif(as.numeric(resources_Npmids),3))]

ggplot(data = resources_Npmids2, aes(x = resource, y = fraction, fill = imex_biogrid,label = label_col)) +
        scale_y_continuous(breaks = seq(0,1.2,0.05))+
        ylab("fraction of papers in IMEx, or BioGRID")+
        xlab(NULL)+
        theme(text = element_text(size = 17),
              axis.text.x = element_text(angle =-25, vjust = 0.6, size = 17), 
              #legend.position="none",
              plot.margin = unit(c(1,1,1,1), "cm"))+
        geom_bar(stat = "identity", width =0.9,  position = "dodge")+
        geom_text(vjust = -0.15, position = position_dodge(width = 0.9))

## ---- fig.width = 17, fig.height= 6--------------------------------------
# enrichment of resources in imex- or biogrid-curated interacting pairs
imex_NpmidsIMEx = length(comp_table_final[imex == 1,unique(pair_id)])
# counting pair_ids in each of the resources
resourcesIMEx = colnames(comp_table_final)[2:12]
resources_Npair_idsIMEx = sapply(resourcesIMEx, function(x) eval(parse(text=paste0("length(comp_table_final[",x," == 1,unique(pair_id)])"))))
# how many pair_ids overlap with IMEx
resources_Npair_ids_in_IMEx = sapply(resourcesIMEx, function(x) eval(parse(text=paste0("length(comp_table_final[imex == 1 & ",x," == 1,unique(pair_id)])"))))
resources_Npair_ids2IMEx = rbind(resources_Npair_idsIMEx, resources_Npair_ids_in_IMEx)
resources_Npair_ids2IMEx = rbind(resources_Npair_ids2IMEx, fraction = c(resources_Npair_ids2IMEx[2,]/resources_Npair_ids2IMEx[1,]))

# enrichment of resources in BioGRID-curated publications
BioGRID_Npair_idsBioGRID = length(comp_table_final[BioGRID == 1,unique(pair_id)])
# counting pair_ids in each of the resources
resourcesBioGRID = colnames(comp_table_final)[2:12]
resources_Npair_idsBioGRID = sapply(resourcesBioGRID, function(x) eval(parse(text=paste0("length(comp_table_final[",x," == 1,unique(pair_id)])"))))
# how many pair_ids overlap with BioGRID
resources_Npair_ids_in_BioGRID = sapply(resourcesBioGRID, function(x) eval(parse(text=paste0("length(comp_table_final[BioGRID == 1 & ",x," == 1,unique(pair_id)])"))))
resources_Npair_ids2BioGRID = rbind(resources_Npair_idsBioGRID, resources_Npair_ids_in_BioGRID)
resources_Npair_ids2BioGRID = rbind(resources_Npair_ids2BioGRID, fraction = c(resources_Npair_ids2BioGRID[2,]/resources_Npair_ids2BioGRID[1,]))

resources_Npair_ids2IMEx = as.data.table(t(resources_Npair_ids2IMEx),keep.rownames = "resource")
resources_Npair_ids2IMEx[, imex_biogrid := "IMEx"]
colnames(resources_Npair_ids2IMEx)[2:3] = c("resources_Npair_ids", "resources_Npair_ids_in_resource")
resources_Npair_ids2BioGRID = as.data.table(t(resources_Npair_ids2BioGRID),keep.rownames = "resource")
resources_Npair_ids2BioGRID[, imex_biogrid := "BioGRID"]
colnames(resources_Npair_ids2BioGRID)[2:3] = c("resources_Npair_ids", "resources_Npair_ids_in_resource")

resources_Npair_ids2 = rbind(resources_Npair_ids2IMEx, resources_Npair_ids2BioGRID)
resources_Npair_ids2[, label_col := paste0(signif(as.numeric(resources_Npair_ids_in_resource),3)," /\n",signif(as.numeric(resources_Npair_ids),3))]

ggplot(data = resources_Npair_ids2, aes(x = resource, y = fraction, fill = imex_biogrid,label = label_col)) +
        scale_y_continuous(breaks = seq(0,1.2,0.05))+
        ylab("fraction of interacting pairs in IMEx or BioGRID")+
        xlab(NULL)+
        theme(text = element_text(size = 17),
              axis.text.x = element_text(angle =-25, vjust = 0.6, size = 17), 
              #legend.position="none",
              plot.margin = unit(c(1,1,1,1), "cm"))+
        geom_bar(stat = "identity", width =0.9,  position = "dodge")+
        geom_text(vjust = -0.15, position = position_dodge(width = 0.9))

## ---- fig.width = 12, fig.height= 6--------------------------------------
Interactions_ids = Interactions_1[imex == 1 | BioGRID == 1 | GO_IPI == 1, unique(pair_id)]
Pathways_ids = Pathways_1[reactome == 1 | OmniPath == 1 | STRING_pathway_inference == 1, unique(pair_id)]
textmining_ids = textmining_1[tm_epmc == 1 | EVEX == 1 | STRING_textmining == 1, unique(pair_id)]
Predictions_ids = Predictions_1[iid_pred == 1 | STRING_phylo_predictions == 1, unique(pair_id)]


area1 = length(Interactions_ids)
area2 = length(Pathways_ids)
area3 = length(textmining_ids)
area4 = length(Predictions_ids)
n12 = sum(!is.na(match(Interactions_ids, Pathways_ids)))
n13 = sum(!is.na(match(Interactions_ids, textmining_ids)))
n14 = sum(!is.na(match(Interactions_ids, Predictions_ids)))
n23 = sum(!is.na(match(Pathways_ids, textmining_ids)))
n24 = sum(!is.na(match(Pathways_ids, Predictions_ids)))
n34 = sum(!is.na(match(textmining_ids, Predictions_ids)))
n123 = length(intersect(intersect((Interactions_ids),(Pathways_ids)),(textmining_ids)))
n124 = length(intersect(intersect((Interactions_ids),(Pathways_ids)),(Predictions_ids)))
n134 = length(intersect(intersect((Interactions_ids),(textmining_ids)),(Predictions_ids)))
n234 = length(intersect(intersect((Pathways_ids),(textmining_ids)),(Predictions_ids)))
n1234 = length(intersect(intersect(intersect((Interactions_ids),(Pathways_ids)),(textmining_ids)), (Predictions_ids)))


venn_Interactions_1 = draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24,
    n34, n123, n124, n134, n234, n1234, 
                          category = c("Interactions", 
                                       "Pathways", 
                                       "Textmining",
                                       "Predictions"), 
                          lty = rep("blank", 4), 
                          fill = c("blue", "red", "green", "grey"), 
                          alpha = rep(0.5, 4), cat.pos = c(350, 25, 160,45), 
                          cat.dist = c(0.08,0.08, -0.08, 0.08), 
                          cat.cex = c(1.4,1.4, 1.4, 1.4), scaled = TRUE, euler.d = TRUE,  margin = 0.05,
                          print.mode =  'raw',
                          cex = 1.5
)

## ---- fig.width = 12, fig.height= 6--------------------------------------
combs_groups = c("Pathways","Textmining","Predictions", "Pathways & Textmining", "Pathways & Predictions", "Textmining & Predictions", "Pathways & Textmining & Predictions")
comp_final_dt <- data.table(group = combs_groups, 
                            overlap = c(n12/area2, n13/area3, n14/area4, n123/n23, n124/n24, n134/n34, n1234/n234),
                            label = c(paste0(signif(n12,3)," /\n",signif(area2,3)), 
                                      paste0(signif(n13,3)," /\n",signif(area3,3)), 
                                      paste0(signif(n14,3)," /\n",signif(area4,3)),
                                      paste0(signif(n123,3)," /\n",signif(n23,3)),
                                      paste0(signif(n124,3)," /\n",signif(n24,3)),
                                      paste0(signif(n134,3)," /\n",signif(n34,3)),
                                      paste0(signif(n1234,3)," /\n",signif(n234,3))),
                            sort = c(3,1,2,6,4,5,7),
                            col = c("slategray1","slategray1","slategray1","skyblue2","skyblue2","skyblue2","skyblue4"))


g <- ggplot(comp_final_dt,aes(x=reorder(group,sort),y=overlap,fill=col))
g <- g + geom_bar(stat="identity")
g <- g + guides(fill=FALSE)
g <- g + ylab("fraction of interacting pairs in IMEx")
g <- g + scale_fill_manual(values=c("skyblue2","skyblue4","slategray1"))
g <- g + geom_text(label = comp_final_dt$label, vjust = -0.5)+
        geom_text(label = signif(c(area1/area1, n12/area2, n13/area3, n14/area4, n123/n23, n124/n24, n134/n34, n1234/n234),3), vjust = 1.4)
g <- g + theme(plot.title = element_text(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.text.y = element_text(size=12),
                 axis.text.x = element_text(size=12,angle =-45, vjust = 0.6),
                 axis.title.y = element_text(size=12),
               axis.title.x = element_blank())
g

