## ----set-options, echo=FALSE---------------------------------------------
options(width = 80)

## ---- tm_upload, cache=TRUE----------------------------------------------
require(data.table)

setwd("./source_files/")
system("tar xf PPI_tm_Nov16.txt.tar.gz")
# I need to clean up several undesired line breaks in the file
system("perl -i -pe 's/Mauviel et\n/Mauviel et/' PPI_tm_Nov16.txt")
system("perl -i -pe 's/According to Bia et\n/According to Bia et/' PPI_tm_Nov16.txt")
system("perl -i -pe 's/Lichty et\n/Lichty et/' PPI_tm_Nov16.txt")
tm_sel <- fread("PPI_tm_Nov16.txt", sep = "\t", header = F, colClasses=c("character","character","character","character","NULL","character","character","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL"),data.table = F)
system("rm PPI_tm_Nov16.txt")
setwd("../")

colnames(tm_sel) <- c("pmcid","pmid","date","upacs","location","symbols")

## ----format_data---------------------------------------------------------
tm_sel$upac_1 <- gsub("\\|\\|.*","",tm_sel$upacs)
tm_sel$upac_2 <- gsub(".*\\|\\|","",tm_sel$upacs)

library(splitstackshape)
tm_long_pt1 <- cSplit(tm_sel, "upac_1", sep = ",", direction = "long")
tm_long <- cSplit(tm_long_pt1, "upac_2", sep = ",", direction = "long")

tm_long$upac_1 <- as.character(tm_long$upac_1)
tm_long$upac_2 <- as.character(tm_long$upac_2)

## ----cache=TRUE----------------------------------------------------------
library(dplyr)
tm_long$pair_id <- apply(select(tm_long,upac_1,upac_2), 1,function(i){
  paste(sort(i),collapse = "_")
})
tm_long$tm <- 1

## ------------------------------------------------------------------------
tm_lite <- unique(select(tm_long,pair_id,pmid,tm))
write.table(tm_lite,"./results/pairs_pmids_tm.txt",col.names=T,row.names=F,quote=F,sep="\t")
system("gzip ./results/pairs_pmids_tm.txt")

