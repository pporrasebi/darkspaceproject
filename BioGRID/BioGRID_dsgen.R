## ----echo=FALSE, cache=TRUE----------------------------------------------
intact_date <- date()

## ----eval=FALSE----------------------------------------------------------
if(!file.exists("./source_files/intact_mitab27.txt")){
  download.file("ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt", destfile = "./source_files/intact_mitab27.txt")
}

## ----eval=FALSE,message=FALSE, warning=FALSE-----------------------------
system("perl ./scripts/MITAB27extractor_v13.pl ./source_files/intact_mitab27.txt ./processed_files/intact_pairs.txt")

## ------------------------------------------------------------------------
intact_full <- read.delim("./processed_files/intact_pairs.txt", header = T, sep = "\t", colClasses = "character")

## ----eval=TRUE,echo=FALSE, cache=TRUE------------------------------------
dip_date <- date()

## ----eval=FALSE,message=FALSE, warning=FALSE-----------------------------
system("perl ./scripts/MITAB25extractor_v12.pl ./source_files/dip20160731.txt ./processed_files/dip20160731_pairs.txt")

## ------------------------------------------------------------------------
dip_full <- read.delim("./processed_files/dip20160731_pairs.txt", header = T, sep = "\t", colClasses = "character")

## ----eval=T--------------------------------------------------------------
intact_sel <- unique(subset(intact_full, select = c("pair_id", "id_a","id_b","pair_id_clean","id_a_clean","id_b_clean","taxid_a","taxid_b","pubid")))

## ------------------------------------------------------------------------
imex_full <- unique(rbind(intact_sel,dip_full))

write.table(imex_full, "./results/imex_full.txt", quote=F, sep ="\t", row.names = F, col.names = T)

## ------------------------------------------------------------------------
imex_pmids <- data.frame(unique(imex_full$pubid))

write.table(imex_pmids, "./results/imex_pmids.txt", quote=F, sep ="\t", row.names = F, col.names = T)

