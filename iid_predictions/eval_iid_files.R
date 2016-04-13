if(!file.exists("iid.human.all.txt.gz")){
        download.file("http://dcv.uhnres.utoronto.ca/downloads/iid.human.all.txt.gz", destfile="./source_files/iid.human.all.txt.gz", method= "curl")
        download_date = date()
}

iidhumanall <- read.table(gzfile("./source_files/iid.human.all.txt.gz"), header = T, sep = "\t")
closeAllConnections()

table(iidhumanall$evidence_type)

# I extract the predicted interactions only (those with no database). 

iid_nodb <- iidhumanall[iidhumanall$dbs=="",]
table(iid_nodb$evidence_type)

library(dplyr)
iid_nodb_sel <- unique(select(iid_nodb,uniprot1,uniprot2))

# I generate the pair_ids

iid_nodb_sel$pair_id <- apply(iid_nodb_sel[,1:2], 1,function(i){
        paste(sort(i),collapse = "_")
})
iid_nodb_sel$iid <- "yes"

# Now I load the IMEx dataset.

system("Rscript ./scripts/IMEx_ds_generator.R --save")

imex_full <- read.delim("./results/imex_full.txt", header=T, sep="\t",colClasses="character")

# I compare both sets

imex_sel <- unique(select(imex_full,pair_id_clean,pair_id_clean,id_a_clean,id_b_clean,taxid_a,taxid_b,pubid))
imex_sel$imex <- "yes"
imex_pairs <- unique(select(imex_sel, pair_id=pair_id_clean,imex))

comp <- unique(merge(iid_nodb_sel,imex_pairs,by="pair_id",all=T))

comp <- mutate(comp, db_pair =
                       ifelse(iid == "yes" & is.na(imex), "iid",
                              ifelse(is.na(iid) & imex == "yes", "imex",
                                     ifelse(iid == "yes" & imex == "yes", "iid & imex",
                                            "check"))))

comp$db_pair <- as.factor(comp$db_pair)

table(comp$db_pair,useNA="ifany")

comp_simple <- unique(select(comp, pair_id,db_pair))

write.table(comp_simple,"./results/pairs_iid_vs_imex.txt",col.names=T,row.names=F,quote=F,sep="\t")
