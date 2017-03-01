IMEx dataset generator
========================================================

#### Load IntAct data



I download the latest version of the IntAct data. The file was downloaded on Mon Apr 11 16:50:39 2016. 


```r
if(!file.exists("./source_files/intact_mitab27.txt")){
  download.file("ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt", destfile = "./source_files/intact_mitab27.txt")
}
```

I need to clean up the intact file.


```r
system("perl ./scripts/MITAB27extractor_v13.pl ./source_files/intact_mitab27.txt ./processed_files/intact_pairs.txt")
```

```r
intact_full <- read.delim("./processed_files/intact_pairs.txt", header = T, sep = "\t", colClasses = "character")
```

#### Load data from DIP


I incorporate the data from DIP as well, taking the latest available release there is. Data was downloaded on Mon Apr 11 16:52:07 2016 and needed to be downloaded manually, since the DIP download options only appear for logged-in users.  

Then I clean up the data as I did with the IntAct dataset. 


```r
system("perl ./scripts/MITAB25extractor_v12.pl ./source_files/dip20160114.txt ./processed_files/dip20160114_pairs.txt")
```

```r
dip_full <- read.delim("./processed_files/dip20160114_pairs.txt", header = T, sep = "\t", colClasses = "character")
```

#### Merge IntAct and DIP data to create an IMEx dataset

The original DIP data is in MITAB 2.5 format and IntAct is 2.7, so I select only 2.5-compatible files from the processed IntAct dataset. 


```r
intact_sel <- unique(subset(intact_full, select = c("pair_id", "id_a","id_b","pair_id_clean","id_a_clean","id_b_clean","taxid_a","taxid_b","pubid")))
```

Now I create a full IMEx dataset. 


```r
imex_full <- unique(rbind(intact_sel,dip_full))

write.table(imex_full, "./results/imex_full.txt", quote=F, sep ="\t", row.names = F, col.names = T)
```

The full dataset contains 534683 interacting pairs. I create a list of PMIDs that have been curated by IMEx-complying databases (IMPORTANT: not all these publications are curated to IMEx standards). 


```r
imex_pmids <- data.frame(unique(imex_full$pubid))

write.table(imex_pmids, "./results/imex_pmids.txt", quote=F, sep ="\t", row.names = F, col.names = T)
```

There are 20168 publications in the full dataset of publications curated in IMEx-complying databases. 
