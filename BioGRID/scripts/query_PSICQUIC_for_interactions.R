## the function to query PSICQUIC for interactions given:
### SPECIES_NAME
### SPECIES_ID (default NA), can be "all" to download all species
### databases (Default -   databases <- c("IntAct", "MINT", "bhf-ucl", "MPIDB", "MatrixDB", "HPIDb","I2D-IMEx","InnateDB-IMEx", "MolCon", "UniProt", "MBInfo")), can be "IMEx" to find all IMEX databases using MI:0959 term
### date - option allows to read already saved files
### detmethod - optional
### pmethod - optional

query_PSICQUIC_for_interactions = function(SPECIES_ID = NA, SPECIES_NAME = NA, databases = NA, date = Sys.Date(), detmethod = NA, pmethod = NA, database.name = NA, return_data = T, show_summary = T, MITAB = "tab27") {
  
  ## checks if databases argument was provided, if not - sets default - all IMEX databases
  if(is.na(databases)[1]){ 
    databases <- c("IntAct", "MINT", "bhf-ucl", "MPIDB", "MatrixDB", 
                   "HPIDb","I2D-IMEx","InnateDB-IMEx", "MolCon", "UniProt", "MBInfo")
  }
  ## checks if databases argument is "IMEx", if yes - sets all IMEX databases using MI:0959 term
  if(databases == "IMEx"){
    databases_url = "http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=STATUS&tags=%22MI:0959%22&format=txt"
    databases_temp = readLines(databases_url)
    databases = gsub("=.+","", databases_temp)
  }
  
  ## Converts SPECIES_NAME to SPECIES_ID if SPECIES_ID is not stated
  if(is.na(SPECIES_ID) & !is.na(SPECIES_NAME)) { 
    source("SPECIES_NAME_TO_ID.R")
    SPECIES_ID = SPECIES_NAME_TO_ID(SPECIES_NAME)$SPECIES_ID }
  else {SPECIES_ID = SPECIES_ID}
  
  ##============================================================================##
  ## Constructs database filename and PSICQUIC_query depending on 
  ## whether detmethod and pmethod were provided as arguments
        # if database.name is NA - generate database filename
if(is.na(database.name)){
        ##============================================================================##
if(SPECIES_ID == "all"){
        if(!is.na(detmethod)){
                if(!is.na(pmethod)){
                        database.name <- paste("./Data/","databaseName_", databases[1], "...", databases[length(databases)],"_","speciesID_",SPECIES_ID,"_",SPECIES_NAME,"_detmethod_",detmethod, "_pmethod_", pmethod, "_", date, sep = "")
                        PSICQUIC_query = paste("detmethod:",detmethod," AND ","pmethod:",pmethod, sep = "")
                }
                if(is.na(pmethod)){
                        database.name <- paste("./Data/","databaseName_", databases[1], "...", databases[length(databases)],"_","speciesID_",SPECIES_ID,"_",SPECIES_NAME,"_detmethod_",detmethod, "_", date, sep = "")
                        PSICQUIC_query = paste("detmethod:",detmethod, sep = "")
                }
        }
        if(is.na(detmethod)){
                database.name <- paste("./Data/","databaseName_", databases[1], "...", databases[length(databases)],"_","speciesID_",SPECIES_ID,"_",SPECIES_NAME,"_",date, sep = "")
                PSICQUIC_query = paste("*", sep = "")
        }  
}
        ##============================================================================## 
if(SPECIES_ID != "all"){
  if(!is.na(detmethod)){
    if(!is.na(pmethod)){
      database.name <- paste("./Data/","databaseName_", databases[1], "...", databases[length(databases)],"_","speciesID_",SPECIES_ID,"_",SPECIES_NAME,"_detmethod_",detmethod, "_pmethod_", pmethod, "_", date, sep = "")
      PSICQUIC_query = paste("species:",SPECIES_ID," AND ","detmethod:",detmethod," AND ","pmethod:",pmethod, sep = "")
    }
    if(is.na(pmethod)){
      database.name <- paste("./Data/","databaseName_", databases[1], "...", databases[length(databases)],"_","speciesID_",SPECIES_ID,"_",SPECIES_NAME,"_detmethod_",detmethod, "_", date, sep = "")
      PSICQUIC_query = paste("species:",SPECIES_ID," AND ","detmethod:",detmethod, sep = "")
    }
  }
  if(is.na(detmethod)){
    database.name <- paste("./Data/","databaseName_", databases[1], "...", databases[length(databases)],"_","speciesID_",SPECIES_ID,"_",SPECIES_NAME,"_",date, sep = "")
    PSICQUIC_query = paste("species:",SPECIES_ID, sep = "")
  }
}
}
        ##============================================================================##
        # if database.name is not NA - use user-specified database filename
if(!is.na(database.name)){
        ##============================================================================##
if(SPECIES_ID == "all"){
        if(!is.na(detmethod)){
                if(!is.na(pmethod)){
                        database.name <- database.name
                        PSICQUIC_query = paste("detmethod:",detmethod," AND ","pmethod:",pmethod, sep = "")
                }
                if(is.na(pmethod)){
                        database.name <- database.name
                        PSICQUIC_query = paste("detmethod:",detmethod, sep = "")
                }
                }
                if(is.na(detmethod)){
                        database.name <- database.name
                        PSICQUIC_query = paste("*", sep = "")
                }
        }
        ##============================================================================##
if(SPECIES_ID != "all"){
        if(!is.na(detmethod)){
                if(!is.na(pmethod)){
                        database.name <- database.name
                        PSICQUIC_query = paste("species:",SPECIES_ID," AND ","detmethod:",detmethod," AND ","pmethod:",pmethod, sep = "")
                        }
                if(is.na(pmethod)){
                        database.name <- database.name
                        PSICQUIC_query = paste("species:",SPECIES_ID," AND ","detmethod:",detmethod, sep = "")
                        }
                }
                if(is.na(detmethod)){
                        database.name <- database.name
                        PSICQUIC_query = paste("species:",SPECIES_ID, sep = "")
                }
        }
        ##============================================================================##
}
        
  ## Checks if databases have been queried today, if not - sends query to the database
  if(!file.exists(database.name)) {
    print("dowloaded using PSICQUIC")
    ## Load PSICQUIC functionality
    suppressPackageStartupMessages(require(PSICQUIC))
    psicquic <- PSICQUIC()
    providers <- providers(psicquic)
    
    # query databases for all known SPECIES_ID protein interactions
    suppressPackageStartupMessages(require(data.table))
    SPECIES_ID_interactome = data.table()
    NO_SPECIES_ID_interactome = character(length = length(databases))
    for(indices in 1:length(databases)) {
      if(databases[indices] %in% providers) {
        ## Query for the number of interactions
        PSICQUIC_query1 = paste0(PSICQUIC_query, "?format=count") 
        N_interactions <- unlist(rawQuery(psicquic, databases[indices], PSICQUIC_query1))
        ## if there are any interactions - Query for the interactions by 1000 at a time
        if(N_interactions > 0){
          N_start = 1
          N_nrows = 2500
          for(n_starts in seq(from = N_start, to = N_interactions, by = N_nrows)){
            PSICQUIC_query2 = paste0(PSICQUIC_query, "?format=",MITAB,"&firstResult=", n_starts,"&maxResults=", N_nrows) 
            SPECIES_ID_interactome_d <- as.data.table(rawQuery(psicquic, databases[indices], PSICQUIC_query2))
            SPECIES_ID_interactome <- rbind(SPECIES_ID_interactome, SPECIES_ID_interactome_d)
          }
        }
        ## if the query finds no interactors
        else {
          NO_SPECIES_ID_interactome[indices] = databases[indices]
        }
      }
      ## if the database is not active
      else {
        NO_SPECIES_ID_interactome[indices] = databases[indices]
      }
    }
    ##============================================================================##
    ## Save dowloaded query result into file 
    fwrite(x = SPECIES_ID_interactome, file = database.name, sep = "\t")
    ##============================================================================##
    if(show_summary == T){
    ## Show what's found
    print(paste0("total number of interactions for",SPECIES_NAME,", detmethod(",detmethod,"), pmethod(",pmethod,"): ", nrow(SPECIES_ID_interactome)), quote = F)
    print(paste("there is no interactions in the databases: ", sep = ""), quote = F)
    print(NO_SPECIES_ID_interactome, quote = F)
    print(paste("the number of interactions per database ", sep = ""), quote = F)
    dbs = as.data.frame(table(SPECIES_ID_interactome$V13, useNA = "ifany"))
    colnames(dbs) = c("database", "N of interactions")
    print(dbs, quote = F)
    ##============================================================================##
    query_log_filename = paste("./Data/logs/","there is no interactions for ",SPECIES_NAME," in the databases ",date, sep = "", ".txt")
    write.table(NO_SPECIES_ID_interactome, query_log_filename, col.names=T,row.names=F,sep="\t",quote=F)
    ##============================================================================##
    }
    if(return_data == T){
    return(SPECIES_ID_interactome)
    }
    if(return_data == F){
            return("done")
    }
  }
  ##============================================================================##
  ## If file exists  - load, show what's found and return
  if(file.exists(database.name)) {
    print("loaded from file")
    SPECIES_ID_interactome = fread(database.name)
    ##============================================================================##  
    if(show_summary == T){
    query_log_filename = paste("./Data/logs/","there is no interactions for ",SPECIES_NAME," in the databases ",date, sep = "", ".txt")
    NO_SPECIES_ID_interactome = read.table(file = query_log_filename, header =T,sep="\t", stringsAsFactors = F)
    
    ## Show what's found
    print(paste("interactions for ",SPECIES_NAME,", detmethod(",detmethod,"), pmethod(",pmethod,"): ", sep = ""), quote = F)
    print(paste0("total number: ", nrow(SPECIES_ID_interactome)), quote = F)
    print(paste("there is no interactions in the databases: ", sep = ""), quote = F)
    print(NO_SPECIES_ID_interactome, quote = F)
    print(paste("the number of interactions per database ", sep = ""), quote = F)
    dbs = as.data.frame(table(SPECIES_ID_interactome$V13, useNA = "ifany"))
    colnames(dbs) = c("database", "N of interactions")
    print(dbs, quote = F)
    }
    ##============================================================================##
    if(return_data == T){
            return(SPECIES_ID_interactome)
    }
    if(return_data == F){
            return("done")
    }
  }
}