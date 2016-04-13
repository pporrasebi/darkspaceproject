# darkspaceproject
Trying to find out how much is out there that needs curation...

### Description and goals

This project aims to try to find out how much published molecular interaction data is out there that has not been curated by molecular interaction databases, particularly IMEx-complying ones. Several different strategies will be used for this purpose and an integrated view of the estimation will be produced in the form of a table. 

### Strategies list

This a list of the different approaches considered in order to get an estimation. Some of these will not provide with potential interactions to curate, but they can be used to rank lists generated with other approaches. 

#### Strategies we are currently working on

* Reactome inferred pairs: Protein pairs are inferred by their association to a reaction as input, output or catalyst. The dataset comprises over 320000 potential pairs, with associated PMIDs for each one of them. 
* IID-predicted: The Integrated Interactions Database (IID) is a predictive meta-database built in Igor Jurisica's lab (http://dcv.uhnres.utoronto.ca/iid/).
  * It comprises data obtained from primary databases (IMEx consortium, BioGRID, HPRD...) plus computational predictions. 
  * We use the subset of data that has been produced by computational prediction and orthologous similarity.
  * The predictive subset comprises over 700000 potential pairs, no associated PMIDs. 
* Text-mining: Collaboration with Senay Kafkas (EPMC). Her pipeline produces as output a list of UniProtKB accessions for pairs of genes/proteins co-occurring in the same sentence where a particular term (selected from a pre-made list) was found. The PMIDs where they were found are indicated, but no score measuring reliability, only the number of PMIDs.
    * The search is done through the full text when the paper is open access (PMC) or only abstract otherwise. 
    * The dataset comprises almost 740000 potential pairs, along with the PMIDs for the papers from which they were inferred. 

#### Strategies we are considering, but have not been implemented yet

* BioGrid data: We need to identify how much information is curated in BioGrid as well, since any predicted interaction will be prioritized down if it is already curated there. 
* Laitor/PESCADOR/MedLine Ranker: (Miguel Andrade and Adriano Barbosa) Text mining tools that allow ranking of the term co-occurrences.
* GO IPIs: 'Inferred from Protein Interaction' annotations made by GO curators. Available through a new PSICQUIC server (EBI-GOA-nonIntAct).
* Denes interaction map.
* STRING database predictions.
* UniProtKB CC subunit lines: They detail complexes internal interactions that might not be curated by IMEx. 
* Structure predictions using LAMs (Local Approximate Models)
* Domain Interaction Prediction: different databases can be used to get this information.
  * DOMINE (http://domine.utdallas.edu/cgi-bin/Domine)
  * 3did (http://3did.irbbarcelona.org)
* Genetic interactions: Although the overlap between protein and genetic interactions is negligible, it might help identifying those predicted interactions that have higher biological interest. 
