# darkspaceproject
Trying to find out how much is out there that needs curation...

### Description and goals

This project aims to try to find out how much published molecular interaction data is out there that has not been curated by molecular interaction databases, particularly IMEx-complying ones. Several different strategies will be used for this purpose and an integrated view of the estimation will be produced in the form of a table. 

### Strategies list

This a list of the different approaches we will take, in no particular order, with a brief comment. Some of these will not provide with potential interactiosn to curate, but just help ranking lists generated with other approaches. 

* Reactome inferred pairs: We already have an R workflow that builts a big file with over 600000 potential interactions inferred from Reactome, along with referenced PMIDs. The file probably needs cleaning an updating and the workflow is not fully automated. 
* FpClass: Interaction prediction tool from Igor Jurisica's lab, does not provide potential PMIDs.
  * We have a sub-project using FpClass to try and identify potential curation errors.
  * There is an FpClass-generated human interactome from May 2015, needs updating. 
  * A new predictive database has been put together by Igor's group, IID (http://dcv.uhnres.utoronto.ca/iid/), but no full download is available for now. 
  * Igor's lab has been contacted to try to get a link to full, updated versions of FpClass-predicted interactomes or IID predicted part. 
* Text-mining strategies: Two different approaches here.
  * EMPC (Senay Kafkas) approach: Produces as output a list of UniProtKB accessions for pairs of genes/proteins co-occurring in the same sentence where a particular term (selected from a pre-made list) was found. The PMIDs where they were found are indicated, but no score measuring reliability, only the number of PMIDs.
    * The search is done through the full text when the paper is open access (PMC) or only abstract otherwise. 
  * Laitor/PESCADOR/MedLine Ranker: Similar to the approach above, but uses text mining tools that allow more sophisticated ranking of the term co-occurrences. 
* GO IPIs: Now they are avaialbable through a new PSICQUIC server (EBI-GOA-nonIntAct).
* Denes interaction map.
* STRING database predictions.
* UniProtKB CC subunit lines: They detail complexes internal interactions that we might not have. 
* Structure predictions using LAMs (Local Approximate Models)
  * Ask Tamas Korcsmaros where to find this kind of information.
* Domain Interaction Prediction: different databases can be used to get this information.
  * DOMINE (http://domine.utdallas.edu/cgi-bin/Domine)
  * 3did (http://3did.irbbarcelona.org)
  * Others (ask Tamas)
* Genetic interactions: Although the overlap between protein and genetic interactions is negligible, it might help identifying those predicted interactions that have higher biological interest. 
* BioGrid data: We need to identify how much information is curated in BioGrid as well, since any predicted interaction will be prioritized down if it is already curated there. 
