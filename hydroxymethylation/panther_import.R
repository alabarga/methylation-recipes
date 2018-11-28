library(stringr)
library(plyr)

read.panther <- function(filename="PTHR8.0_chicken") {
  ## Panther chicken classification file from
  ## ftp://chelsea.usc.edu/sequence_classifications/current_release/
  
  #Pathway accession
  #Pathway name
  #Pathway component accession
  #Pathway component name
  #UniProt ID
  #Protein definition
  #Confidence code
  #Evidence
  #Evidence type (e.g., PubMed, OMIM)
  #PANTHER subfamily ID (the SF to which the sequence belongs to)
  #PANTHER subfamily name
  
  panther <- read.delim(filename, sep="\t", head=F,
                        stringsAsFactors=F)

  colnames(panther) <- c("pathway.acc","pathway.name", "pathway.comp.acc",
                            "pathway.comp.name","uniprot","protein.definition","confidence.code","evidence.id",'evidence.type',"panther.subfamily","panther.subfamily.name")

  ## accession numbers to Entrez, UniProt, Ensembl
  panther$ensembl.id <- str_match(panther$uniprot, "Ensembl=(.*)\\|")[,2]
  panther$uniprot.id <- str_match(panther$uniprot, "UniProtKB=(.*)$")[,2]
  panther$organism <- str_match(panther$uniprot, "^(\\w*)\\|")[,2]

  return(panther)

}

## Extract pathway information from panther classification
get.pathways <- function(panther, acc) {
  
  ix.uni <- which(panther$data$uniprot.id == acc)
  ix.ens <- which(panther$data$ensembl.id == acc)
  ix.ent <- which(panther$data$entrez.id == acc)
  
  pathways <- character(0)
  
  if (length(ix.uni) > 0) {
    ix <- ix.uni
  } else if (length(ix.ens > 0)) {
    ix <- ix.ens
  } else if (length(ix.ent > 0)){
    ix <- ix.ent
  } 
  panther$panther.pathway[[ix]]
  
}
