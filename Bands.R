#input = character vector HUGO gene symbols
library(org.Hs.eg.db)
Bands <- function(input){
  output <- as.data.frame(input)
  # Annotate Entrez Gene ID
  mapped_GeneID <- as.data.frame(org.Hs.egSYMBOL2EG[mappedkeys(org.Hs.egSYMBOL2EG)])
  for (i in 1:length(input)){
    index = which(input[i] == mapped_GeneID$symbol)
    if (length(index) > 0)
    {output[i, 2] <- mapped_GeneID$gene_id[index[1]]}
    else {output[i,2] <- NA}
  }
  colnames(output)[2] <- "Gene ID"
  # Annotate Band
  mapped_CHRBAND <- as.data.frame(org.Hs.egMAP[mappedkeys(org.Hs.egMAP)])
  for (i in 1:length(input)){
    index = which(output$`Gene ID`[i] == mapped_CHRBAND$gene_id)
    if (length(index) > 0)
    {output[i,3] <- mapped_CHRBAND$cytogenetic_location[index[1]]}
    else {output[i,3] <- NA}
  }
  colnames(output)[3] <- "Band"
  return(output)
}