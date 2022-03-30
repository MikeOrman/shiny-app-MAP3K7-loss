#Formats a manhattan plot object (variable class = GRanges) 
#that can be plotted using the karyoploteR function
#input = dataframe, where column 1 = Hugo Names, and column 2 = p-values (or other quantification you wish to map)
library(org.Hs.eg.db)
library(GenomicRanges)
library(karyoploteR)
library(dplyr)
library(tidyr)
manhattan.object <- function(input){
  main <- input
  #Annotate Entrez Gene ID
  #Obtain dataframe of Gene ID mapped to HGNC Symbol
  mapped_HGNC <- as.data.frame(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)])
  #Map Entrez Gene ID to HGNC Symbols in main dataframe
  for (i in 1:nrow(main)){
    index = as.numeric(which(main[i,1] == mapped_HGNC$symbol))
    if (length(index) >= 1)
    {main[i,3] <- mapped_HGNC$gene_id[index][1]}
    else {main[i,3] <- NA}
  }
  colnames(main) <- c("Hugo_Symbol", "pval", "Gene ID")
  main <- na.omit(main)
  #Annotate Chromosome Number
  mapped_CHR <- as.data.frame(org.Hs.egCHR[mappedkeys(org.Hs.egCHR)])
  #Map Entrez Gene ID to Chromosome Number in main dataframe
  for (i in 1:nrow(main)){
    index = (which(main$`Gene ID`[i] == mapped_CHR$gene_id))
    if (length(index) >= 1)
    {main[i,4] <- as.character(mapped_CHR$chromosome[index][1])}
    else {main[i,4] <- NA}
  }
  colnames(main)[4] <- "Chr"
  main <- na.omit(main)
  #Annotate Starting Position
  mapped_CHRLOC <- as.data.frame(org.Hs.egCHRLOC[mappedkeys(org.Hs.egCHRLOC)])
  for (i in 1:nrow(main)){
    index = as.numeric(which(main$`Gene ID`[i] == mapped_CHRLOC$gene_id))
    if (length(index) >= 1)
    {main[i,5] <- abs(mapped_CHRLOC$start_location[index][1])}
    else {main[i,5] <- NA}
  }
  colnames(main)[5] <- "Starting Position"
  main <- na.omit(main)
  #Annotate Ending Position
  mapped_CHRLOCEND <- as.data.frame(org.Hs.egCHRLOCEND[mappedkeys(org.Hs.egCHRLOCEND)])
  for (i in 1:nrow(main)){
    index = as.numeric(which(main$`Gene ID`[i] == mapped_CHRLOCEND$gene_id))
    if (length(index) >= 1)
    {main[i,6] <- abs(mapped_CHRLOCEND$end_location[index][1])}
    else {main[i,6] <- NA}
  }
  colnames(main)[6] <- "Ending Position"
  main <- na.omit(main)
  #Annotate Strand
  for (i in 1:nrow(main)){
    index = as.numeric(which(main$`Gene ID`[i] == mapped_CHRLOC$gene_id))
    if (mapped_CHRLOC$start_location[index] > 0)
    {main[i,7] <- "+"}
    else {main[i,7] <- "-"}
  }
  colnames(main)[7] <- "Strand"
  main <- na.omit(main)
  #Create GR object
  #Order data frame
  main1 <- main[order(main$Chr),]
  #Add "chr" string to chromosome number
  #Map Entrez Gene ID to Chromosome Number in main dataframe
  for (i in 1:nrow(main1)){
    main1[i,4] <- paste("chr",as.character(main1[i,4]), sep = "")
  }
  #Make RLE object holding each gene's chromosome number
  Chr.table <- main1 %>% dplyr::count(Chr)
  Chr <- Chr.table[["Chr"]]
  Chr.counts <- Chr.table[["n"]]
  Chr.RLE <- Rle(Chr, Chr.counts)
  #Make IRange object holding each gene's base pair range on its' respective chromosome
  BP.IRange <- IRanges(start = main1$`Starting Position`, end = main1$`Ending Position`)
  #Create GR object for plotting
  gr <- GRanges(seqnames = Chr.RLE, ranges = BP.IRange)
  names(gr) <- main1$Hugo_Symbol
  gr$pval <- main1$pval
  return(gr)
}