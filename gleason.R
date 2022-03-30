# Computes enrichment of high-risk gleason scores in primary subype patients stratified by gene expression
# of the gene in the summary table
# Input = input table from Integrated Analysis: Genomic + Transcriptomic script
clinical.sample <- read.table("TCGA clinical sample.txt", header = TRUE, sep = "\t", check.names = FALSE)
gleason.altered <- function(input){
  last.col <- ncol(input)+1
# For each gene in the summary table:
  for (k in 1:nrow(input)){
    gene <- input$`Hugo Symbol`[k]
    # Obtain sample names for primary subtype
    sample.names <- colnames(TCGA.coloss)[2:ncol(TCGA.coloss)]
    # Create data frame of rows = samples | columns = gleason score, risk group, Gene A expression
    data <- data.frame()
    for (i in 1:length(sample.names)){
      index <- which(sample.names[i] == clinical.sample$SAMPLE_ID)
      data[i,1] <- clinical.sample$REVIEWED_GLEASON_SUM[index]
      if (data[i,1] > 7) {data[i,2] <- "High"}
      if (data[i,1] == 7) {data[i,2] <- "Intermediate"}
      if (data[i,1] < 7) {data[i,2] <- "Low"}
      }
    for (i in 1:length(sample.names)){
      index <- which(sample.names[i] == colnames(TCGA.mRNA))
      index2 <- which(gene == TCGA.mRNA$Hugo_Symbol)[1]
      data[i,3] <- TCGA.mRNA[index2,index]
    }
    colnames(data) <- c("Score", "Risk Group", "Expression")
# Order dataset by gene expression. Subset into high and low expressors.
    data <- data[order(data[,2], decreasing = TRUE),]
    # Stratify patients into n=3 quantiles.
    quantiles <- quantile(data$Expression, c(0.5))
    high.expressing <- data[data$Expression >= quantiles[1],]
    low.expressing <- data[data$Expression <= quantiles[1],]
# Compute high-risk enrichment
    # Construct contingency table
    high.expressing.high.risk <- nrow(high.expressing[high.expressing$`Risk Group` == "High",])
    high.expressing.not.high.risk <- nrow(high.expressing) - high.expressing.high.risk
    low.expressing.high.risk <- nrow(low.expressing[low.expressing$`Risk Group` == "High",])
    low.expressing.not.high.risk <- nrow(low.expressing) - low.expressing.high.risk
    contignecy.table <- matrix(data = c(high.expressing.high.risk, high.expressing.not.high.risk,
                                        low.expressing.high.risk, low.expressing.not.high.risk), 
                              ncol = 2, nrow = 2, byrow = TRUE)
    # Compute high-risk gleason enrichment between high and low expressers
    fishers.test <- fisher.test(contignecy.table)
    pval <- fishers.test$p.value
# Add high risk pval to summary table
    input[k,last.col] <- pval
  }
  colnames(input)[last.col] <- "Subtype High-Risk Enrichment p-val"
  return(input)
}
gleason.unaltered <- function(input){
  last.col <- ncol(input)+1
  # For each gene in the summary table:
  for (k in 1:nrow(input)){
    gene <- input$`Hugo Symbol`[k]
    # Obtain sample names for primary subtype
    sample.names <- colnames(TCGA.WT)[2:ncol(TCGA.WT)]
    # Create data frame of rows = samples | columns = gleason score, risk group, Gene A expression
    data <- data.frame()
    for (i in 1:length(sample.names)){
      index <- which(sample.names[i] == clinical.sample$SAMPLE_ID)
      data[i,1] <- clinical.sample$REVIEWED_GLEASON_SUM[index]
      if (data[i,1] > 7) {data[i,2] <- "High"}
      if (data[i,1] == 7) {data[i,2] <- "Intermediate"}
      if (data[i,1] < 7) {data[i,2] <- "Low"}
    }
    for (i in 1:length(sample.names)){
      index <- which(sample.names[i] == colnames(TCGA.mRNA))
      index2 <- which(gene == TCGA.mRNA$Hugo_Symbol)[1]
      data[i,3] <- TCGA.mRNA[index2,index]
    }
    colnames(data) <- c("Score", "Risk Group", "Expression")
    # Order dataset by gene expression. Subset into high and low expressors.
    data <- data[order(data[,2], decreasing = TRUE),]
    # Stratify patients into n=3 quantiles.
    quantiles <- quantile(data$Expression, c(0.5))
    high.expressing <- data[data$Expression >= quantiles[1],]
    low.expressing <- data[data$Expression <= quantiles[1],]
    # Compute high-risk enrichment
    # Construct contingency table
    high.expressing.high.risk <- nrow(high.expressing[high.expressing$`Risk Group` == "High",])
    high.expressing.not.high.risk <- nrow(high.expressing) - high.expressing.high.risk
    low.expressing.high.risk <- nrow(low.expressing[low.expressing$`Risk Group` == "High",])
    low.expressing.not.high.risk <- nrow(low.expressing) - low.expressing.high.risk
    contignecy.table <- matrix(data = c(high.expressing.high.risk, high.expressing.not.high.risk,
                                        low.expressing.high.risk, low.expressing.not.high.risk), 
                               ncol = 2, nrow = 2, byrow = TRUE)
    # Compute high-risk gleason enrichment between high and low expressers
    fishers.test <- fisher.test(contignecy.table)
    pval <- fishers.test$p.value
    # Add high risk pval to summary table
    input[k,last.col] <- pval
  }
  colnames(input)[last.col] <- "WT High-Risk Enrichment p-val"
  return(input)
}