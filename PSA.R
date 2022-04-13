# Computes enrichment of high-risk gleason scores in primary patients stratified by gene expression
# Genes in binary matrix and genes in mRNA data must be harmonized
# Input = binary alteration matrix
library(rstatix)
clinical.sample <- read.table("TCGA clinical sample.txt", header = TRUE, sep = "\t", check.names = FALSE)
clinical.patient <- read.table("TCGA clinical patient.txt", header = TRUE, sep = "\t", check.names = FALSE)
PSA <- function(input){
  output <- data.frame()
  # For each gene in the summary table:
  for (k in 1:nrow(input)){
    gene <- input$`Hugo_Symbol`[k]
    # Plot subtype high-risk enrichment
    sample.names <- colnames(input)[2:ncol(input)]
    # Create data frame of rows = samples | columns = gleason score, risk group, Gene A expression
    data <- data.frame()
    for (i in 1:length(sample.names)){
      patient.ID.index <- which(sample.names[i] == clinical.sample$SAMPLE_ID)
      patient.ID <- clinical.sample$PATIENT_ID[patient.ID.index]
      patient.index <- which(patient.ID == clinical.patient$PATIENT_ID)
      patient.PSA <- (clinical.patient$PSA_MOST_RECENT_RESULTS[patient.index])
      if (patient.PSA == "[Not Available]") {patient.PSA <- NA}
      data[i,1] <- patient.ID
      data[i,2] <- as.numeric(patient.PSA)
      mRNA.sample.index <- which(sample.names[i] == colnames(TCGA.mRNA))
      mRNA.gene.index <- which(gene == TCGA.mRNA$Hugo_Symbol)
      data[i,3] <- TCGA.mRNA[mRNA.gene.index,mRNA.sample.index]
    }
    data <- na.omit(data)
    colnames(data) <- c("Patient ID", "PSA", paste(gene, "mRNA Expression", sep = " "))
    # Order dataset by gene expression. Subset into high and low expressers.
    data <- data[order(data[,3], decreasing = TRUE),]
    # Stratify patients into n=2 quantiles.
    quantiles <- quantile(data[,3], c(0.50))
    high.expressing <- data[data[,3] >= quantiles[1],]
    high.expressing[,4] <- paste("High", gene, "Expressers", sep = " ")
    colnames(high.expressing)[4] <- "Group"
    low.expressing <- data[data[,3] <= quantiles[1],]
    low.expressing[,4] <- paste("Low", gene, "Expressers", sep = " ")
    colnames(low.expressing)[4] <- "Group"
    data <- rbind(low.expressing, high.expressing)
    # Calculate t-test
    t.test <- data %>% t_test(PSA ~ Group, paired = FALSE)
    # Add result to output table
    output[k,1] <- gene
    output[k,2] <- t.test$statistic
    output[k,3] <- t.test$p
  }
  output[,4] <- p.adjust(output[,3], method = "fdr")
  colnames(output) <- c("Hugo Symbol", "PSA t-test statistic (High Expressers - Low Expressers", "p", "FDR")
  return(output)
}