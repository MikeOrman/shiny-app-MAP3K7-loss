# dataset = primary.alteration or met.alteration data frame
# mRNA = formatted mRNA data
#--------------------------------------COMPUTE GENOME-WIDE DGE FDR--------------
genome_wide_DGE <- function(dataset, mRNA){
  DGE <- data.frame()
  for (i in 1:nrow(mRNA)){
    gene <- mRNA$Hugo_Symbol[i]
    DGE[i,1] <- gene
    # Subset patients in dataset by alteration of gene mRNA data frame
    altered <- alteration_subset(dataset, c(gene, "altered"))
    unaltered <- alteration_subset(dataset, c(gene, "unaltered"))
    #---------------------------COMPUTE DGE BY STUDENT'S T-TEST-----------------
    if (class(altered) == "data.frame" & class(unaltered) == "data.frame"){
      altered.patients <- colnames(altered)[2:ncol(altered)]
      unaltered.patients <- colnames(unaltered)[2:ncol(unaltered)]
      # Subset mRNA data
      altered.GEX.patients <- intersect(altered.patients, colnames(mRNA)[2:ncol(mRNA)])
      unaltered.GEX.patients <- intersect(unaltered.patients, colnames(mRNA)[2:ncol(mRNA)])
      index <- c(TRUE)
      for (k in 2:ncol(mRNA)){
        if (sum(colnames(mRNA)[k] == altered.GEX.patients) == 1){
          index[k] <- TRUE} else {index[k] <- FALSE}
      }
      altered.GEX <- mRNA[, index]
      index <- c(TRUE)
      for (k in 2:ncol(mRNA)){
        if (sum(colnames(mRNA)[k] == unaltered.GEX.patients) == 1){
          index[k] <- TRUE} else {index[k] <- FALSE}
      }
      unaltered.GEX <- mRNA[, index]
      # Ensure at least 3 samples are present for DGE calculation
      if (class(altered.GEX) == "data.frame" & class(unaltered.GEX) == "data.frame"){
        if (ncol(altered.GEX) >= 4 & ncol(unaltered.GEX) >= 4) {
          # Correct row names
          rownames(altered.GEX) <- 1:nrow(altered.GEX)
          rownames(unaltered.GEX) <- 1:nrow(unaltered.GEX)
          # Find index for concordant gene being tested
          index <- which(gene == altered.GEX$Hugo_Symbol)[1]
          # Compute fold change
          DGE[i,2] <- (mean(as.numeric(altered.GEX[index,2:ncol(altered.GEX)])) / 
                                    mean(as.numeric(unaltered.GEX[index,2:ncol(unaltered.GEX)])))
          # Compute pval
          DGE[i,3] <- t.test(as.numeric(altered.GEX[index,2:ncol(altered.GEX)]), 
                                        (as.numeric(unaltered.GEX[index,2:ncol(unaltered.GEX)])))$p.value
        }
      }
    }
  }
  DGE[,4] <- p.adjust(DGE[,3], method = "fdr")
  colnames(DGE) <- c("Hugo Symbol", "FC (altered:unaltered)", "p-val", "FDR")
  return(DGE)
}
# summary = summary table from genomic analysis
# input = output from genome_wide_DGE function
#----------ADD DGE STATS TO SUMMARY TABLE | FILTER FOR CONCORDANCE--------------
# Input = output from the above function
concordant_DGE <- function(summary, input){
# Initialize concordant DGE data frame
  concordant.DGE <- data.frame()
  for (i in 1:nrow(summary)){
    concordant.DGE[i,1] <- summary$`Hugo Symbol`[i]
    concordant.DGE[i,2] <- summary$Alteration[i]
    # Find index in genome-wide table
    index <- which(concordant.DGE[i,1] == input[,1])[1]
    if (concordant.DGE[i,2] == "LOF" & input[index,2] < 1 & is.na(input[index,2]) == FALSE){
      concordant.DGE[i,3] <- input[index,2]
      concordant.DGE[i,4] <- input[index,4]}
    if (concordant.DGE[i,2] == "GOF" & input[index,2] > 1 & is.na(input[index,2]) == FALSE){
      concordant.DGE[i,3] <- input[index,2]
      concordant.DGE[i,4] <- input[index,4]}
    }
  colnames(concordant.DGE) <- c("Hugo Symbol", "Alteration", "Concordant FC (altered:unaltered)", "FDR")
  concordant.DGE <- na.omit(concordant.DGE)
  return(concordant.DGE)
  }