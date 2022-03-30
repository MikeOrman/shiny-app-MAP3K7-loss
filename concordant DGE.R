# Computes concordant gene expression given the summary table from the genomic analysis output
# summary = formatted summary table from Expression Analysis.Rmd
# background = CNA profile of background context for testing concordant expression
# mRNA = formatted mRNA data from Expression Analysis.Rmd
# A value of "NA" in the output means there were not enough samples to calculate concordant DGE on the desired background
library(ggplot2)
library(ggpubr)
concordant_DGE <- function(summary, background, mRNA){
# INITIALIZE DATA FRAMES
  concordant.DGE <- data.frame()
# FOR EACH ROW OF THE SUMMARY TABLE:
  for (i in 1:nrow(summary)){
    concordant.DGE[i,1] <- summary$`Hugo Symbol`[i]
#---------------------------SUBSET PATIENTS-------------------------------------
  if (summary$Alteration[i] == "LOF"){
    concordant.DGE[i,2] <- "LOF"
    altered <- subtype_subset(background, c(summary$`Hugo Symbol`[i], "loss"))
    unaltered <- subtype_subset(background, c(summary$`Hugo Symbol`[i], "WT"))
  }
  if (summary$Alteration[i] == "GOF"){
    concordant.DGE[i,2] <- "GOF"
    altered <- subtype_subset(background, c(summary$`Hugo Symbol`[i], "gain"))
    unaltered <- subtype_subset(background, c(summary$`Hugo Symbol`[i], "WT"))
  }
#---------------------------COMPUTE DGE BY STUDENT'S T-TEST---------------------
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
          # Find index for concordant gene being tested
          index <- which(summary$`Hugo Symbol`[i] == altered.GEX$Hugo_Symbol)[1]
          # Compute concordant fold change
          concordant.DGE[i,3] <- (mean(as.numeric(altered.GEX[index,2:ncol(altered.GEX)])) / 
                             mean(as.numeric(unaltered.GEX[index,2:ncol(unaltered.GEX)])))
          # Compute concordant pval
          concordant.DGE[i,4] <- t.test(as.numeric(altered.GEX[index,2:ncol(altered.GEX)]), 
                                 (as.numeric(unaltered.GEX[index,2:ncol(unaltered.GEX)])))$p.value
        }
      }
    }
  }
  colnames(concordant.DGE) <- c("Hugo Symbol", "Alteration", "Concordant FC (altered:unaltered)", "p-val")
  return(concordant.DGE)
}