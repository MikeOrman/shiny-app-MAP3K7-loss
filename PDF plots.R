# Output plots to PDF. Can be run independently of the "Shiny app V3.Rmd"
source("plotting functions.R")
# Rank relative to other genes
gene.rank.dotplot(summary.full)
# Cytogenetic location relative to other genes
karyoplots(summary.full)
# Genomic
genomic.plots(summary.full)
# Concordant DGE
concordant.DGE.boxplot(summary.full)
# Tumor Grade
gleason.barplot(summary.full)
#---------------------------------------SURVIVAL BENEFIT------------------------
# Unfortunately, the survfit function does not work when called inside of a custom function. As a work-around, 
# I have placed the full function here
input <- summary.full
  clinical.sample <- read.table("TCGA clinical sample.txt", header = TRUE, sep = "\t", check.names = FALSE)
  clinical.patient <- read.table("TCGA clinical patient.txt", header = TRUE, sep = "\t", check.names = FALSE)
  pdf("Survival plots.pdf")
  # For each gene in the summary table:
  for (k in 1:nrow(input)){
    gene <- input$`Hugo Symbol`[k]
    # Obtain sample names for primary subtype
    sample.names <- colnames(TCGA.coloss)[2:ncol(TCGA.coloss)]
    sample.names.WT <- colnames(TCGA.WT)[2:ncol(TCGA.WT)]
    # Create  coloss data frame of rows = samples | columns = Patient ID, Group, Time, Status
    data <- data.frame()
    for (i in 1:length(sample.names)){
      index <- which(sample.names[i] == clinical.sample$SAMPLE_ID)
      sample <- clinical.sample$PATIENT_ID[index]
      index2 <- which(sample == clinical.patient$PATIENT_ID)
      data[i,1] <- clinical.patient$PATIENT_ID[index2]
      if (clinical.patient$DFS_MONTHS[index2] != "[Not Available]") {data[i,2] <- clinical.patient$DFS_MONTHS[index2]}
      if (clinical.patient$DFS_STATUS[index2] != "[Not Available]") {data[i,3] <- clinical.patient$DFS_STATUS[index2]}
    }
    for (i in 1:length(sample.names)){
      index <- which(sample.names[i] == colnames(TCGA.mRNA))
      index2 <- which(gene == TCGA.mRNA$Hugo_Symbol)
      data[i,4] <- TCGA.mRNA[index2,index]
    }
    data <- na.omit(data)
    colnames(data) <- c("Patient ID", "DFS Months", "DFS Status", "Expression")
    # Rename status column
    for (i in 1:nrow(data)){
      if (data[i,3] == "0:DiseaseFree") {data[i,3] <- 0}
      if (data[i,3] == "1:Recurred/Progressed") {data[i,3] <- 1}
    }
    # Order dataset by gene expression. Subset into high and low expressors
    data <- data[order(data$Expression, decreasing = TRUE),]
    # Stratify patients into n=2 quantiles.
    quantiles <- quantile(data$Expression, 0.50)
    high.expressing <- data[data$Expression >= quantiles[1],]
    high.expressing[,5] <- paste("High", paste(gene, "mRNA Expression", sep = " "), sep = " ")
    colnames(high.expressing)[5] <- "Group"
    low.expressing <- data[data$Expression <= quantiles[1],]
    low.expressing[,5] <- paste("Low", paste(gene, "mRNA Expression", sep = " "), sep = " ")
    colnames(low.expressing)[5] <- "Group"
    # Prepare dataframe for survival analysis
    high.expressing <- high.expressing[,c(1, 5, 2, 3)]
    low.expressing <- low.expressing[,c(1, 5, 2, 3)]
    surv_data <- rbind(high.expressing, low.expressing)
    surv_data$`DFS Months` <- as.numeric(surv_data$`DFS Months`)
    surv_data$`DFS Status` <- as.numeric(surv_data$`DFS Status`)
    # Make coloss surv object
    surv_object <- Surv(time = surv_data$`DFS Months`, event = surv_data$`DFS Status`)
    # Create  WT data frame of rows = samples | columns = Patient ID, Group, Time, Status
    data <- data.frame()
    for (i in 1:length(sample.names.WT)){
      index <- which(sample.names.WT[i] == clinical.sample$SAMPLE_ID)
      sample <- clinical.sample$PATIENT_ID[index]
      index2 <- which(sample == clinical.patient$PATIENT_ID)
      data[i,1] <- clinical.patient$PATIENT_ID[index2]
      if (clinical.patient$DFS_MONTHS[index2] != "[Not Available]") {data[i,2] <- clinical.patient$DFS_MONTHS[index2]}
      if (clinical.patient$DFS_STATUS[index2] != "[Not Available]") {data[i,3] <- clinical.patient$DFS_STATUS[index2]}
    }
    for (i in 1:length(sample.names.WT)){
      index <- which(sample.names.WT[i] == colnames(TCGA.mRNA))
      index2 <- which(gene == TCGA.mRNA$Hugo_Symbol)
      data[i,4] <- TCGA.mRNA[index2,index]
    }
    data <- na.omit(data)
    colnames(data) <- c("Patient ID", "DFS Months", "DFS Status", "Expression")
    # Rename status column
    for (i in 1:nrow(data)){
      if (data[i,3] == "0:DiseaseFree") {data[i,3] <- 0}
      if (data[i,3] == "1:Recurred/Progressed") {data[i,3] <- 1}
    }
    # Order dataset by gene expression. Subset into high and low expressors
    data <- data[order(data$Expression, decreasing = TRUE),]
    # Stratify patients into n=2 quantiles.
    quantiles <- quantile(data$Expression, 0.50)
    high.expressing <- data[data$Expression >= quantiles[1],]
    high.expressing[,5] <- paste("High", paste(gene, "mRNA Expression", sep = " "), sep = " ")
    colnames(high.expressing)[5] <- "Group"
    low.expressing <- data[data$Expression <= quantiles[1],]
    low.expressing[,5] <- paste("Low", paste(gene, "mRNA Expression", sep = " "), sep = " ")
    colnames(low.expressing)[5] <- "Group"
    # Prepare dataframe for survival analysis
    high.expressing <- high.expressing[,c(1, 5, 2, 3)]
    low.expressing <- low.expressing[,c(1, 5, 2, 3)]
    surv_data.WT <- rbind(high.expressing, low.expressing)
    surv_data.WT$`DFS Months` <- as.numeric(surv_data.WT$`DFS Months`)
    surv_data.WT$`DFS Status` <- as.numeric(surv_data.WT$`DFS Status`)
    # Make WT surv object
    surv_object.WT <- Surv(time = surv_data.WT$`DFS Months`, event = surv_data.WT$`DFS Status`)
    # Make plots
    fit <- survival::survfit(surv_object ~ Group, data = surv_data)
    fit.WT <- survival::survfit(surv_object.WT ~ Group, data = surv_data.WT)
    print((ggsurvplot(fit, pval = TRUE, pval.method = TRUE,
                      risk.table = TRUE,
                      palette = c("#E7B800", "#2E9FDF"),
                      risk.table.y.text = FALSE,
                      ylab = "DFS Ratio",
                      xlab = "Months",
                      title = "Primary Subtype Patients")))
    print((ggsurvplot(fit.WT, pval = TRUE, pval.method = TRUE,
                      risk.table = TRUE,
                      palette = c("#E7B800", "#2E9FDF"),
                      risk.table.y.text = FALSE,
                      ylab = "DFS Ratio",
                      xlab = "Months",
                      title = "Primary WT Patients")))
  }
  dev.off()