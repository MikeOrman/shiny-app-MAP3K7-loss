library(ggplot2)
library(ggrepel)
source("manhattan.R")
library(ggpubr)
library(gridExtra)
library(survival)
library(survminer)
# Note: Run this after running all on "shiny app V3 - MAP3K7 loss.Rmd"
#---------------------------------------LOAD WORKSPACE--------------------------
load("~/Github/Projects/shiny app MAP3K7 loss/shiny app V3 - MAP3K7 loss.RData")
#---------------------------------------RANKED LIST-----------------------------
gene.rank.dotplot <- function(input){
  pdf("Gene Rank.pdf")
  h = -log10(0.05)
  for (i in 1:nrow(input)){
    gene <- input$`Hugo Symbol`[i]
    gene.layer <- subset(input,  `Hugo Symbol` == gene)
    print((ggplot(input, aes(`Gene Rank`, `-log10(Molecular FDR)`, label = `Hugo Symbol`))
      + geom_point() 
      + geom_point(data = gene.layer, colour="red")
      + geom_text_repel(data = subset(input, `Hugo Symbol` == gene), color = "red", nudge_x = 0, nudge_y = .3)
      + geom_hline(yintercept = h , linetype = "solid", color = "red", size = 0.8)
      + geom_text(aes(0,h,label = "FDR = 0.05", vjust = -1, hjust = -6))))
    }
  dev.off()
}
#---------------------------------------KARYOPLOT-------------------------------
karyoplots <- function(input){
  # Create manhattan object
  mp <- manhattan.object(input[,c(1,4)])
  pdf("Karyoplots.pdf")
  for (i in 1:nrow(input)){
    if (is.na(input$`Cytogenetic Band`[i]) == FALSE){
      gene <- input$`Hugo Symbol`[i]
      if (sum(mp@ranges@NAMES == gene) > 0){
        # Create manhattan object for labelling gene
        mp.label <- mp[gene]
        # Plot manhattan object
        kp <- kpPlotManhattan(karyoplot = plotKaryotype(plot.type = 2), 
                              data = mp,
                              points.col = "Navy",
                              highlight = gene,
                              highlight.col = "red",
                              genomewide.col = "white",
                              suggestive.col = "white")
        kpText(kp, mp.label, labels = names(mp.label), y = 1)
        print(kp)
      }
    }
  }
  dev.off()
  }
#---------------------------------------GENETIC ALTERATIONS---------------------
genomic.plots <- function(summary.full){
  pdf(file = "Alteration bar plots.pdf")
  for (i in 1:nrow(summary.full)){
    gene <- summary.full$`Hugo Symbol`[i]
    if (summary.full$Alteration[i] == "LOF"){
      primary.WT <- primary.WT.LOF.frequency$`LOF alteration rate: Primary WT patients`[which(primary.WT.LOF.frequency$Hugo_Symbol == gene)]
      primary.coloss <- primary.coloss.LOF.frequency$`LOF alteration rate: Primary coloss patients`[which(primary.coloss.LOF.frequency$Hugo_Symbol == gene)]
      met.WT <- met.WT.LOF.frequency$`LOF alteration rate: Metastatic WT patients`[which(met.coloss.LOF.frequency$Hugo_Symbol == gene)]
      met.coloss <- met.coloss.LOF.frequency$`LOF alteration rate: Metastatic coloss patients`[which(met.coloss.LOF.frequency$Hugo_Symbol == gene)]
      x <- c("Primary WT", "Primary Subtype", "Metastatic WT", "Metastatic Subtype")
      y <- c(primary.WT, primary.coloss, met.WT, met.coloss)
      genomic.plot.data <- data.frame(x, y)
      genomic.plot.data$x <- factor(genomic.plot.data$x, levels = genomic.plot.data$x)
      stats <- t(c("Freq", "Primary Subtype", "Metastatic Subtype",
                   paste("p=", formatC(summary.full$`Alteration Frequency p-val`[i], format = "e", digits = 2), sep = "")))
      stats <- as.data.frame(stats)
      colnames(stats) <- c(".y.", "group1", "group2", "p")
      print((ggplot(genomic.plot.data, aes(x, y)) + geom_bar(stat='identity') + ggtitle(paste(gene, "LOF Alterations", sep = " ")) 
        + xlab("Genetic Background") + ylab("Alteration Frequency") + ylim(0,1)
        + stat_pvalue_manual(stats, label = "p", tip.length = 0.01, 
                             y.position = max(y) + .0625)))
      }
    if (summary.full$Alteration[i] == "GOF"){
      primary.WT <- primary.WT.GOF.frequency$`GOF alteration rate: Primary WT patients`[which(primary.WT.GOF.frequency$Hugo_Symbol == gene)]
      primary.coloss <- primary.coloss.GOF.frequency$`GOF alteration rate: Primary coloss patients`[which(primary.coloss.GOF.frequency$Hugo_Symbol == gene)]
      met.WT <- met.WT.GOF.frequency$`GOF alteration rate: Metastatic WT patients`[which(met.coloss.GOF.frequency$Hugo_Symbol == gene)]
      met.coloss <- met.coloss.GOF.frequency$`GOF alteration rate: Metastatic coloss patients`[which(met.coloss.GOF.frequency$Hugo_Symbol == gene)]
      x <- c("Primary WT", "Primary Subtype", "Metastatic WT", "Metastatic Subtype")
      y <- c(primary.WT, primary.coloss, met.WT, met.coloss)
      genomic.plot.data <- data.frame(x, y)
      genomic.plot.data$x <- factor(genomic.plot.data$x, levels = genomic.plot.data$x)
      stats <- t(c("Freq", "Primary Subtype", "Metastatic Subtype",
                   paste("p=", formatC(summary.full$`Alteration Frequency p-val`[i], format = "e", digits = 2), sep = "")))
      stats <- as.data.frame(stats)
      colnames(stats) <- c(".y.", "group1", "group2", "p")
      print((ggplot(genomic.plot.data, aes(x, y)) + geom_bar(stat='identity') + ggtitle(paste(gene, "GOF Alterations", sep = " ")) 
             + xlab("Genetic Background") + ylab("Alteration Frequency") + ylim(0,1)
             + stat_pvalue_manual(stats, label = "p", tip.length = 0.01, 
                                  y.position = max(y) + .0625)))
    }
  }
  dev.off()
}
#---------------------------------------CONCORDANT DGE--------------------------
concordant.DGE.boxplot <- function(summary.full){
  DGE.boxplot <- data.frame()
  pdf(file = "Concordant gene expression plots.PDF")
  for (i in 1:nrow(summary.full)){
    gene <- summary.full$`Hugo Symbol`[i]
    DGE.boxplot[i,1] <- summary.full$`Hugo Symbol`[i]
    #---------------------------SUBSET PATIENTS-------------------------------------
    if (summary.full$Alteration[i] == "LOF"){
      DGE.boxplot[i,2] <- "LOF"
      primary.altered <- subtype_subset(TCGA.coloss, c(summary.full$`Hugo Symbol`[i], "loss"))
      primary.unaltered <- subtype_subset(TCGA.coloss, c(summary.full$`Hugo Symbol`[i], "WT"))
      primary.altered.name <- paste(summary.full$`Hugo Symbol`[i], "Altered", sep = " ")
      primary.unaltered.name <- paste(summary.full$`Hugo Symbol`[i], "WT", sep = " ")
      met.altered <- subtype_subset(met.coloss, c(summary.full$`Hugo Symbol`[i], "loss"))
      met.unaltered <- subtype_subset(met.coloss, c(summary.full$`Hugo Symbol`[i], "WT"))
      met.altered.name <- paste(summary.full$`Hugo Symbol`[i], "Altered", sep = " ")
      met.unaltered.name <- paste(summary.full$`Hugo Symbol`[i], "WT", sep = " ")
    }
    if (summary.full$Alteration[i] == "GOF"){
      DGE.boxplot[i,2] <- "GOF"
      primary.altered <- subtype_subset(TCGA.coloss, c(summary.full$`Hugo Symbol`[i], "gain"))
      primary.unaltered <- subtype_subset(TCGA.coloss, c(summary.full$`Hugo Symbol`[i], "WT"))
      primary.altered.name <- paste(summary.full$`Hugo Symbol`[i], "Altered", sep = " ")
      primary.unaltered.name <- paste(summary.full$`Hugo Symbol`[i], "WT", sep = " ")
      met.altered <- subtype_subset(met.coloss, c(summary.full$`Hugo Symbol`[i], "gain"))
      met.unaltered <- subtype_subset(met.coloss, c(summary.full$`Hugo Symbol`[i], "WT"))
      met.altered.name <- paste(summary.full$`Hugo Symbol`[i], "Altered", sep = " ")
      met.unaltered.name <- paste(summary.full$`Hugo Symbol`[i], "WT", sep = " ")
    }
    #---------------------------BOXPLOTS OF CONCORDANT DGE--------------------------
    primary.unaltered.patients <- colnames(primary.unaltered)[2:ncol(primary.unaltered)]
    primary.altered.patients <- colnames(primary.altered)[2:ncol(primary.altered)]
    met.unaltered.patients <- colnames(met.unaltered)[2:ncol(met.unaltered)]
    met.altered.patients <- colnames(met.altered)[2:ncol(met.altered)]
    # Subset TCGA mRNA data  
      primary.altered.GEX.patients <- intersect(primary.altered.patients, colnames(TCGA.mRNA)[2:ncol(TCGA.mRNA)])
      primary.unaltered.GEX.patients <- intersect(primary.unaltered.patients, colnames(TCGA.mRNA)[2:ncol(TCGA.mRNA)])
      index <- c(TRUE)
      for (k in 2:ncol(TCGA.mRNA)){
        if (sum(colnames(TCGA.mRNA)[k] == primary.altered.GEX.patients) == 1){
          index[k] <- TRUE} else {index[k] <- FALSE}
      }
      primary.altered.GEX <- TCGA.mRNA[, index]
      index <- c(TRUE)
      for (k in 2:ncol(TCGA.mRNA)){
        if (sum(colnames(TCGA.mRNA)[k] == primary.unaltered.GEX.patients) == 1){
          index[k] <- TRUE} else {index[k] <- FALSE}
      }
      primary.unaltered.GEX <- TCGA.mRNA[, index]
      # Subset Abida mRNA data
      met.altered.GEX.patients <- intersect(met.altered.patients, colnames(Abida.mRNA)[2:ncol(Abida.mRNA)])
      met.unaltered.GEX.patients <- intersect(met.unaltered.patients, colnames(Abida.mRNA)[2:ncol(Abida.mRNA)])
      index <- c(TRUE)
      for (k in 2:ncol(Abida.mRNA)){
        if (sum(colnames(Abida.mRNA)[k] == met.altered.GEX.patients) == 1){
          index[k] <- TRUE} else {index[k] <- FALSE}
      }
      met.altered.GEX <- Abida.mRNA[, index]
      index <- c(TRUE)
      for (k in 2:ncol(Abida.mRNA)){
        if (sum(colnames(Abida.mRNA)[k] == met.unaltered.GEX.patients) == 1){
          index[k] <- TRUE} else {index[k] <- FALSE}
      }
      met.unaltered.GEX <- Abida.mRNA[, index]
      # Combine primary data into single data frame
      index <- which(summary.full$`Hugo Symbol`[i] == primary.altered.GEX$Hugo_Symbol)[1]
      primary.altered.vec <- as.numeric(primary.altered.GEX[index,2:ncol(primary.altered.GEX)])
      primary.unaltered.vec <- as.numeric(primary.unaltered.GEX[index,2:ncol(primary.unaltered.GEX)])
      data.primary <- data.frame()
      data.primary[1:length(primary.unaltered.vec),1] <- primary.unaltered.name
      data.primary[1:length(primary.unaltered.vec),2] <- primary.unaltered.vec
      data.primary[((length(primary.unaltered.vec)+1):(length(primary.unaltered.vec)+length(primary.altered.vec))),1] <- primary.altered.name
      data.primary[((length(primary.unaltered.vec)+1):(length(primary.unaltered.vec)+length(primary.altered.vec))),2] <- primary.altered.vec
      data.primary[,3] <- "Primary Subtype"
      colnames(data.primary) <- c("Alteration Status",
                                  paste(gene, "Log2 Normalized Expression", sep = " "),
                                  "Background")
      data.primary$`Alteration Status` <- factor(data.primary$`Alteration Status`, 
                                                 levels = c((paste(gene, "WT", sep = " ")), (paste(gene, "Altered", sep = " "))))
      # Combine metastatic data into single data frame
      index <- which(summary.full$`Hugo Symbol`[i] == met.altered.GEX$Hugo_Symbol)[1]
      met.altered.vec <- as.numeric(met.altered.GEX[index,2:ncol(met.altered.GEX)])
      met.unaltered.vec <- as.numeric(met.unaltered.GEX[index,2:ncol(met.unaltered.GEX)])
      data.met <- data.frame()
      data.met[1:length(met.unaltered.vec),1] <- met.unaltered.name
      data.met[1:length(met.unaltered.vec),2] <- met.unaltered.vec
      data.met[((length(met.unaltered.vec)+1):(length(met.unaltered.vec)+length(met.altered.vec))),1] <- met.altered.name
      data.met[((length(met.unaltered.vec)+1):(length(met.unaltered.vec)+length(met.altered.vec))),2] <- met.altered.vec
      data.met[,3] <- "Metastatic Subtype"
      colnames(data.met) <- c("Alteration Status",
                              paste(gene, "Log2 Normalized Expression", sep = " "),
                              "Background")
      data.met$`Alteration Status` <- factor(data.met$`Alteration Status`, 
                                                 levels = c((paste(gene, "WT", sep = " ")), (paste(gene, "Altered", sep = " "))))
      # Bind primary and metastatic data
      data <- rbind(data.primary, data.met)
      data$Background <- factor(data$Background,
                                levels = c("Primary Subtype", "Metastatic Subtype"))
      # Make boxplot
      print((ggplot(data, aes(`Alteration Status`, data[,2], fill = `Alteration Status`))
             + geom_boxplot() + geom_jitter() 
             + ylab(paste(gene, "log2(Normalized mRNA Expression)"))
             + facet_wrap(~Background)) + stat_compare_means(method = "t.test", vjust = -1, hjust = 0))
      }
  dev.off()
}
#---------------------------------------TUMOR GRADE-----------------------------
gleason.barplot <- function(input){
  pdf("Gleason scores.pdf")
  clinical.sample <- read.table("TCGA clinical sample.txt", header = TRUE, sep = "\t", check.names = FALSE)
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
      index2 <- which(gene == TCGA.mRNA$Hugo_Symbol)
      data[i,3] <- TCGA.mRNA[index2,index]
      }
    colnames(data) <- c("Score", "Risk Group", "Expression")
    # Order dataset by gene expression. Subset into high and low expressors.
    data <- data[order(data[,2], decreasing = TRUE),]
    # Stratify patients into n=2 quantiles.
    quantiles <- quantile(data$Expression, c(0.50))
    high.expressing <- data[data$Expression > quantiles[1],]
    high.expressing[,4] <- "High Expressers"
    colnames(high.expressing)[4] <- "mRNA Expression"
    low.expressing <- data[data$Expression <= quantiles[1],]
    low.expressing[,4] <- "Low Expressers"
    colnames(low.expressing)[4] <- "mRNA Expression"
    data <- rbind(low.expressing, high.expressing)
    data$`mRNA Expression` <- factor(data$`mRNA Expression`, levels = c("Low Expressers", "High Expressers"))
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
    stats <- t(c("Count", "Low Expressers", "High Expressers",
                 paste("p=", formatC(pval, format = "e", digits = 2), sep = "")))
    stats <- as.data.frame(stats)
    colnames(stats) <- c(".y.", "group1", "group2", "p")
    print((ggplot(data, aes(`mRNA Expression`)) + geom_bar(aes(fill=`Risk Group`)) 
           + labs(title = "Primary Subtype Patients", x = paste(gene, "mRNA Expression"), y = "Count")
           + ylim(0,(nrow(data)*.60)) + stat_pvalue_manual(stats, label = "p", y.position = nrow(data)*.55)))
  }
  dev.off()
}
#---------------------------------------ALL PLOTS-------------------------------
all.plots <- function(summary.full){
  pdf("All plots.pdf")
  par(mfcol = c(2,2))
  clinical.sample <- read.table("TCGA clinical sample.txt", header = TRUE, sep = "\t", check.names = FALSE)
  for (i in 1:nrow(summary.full)){
    print(i)
    gene <- summary.full$`Hugo Symbol`[i]
    #---------------------------------------GENETIC ALTERATIONS---------------------
    if (summary.full$Alteration[i] == "LOF"){
      primary.WT <- primary.WT.LOF.frequency$`LOF alteration rate: Primary WT patients`[which(primary.WT.LOF.frequency$Hugo_Symbol == gene)]
      primary.coloss <- primary.coloss.LOF.frequency$`LOF alteration rate: Primary coloss patients`[which(primary.coloss.LOF.frequency$Hugo_Symbol == gene)]
      met.WT <- met.WT.LOF.frequency$`LOF alteration rate: Metastatic WT patients`[which(met.coloss.LOF.frequency$Hugo_Symbol == gene)]
      met.coloss <- met.coloss.LOF.frequency$`LOF alteration rate: Metastatic coloss patients`[which(met.coloss.LOF.frequency$Hugo_Symbol == gene)]
      x <- c("Primary WT", "Primary Subtype", "Metastatic WT", "Metastatic Subtype")
      y <- c(primary.WT, primary.coloss, met.WT, met.coloss)
      genomic.plot.data <- data.frame(x, y)
      genomic.plot.data$x <- factor(genomic.plot.data$x, levels = genomic.plot.data$x)
      stats <- t(c("Freq", "Primary Subtype", "Metastatic Subtype",
                   paste("p=", formatC(summary.full$`Alteration Frequency p-val`[i], format = "e", digits = 2), sep = "")))
      stats <- as.data.frame(stats)
      colnames(stats) <- c(".y.", "group1", "group2", "p")
      plot1 <- ((ggplot(genomic.plot.data, aes(x, y)) + geom_bar(stat='identity') + ggtitle(paste(gene, "LOF Alterations", sep = " ")) 
             + xlab("Genetic Background") + ylab("Alteration Frequency") + ylim(0,1)
             + stat_pvalue_manual(stats, label = "p", tip.length = 0.01, 
                                  y.position = 0.9)))
    }
    if (summary.full$Alteration[i] == "GOF"){
      primary.WT <- primary.WT.GOF.frequency$`GOF alteration rate: Primary WT patients`[which(primary.WT.GOF.frequency$Hugo_Symbol == gene)]
      primary.coloss <- primary.coloss.GOF.frequency$`GOF alteration rate: Primary coloss patients`[which(primary.coloss.GOF.frequency$Hugo_Symbol == gene)]
      met.WT <- met.WT.GOF.frequency$`GOF alteration rate: Metastatic WT patients`[which(met.coloss.GOF.frequency$Hugo_Symbol == gene)]
      met.coloss <- met.coloss.GOF.frequency$`GOF alteration rate: Metastatic coloss patients`[which(met.coloss.GOF.frequency$Hugo_Symbol == gene)]
      x <- c("Primary WT", "Primary Subtype", "Metastatic WT", "Metastatic Subtype")
      y <- c(primary.WT, primary.coloss, met.WT, met.coloss)
      genomic.plot.data <- data.frame(x, y)
      genomic.plot.data$x <- factor(genomic.plot.data$x, levels = genomic.plot.data$x)
      stats <- t(c("Freq", "Primary Subtype", "Metastatic Subtype",
                   paste("p=", formatC(summary.full$`Alteration Frequency p-val`[i], format = "e", digits = 2), sep = "")))
      stats <- as.data.frame(stats)
      colnames(stats) <- c(".y.", "group1", "group2", "p")
      plot1 <- ((ggplot(genomic.plot.data, aes(x, y)) + geom_bar(stat='identity') + ggtitle(paste(gene, "GOF Alterations", sep = " ")) 
             + xlab("Genetic Background") + ylab("Alteration Frequency") + ylim(0,1)
             + stat_pvalue_manual(stats, label = "p", tip.length = 0.01, 
                                  y.position = 0.9)))
    }
    #---------------------------------------CONCORDANT DGE--------------------------
    DGE.boxplot <- data.frame()
    DGE.boxplot[i,1] <- summary.full$`Hugo Symbol`[i]
    # Redefine TCGA coloss and met coloss variables
    primary.coloss <- subtype_subset(subtype_subset(primary, c("MAP3K7", "loss")), c("CHD1", "loss"))
    met.coloss <- subtype_subset(subtype_subset(met, c("MAP3K7", "loss")), c("CHD1", "loss"))
    if (summary.full$Alteration[i] == "LOF"){
      DGE.boxplot[i,2] <- "LOF"
      primary.altered <- subtype_subset(TCGA.coloss, c(summary.full$`Hugo Symbol`[i], "loss"))
      primary.unaltered <- subtype_subset(TCGA.coloss, c(summary.full$`Hugo Symbol`[i], "WT"))
      primary.altered.name <- paste(summary.full$`Hugo Symbol`[i], "Altered", sep = " ")
      primary.unaltered.name <- paste(summary.full$`Hugo Symbol`[i], "WT", sep = " ")
      met.altered <- subtype_subset(met.coloss, c(summary.full$`Hugo Symbol`[i], "loss"))
      met.unaltered <- subtype_subset(met.coloss, c(summary.full$`Hugo Symbol`[i], "WT"))
      met.altered.name <- paste(summary.full$`Hugo Symbol`[i], "Altered", sep = " ")
      met.unaltered.name <- paste(summary.full$`Hugo Symbol`[i], "WT", sep = " ")
    }
    if (summary.full$Alteration[i] == "GOF"){
      DGE.boxplot[i,2] <- "GOF"
      primary.altered <- subtype_subset(TCGA.coloss, c(summary.full$`Hugo Symbol`[i], "gain"))
      primary.unaltered <- subtype_subset(TCGA.coloss, c(summary.full$`Hugo Symbol`[i], "WT"))
      primary.altered.name <- paste(summary.full$`Hugo Symbol`[i], "Altered", sep = " ")
      primary.unaltered.name <- paste(summary.full$`Hugo Symbol`[i], "WT", sep = " ")
      met.altered <- subtype_subset(met.coloss, c(summary.full$`Hugo Symbol`[i], "gain"))
      met.unaltered <- subtype_subset(met.coloss, c(summary.full$`Hugo Symbol`[i], "WT"))
      met.altered.name <- paste(summary.full$`Hugo Symbol`[i], "Altered", sep = " ")
      met.unaltered.name <- paste(summary.full$`Hugo Symbol`[i], "WT", sep = " ")
    }
    primary.unaltered.patients <- colnames(primary.unaltered)[2:ncol(primary.unaltered)]
    primary.altered.patients <- colnames(primary.altered)[2:ncol(primary.altered)]
    met.unaltered.patients <- colnames(met.unaltered)[2:ncol(met.unaltered)]
    met.altered.patients <- colnames(met.altered)[2:ncol(met.altered)]
    # Subset TCGA mRNA data  
    primary.altered.GEX.patients <- intersect(primary.altered.patients, colnames(TCGA.mRNA)[2:ncol(TCGA.mRNA)])
    primary.unaltered.GEX.patients <- intersect(primary.unaltered.patients, colnames(TCGA.mRNA)[2:ncol(TCGA.mRNA)])
    index <- c(TRUE)
    for (k in 2:ncol(TCGA.mRNA)){
      if (sum(colnames(TCGA.mRNA)[k] == primary.altered.GEX.patients) == 1){
        index[k] <- TRUE} else {index[k] <- FALSE}
    }
    primary.altered.GEX <- TCGA.mRNA[, index]
    index <- c(TRUE)
    for (k in 2:ncol(TCGA.mRNA)){
      if (sum(colnames(TCGA.mRNA)[k] == primary.unaltered.GEX.patients) == 1){
        index[k] <- TRUE} else {index[k] <- FALSE}
    }
    primary.unaltered.GEX <- TCGA.mRNA[, index]
    # Subset Abida mRNA data
    met.altered.GEX.patients <- intersect(met.altered.patients, colnames(Abida.mRNA)[2:ncol(Abida.mRNA)])
    met.unaltered.GEX.patients <- intersect(met.unaltered.patients, colnames(Abida.mRNA)[2:ncol(Abida.mRNA)])
    index <- c(TRUE)
    for (k in 2:ncol(Abida.mRNA)){
      if (sum(colnames(Abida.mRNA)[k] == met.altered.GEX.patients) == 1){
        index[k] <- TRUE} else {index[k] <- FALSE}
    }
    met.altered.GEX <- Abida.mRNA[, index]
    index <- c(TRUE)
    for (k in 2:ncol(Abida.mRNA)){
      if (sum(colnames(Abida.mRNA)[k] == met.unaltered.GEX.patients) == 1){
        index[k] <- TRUE} else {index[k] <- FALSE}
    }
    met.unaltered.GEX <- Abida.mRNA[, index]
    # Combine primary data into single data frame
    index <- which(summary.full$`Hugo Symbol`[i] == primary.altered.GEX$Hugo_Symbol)[1]
    primary.altered.vec <- as.numeric(primary.altered.GEX[index,2:ncol(primary.altered.GEX)])
    primary.unaltered.vec <- as.numeric(primary.unaltered.GEX[index,2:ncol(primary.unaltered.GEX)])
    data.primary <- data.frame()
    data.primary[1:length(primary.unaltered.vec),1] <- primary.unaltered.name
    data.primary[1:length(primary.unaltered.vec),2] <- primary.unaltered.vec
    data.primary[((length(primary.unaltered.vec)+1):(length(primary.unaltered.vec)+length(primary.altered.vec))),1] <- primary.altered.name
    data.primary[((length(primary.unaltered.vec)+1):(length(primary.unaltered.vec)+length(primary.altered.vec))),2] <- primary.altered.vec
    data.primary[,3] <- "Primary Subtype"
    colnames(data.primary) <- c("Alteration Status",
                                paste(gene, "Log2 Normalized Expression", sep = " "),
                                "Background")
    data.primary$`Alteration Status` <- factor(data.primary$`Alteration Status`, 
                                               levels = c((paste(gene, "WT", sep = " ")), (paste(gene, "Altered", sep = " "))))
    # Combine metastatic data into single data frame
    index <- which(summary.full$`Hugo Symbol`[i] == met.altered.GEX$Hugo_Symbol)[1]
    met.altered.vec <- as.numeric(met.altered.GEX[index,2:ncol(met.altered.GEX)])
    met.unaltered.vec <- as.numeric(met.unaltered.GEX[index,2:ncol(met.unaltered.GEX)])
    data.met <- data.frame()
    data.met[1:length(met.unaltered.vec),1] <- met.unaltered.name
    data.met[1:length(met.unaltered.vec),2] <- met.unaltered.vec
    data.met[((length(met.unaltered.vec)+1):(length(met.unaltered.vec)+length(met.altered.vec))),1] <- met.altered.name
    data.met[((length(met.unaltered.vec)+1):(length(met.unaltered.vec)+length(met.altered.vec))),2] <- met.altered.vec
    data.met[,3] <- "Metastatic Subtype"
    colnames(data.met) <- c("Alteration Status",
                            paste(gene, "Log2 Normalized Expression", sep = " "),
                            "Background")
    data.met$`Alteration Status` <- factor(data.met$`Alteration Status`, 
                                           levels = c((paste(gene, "WT", sep = " ")), (paste(gene, "Altered", sep = " "))))
    # Bind primary and metastatic data
    data <- rbind(data.primary, data.met)
    data$Background <- factor(data$Background,
                              levels = c("Primary Subtype", "Metastatic Subtype"))
    # Make boxplot
    plot2 <- ((ggplot(data, aes(`Alteration Status`, data[,2], fill = `Alteration Status`))
           + geom_boxplot() + geom_jitter() 
           + ylab(paste(gene, "log2(Normalized mRNA Expression)")) + ylim(0,14)
           + facet_wrap(~Background)) + stat_compare_means(method = "t.test", vjust = 0.5, hjust = -0.2, label = "p"))
    #---------------------------------------TUMOR GRADE-----------------------------
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
      index2 <- which(gene == TCGA.mRNA$Hugo_Symbol)
      data[i,3] <- TCGA.mRNA[index2,index]
      }
    colnames(data) <- c("Score", "Risk Group", "Expression")
    # Order dataset by gene expression. Subset into high and low expressors.
    data <- data[order(data[,2], decreasing = TRUE),]
    # Stratify patients into n=2 quantiles.
    quantiles <- quantile(data$Expression, c(0.50))
    high.expressing <- data[data$Expression > quantiles[1],]
    high.expressing[,4] <- "High Expressers"
    colnames(high.expressing)[4] <- "mRNA Expression"
    low.expressing <- data[data$Expression <= quantiles[1],]
    low.expressing[,4] <- "Low Expressers"
    colnames(low.expressing)[4] <- "mRNA Expression"
    data <- rbind(low.expressing, high.expressing)
    data$`mRNA Expression` <- factor(data$`mRNA Expression`, levels = c("Low Expressers", "High Expressers"))
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
    stats <- t(c("Count", "Low Expressers", "High Expressers",
                 paste("p=", formatC(pval, format = "e", digits = 2), sep = "")))
    stats <- as.data.frame(stats)
    colnames(stats) <- c(".y.", "group1", "group2", "p")
    plot3 <- (ggplot(data, aes(`mRNA Expression`)) + geom_bar(aes(fill=`Risk Group`)) 
              + labs(title = "Primary Subtype Patients", x = paste(gene, "mRNA Expression"), y = "Count")
              + ylim(0,25) + stat_pvalue_manual(stats, label = "p", y.position = 23))
    #---------------------------------------COMBINE ALL PLOTS-------------------
    grid.arrange(plot1, plot2, plot1, plot2)
    }
  dev.off()
  }