library(ggplot2)
library(ggrepel)
source("manhattan.R")
library(ggpubr)
library(rstatix)
library(gridExtra)
library(survival)
library(survminer)
# Note: Run this after running all on "shiny app V4 - MAP3K7 altered.Rmd"
query.gene <- "MAP3K7"
query.gene.alteration <- "MAP3K7 LOF"
#---------------------------------------LOAD WORKSPACE--------------------------
load("~/Github/Projects/shiny app MAP3K7 loss/shiny app V4 - MAP3K7 loss.RData")
#---------------------------------------RANKED LIST-----------------------------
gene.rank.dotplot <- function(input){
  pdf("Gene Rank.pdf")
  input$`Combined FDR` <- -log10(input$`Combined FDR`)
  # Mark plot with FDR cutoff corresponding to top 5% hits
  FDR <- quantile(input$`Combined FDR`, c(0.95))
  h = round(FDR, 2)
  for (i in 1:nrow(input)){
    gene <- input$`Hugo Symbol`[i]
    gene.layer <- subset(input,  `Hugo Symbol` == gene)
    print((ggplot(input, aes(`Gene Rank`, `Combined FDR`, label = `Hugo Symbol`))
      + geom_point()
      + xlim(0, nrow(input))
      + ylab("-log10 Combined FDR")
      + geom_point(data = gene.layer, colour="red")
      + geom_text_repel(data = subset(input, `Hugo Symbol` == gene), color = "red", nudge_x = 0, nudge_y = .3)
      + geom_hline(yintercept = h , linetype = "solid", color = "red", size = 0.8)
      + geom_text(aes(0,h,label = paste("FDR=", formatC(10^(-FDR), format = "e", digits = 2), sep = "")),
                      nudge_x = .75*nrow(input), nudge_y = 0.05*h)))
    }
  dev.off()
}
#---------------------------------------KARYOPLOT-------------------------------
karyoplots <- function(input){
  # Create manhattan object
  mp.input <- input[,c(1,4)]
  # Clear significance. This will plot genes at the same height, regardless of P-value.
  mp.input[,2] <- 1
  mp <- manhattan.object(mp.input)
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
genomic.plots <- function(input){
  pdf(file = "Alteration bar plots.pdf")
  for (i in 1:nrow(input)){
    gene <- input$`Hugo Symbol`[i]
    if (input$Alteration[i] == "LOF"){
      primary.WT <- primary.WT.LOF.frequency$`LOF alteration rate: Primary WT patients`[which(primary.WT.LOF.frequency$Hugo_Symbol == gene)]
      primary.ST <- primary.ST.LOF.frequency$`LOF alteration rate: Primary ST patients`[which(primary.ST.LOF.frequency$Hugo_Symbol == gene)]
      met.WT <- met.WT.LOF.frequency$`LOF alteration rate: Metastatic WT patients`[which(met.WT.LOF.frequency$Hugo_Symbol == gene)]
      met.ST <- met.ST.LOF.frequency$`LOF alteration rate: Metastatic ST patients`[which(met.ST.LOF.frequency$Hugo_Symbol == gene)]
      x <- c(paste("Primary", query.gene, "WT"), 
             paste("Primary", query.gene, "LOF"), 
             paste("Metastatic", query.gene, "WT"), 
             paste("Metastatic", query.gene, "LOF"))
      y <- c(primary.WT, primary.ST, met.WT, met.ST)
      genomic.plot.data <- data.frame(x, y)
      genomic.plot.data$x <- factor(genomic.plot.data$x, levels = x)
      stats <- t(c("Alteration Frequency", paste("Primary", query.gene, "LOF"), paste("Metastatic", query.gene, "LOF"),
                   paste("p=", formatC(input$`Subtype Metastatic Enrichment pval`[i], format = "e", digits = 2), sep = "")))
      stats <- as.data.frame(stats)
      colnames(stats) <- c(".y.", "group1", "group2", "p")
      print((ggplot(genomic.plot.data, aes(x, y))
             + geom_bar(stat='identity', fill = c("black", "blue", "black", "blue"))
             + ggtitle(paste(gene, "LOF Alterations", sep = " "))
             + theme(axis.title.x = element_blank())
             + xlab("Patient Background")
             + ylab("Alteration Frequency") + ylim(0,1)
             + stat_pvalue_manual(stats, label = "p", tip.length = 0.01, y.position = max(y) + .0625)))
      }
    if (input$Alteration[i] == "GOF"){
      primary.WT <- primary.WT.GOF.frequency$`GOF alteration rate: Primary WT patients`[which(primary.WT.GOF.frequency$Hugo_Symbol == gene)]
      primary.ST <- primary.ST.GOF.frequency$`GOF alteration rate: Primary ST patients`[which(primary.ST.GOF.frequency$Hugo_Symbol == gene)]
      met.WT <- met.WT.GOF.frequency$`GOF alteration rate: Metastatic WT patients`[which(met.WT.GOF.frequency$Hugo_Symbol == gene)]
      met.ST <- met.ST.GOF.frequency$`GOF alteration rate: Metastatic ST patients`[which(met.ST.GOF.frequency$Hugo_Symbol == gene)]
      x <- c(paste("Primary", query.gene, "WT"), 
             paste("Primary", query.gene, "GOF"), 
             paste("Metastatic", query.gene, "WT"), 
             paste("Metastatic", query.gene, "GOF"))      
      y <- c(primary.WT, primary.ST, met.WT, met.ST)
      genomic.plot.data <- data.frame(x, y)
      genomic.plot.data$x <- factor(genomic.plot.data$x, levels = x)
      stats <- t(c("Alteration Frequency", paste("Primary", query.gene, "GOF"), paste("Metastatic", query.gene, "GOF"),
                   paste("p=", formatC(input$`Subtype Metastatic Enrichment pval`[i], format = "e", digits = 2), sep = "")))
      stats <- as.data.frame(stats)
      colnames(stats) <- c(".y.", "group1", "group2", "p")
      print((ggplot(genomic.plot.data, aes(x, y)) 
             + geom_bar(stat='identity', fill = c("black", "blue", "black", "blue")) 
             + ggtitle(paste(gene, "GOF Alterations", sep = " ")) 
             + theme(axis.title.x = element_blank())
             + xlab("Patient Background")
             + ylab("Alteration Frequency") + ylim(0,1)
             + stat_pvalue_manual(stats, label = "p", tip.length = 0.01, y.position = max(y) + .0625)))
    }
  }
  dev.off()
}
#---------------------------------------CONCORDANT DGE--------------------------
concordant.DGE.boxplot <- function(input){
  DGE.boxplot <- data.frame()
  pdf(file = "Concordant gene expression plots.PDF")
  for (i in 1:nrow(input)){
    gene <- input$`Hugo Symbol`[i]
    DGE.boxplot[i,1] <- input$`Hugo Symbol`[i]
    #---------------------------SUBSET PATIENTS-------------------------------------
    if (input$Alteration[i] == "LOF"){
      DGE.boxplot[i,2] <- "LOF"
      primary.altered <- alteration_subset(TCGA.LOF, c(input$`Hugo Symbol`[i], "altered"))
      primary.unaltered <- alteration_subset(TCGA.LOF, c(input$`Hugo Symbol`[i], "unaltered"))
      primary.altered.name <- paste(input$`Hugo Symbol`[i], "LOF", sep = " ")
      primary.unaltered.name <- paste(input$`Hugo Symbol`[i], "WT", sep = " ")
      met.altered <- alteration_subset(met.LOF, c(input$`Hugo Symbol`[i], "altered"))
      met.unaltered <- alteration_subset(met.LOF, c(input$`Hugo Symbol`[i], "unaltered"))
      met.altered.name <- paste(input$`Hugo Symbol`[i], "LOF", sep = " ")
      met.unaltered.name <- paste(input$`Hugo Symbol`[i], "WT", sep = " ")
    }
    if (input$Alteration[i] == "GOF"){
      DGE.boxplot[i,2] <- "GOF"
      primary.altered <- alteration_subset(TCGA.ST, c(input$`Hugo Symbol`[i], "altered"))
      primary.unaltered <- alteration_subset(TCGA.ST, c(input$`Hugo Symbol`[i], "unaltered"))
      primary.altered.name <- paste(input$`Hugo Symbol`[i], "GOF", sep = " ")
      primary.unaltered.name <- paste(input$`Hugo Symbol`[i], "WT", sep = " ")
      met.altered <- alteration_subset(met.ST, c(input$`Hugo Symbol`[i], "altered"))
      met.unaltered <- alteration_subset(met.ST, c(input$`Hugo Symbol`[i], "unaltered"))
      met.altered.name <- paste(input$`Hugo Symbol`[i], "GOF", sep = " ")
      met.unaltered.name <- paste(input$`Hugo Symbol`[i], "WT", sep = " ")
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
      index <- which(input$`Hugo Symbol`[i] == primary.altered.GEX$Hugo_Symbol)[1]
      primary.altered.vec <- as.numeric(primary.altered.GEX[index,2:ncol(primary.altered.GEX)])
      primary.unaltered.vec <- as.numeric(primary.unaltered.GEX[index,2:ncol(primary.unaltered.GEX)])
      data.primary <- data.frame()
      data.primary[1:length(primary.unaltered.vec),1] <- primary.unaltered.name
      data.primary[1:length(primary.unaltered.vec),2] <- primary.unaltered.vec
      data.primary[((length(primary.unaltered.vec)+1):(length(primary.unaltered.vec)+length(primary.altered.vec))),1] <- primary.altered.name
      data.primary[((length(primary.unaltered.vec)+1):(length(primary.unaltered.vec)+length(primary.altered.vec))),2] <- primary.altered.vec
      data.primary[,3] <- paste("Primary", query.gene.alteration)
      colnames(data.primary) <- c("Alteration Status",
                                  paste(gene, "Log2 Normalized Expression", sep = " "),
                                  "Patient Background")
# Combine metastatic data into single data frame
      index <- which(input$`Hugo Symbol`[i] == met.altered.GEX$Hugo_Symbol)[1]
      met.altered.vec <- as.numeric(met.altered.GEX[index,2:ncol(met.altered.GEX)])
      met.unaltered.vec <- as.numeric(met.unaltered.GEX[index,2:ncol(met.unaltered.GEX)])
      data.met <- data.frame()
      data.met[1:length(met.unaltered.vec),1] <- met.unaltered.name
      data.met[1:length(met.unaltered.vec),2] <- met.unaltered.vec
      data.met[((length(met.unaltered.vec)+1):(length(met.unaltered.vec)+length(met.altered.vec))),1] <- met.altered.name
      data.met[((length(met.unaltered.vec)+1):(length(met.unaltered.vec)+length(met.altered.vec))),2] <- met.altered.vec
      data.met[,3] <- paste("Metastatic", query.gene.alteration)
      colnames(data.met) <- c("Alteration Status",
                                    paste(gene, "Log2 Normalized Expression", sep = " "),
                                    "Patient Background")
# Bind primary and metastatic data. Define order of factor levels for plotting.
      data <- rbind(data.primary, data.met)
      if (input$Alteration[i] == "LOF"){
        data$`Alteration Status` <- factor(data$`Alteration Status`, 
                                           levels = c(paste(gene, "WT", sep = " "), 
                                                      paste(gene, "LOF", sep = " ")))
      }
      if (input$Alteration[i] == "GOF"){
        data$`Alteration Status` <- factor(data$`Alteration Status`, 
                                           levels = c(paste(gene, "WT", sep = " "), 
                                                      paste(gene, "GOF", sep = " ")))
      }
      data$`Patient Background` <- factor(data$`Patient Background`,
                                          levels = c(paste("Primary", query.gene.alteration), 
                                                     paste("Metastatic", query.gene.alteration)))
# Format statistics table for plotting
      facet <- c(paste("Primary", query.gene.alteration, sep = " "), 
                 paste("Metastatic", query.gene.alteration, sep = " "))
      y <- paste(gene, "log2(Normalized mRNA Expression)")
      if (input$Alteration[i] == "LOF"){
        group1 <- paste(gene, "WT", sep = " ")
        group2 <- paste(gene, "LOF", sep = " ")
      }
      if (input$Alteration[i] == "GOF"){
        group1 <- paste(gene, "WT", sep = " ")
        group2 <- paste(gene, "GOF", sep = " ")
      }
      FDR <- c(paste("p=", formatC(input$`Primary Expression pval`[input$`Hugo Symbol`==gene], 
                                     format = "e", digits = 2), sep = ""),
               paste("p=", formatC(input$`Metastatic Expression pval`[input$`Hugo Symbol`==gene], 
                                     format = "e", digits = 2), sep = ""))
      xmin <- c(1,1)
      xmax <- c(2,2)
      primary.y.position <- max(data[data$`Patient Background` == paste("Primary", query.gene.alteration, sep = " "),2]) + 1
      met.y.position <- max(data[data$`Patient Background` == paste("Metastatic", query.gene.alteration, sep = " "),2]) + 1
      y.positions <- c(primary.y.position, met.y.position)
      stats <- data.frame(facet, "Alteration Status", y, group1, group2, FDR, xmin, xmax, y.positions)
      colnames(stats) <- c("Patient Background", "Alteration Status", ".y.", "group1", "group2", "p", "xmin", "xmax", "y.position")
# Make boxplot
      print((ggplot(data, aes(`Alteration Status`, data[,2], color = `Alteration Status`))
             + scale_color_manual(values=c("Black", "Blue"))
             + geom_boxplot() + geom_jitter() 
             + theme(axis.title.x = element_blank())
             + ylab(paste(gene, "log2 Normalized mRNA Expression"))
             + ylim(0, max(data[,2])+2)
             + facet_wrap(~`Patient Background`)
             + stat_pvalue_manual(stats, label = "p")))
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
# Plot subtype high-risk enrichment
    sample.names <- colnames(TCGA.ST)[2:ncol(TCGA.ST)]
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
    # Order dataset by gene expression. Subset into high and low expressers.
    data <- data[order(data[,2], decreasing = TRUE),]
    # Stratify patients into n=2 quantiles.
    quantiles <- quantile(data$Expression, c(0.50))
    high.expressing <- data[data$Expression >= quantiles[1],]
    high.expressing[,4] <- "High Expressers"
    colnames(high.expressing)[4] <- "mRNA Expression"
    low.expressing <- data[data$Expression <= quantiles[1],]
    low.expressing[,4] <- "Low Expressers"
    colnames(low.expressing)[4] <- "mRNA Expression"
    data <- rbind(low.expressing, high.expressing)
    data$`mRNA Expression` <- factor(data$`mRNA Expression`, levels = c("Low Expressers", "High Expressers"))
    # Index high-risk enrichment from input table
    FDR.index <- which(gene == input$`Hugo Symbol`)
    FDR <- input$`Subtype High Grade Enrichment pval`[FDR.index]
    stats <- t(c("Count", "Low Expressers", "High Expressers",
                 paste("p=", formatC(FDR, format = "e", digits = 2), sep = "")))
    stats <- as.data.frame(stats)
    colnames(stats) <- c(".y.", "group1", "group2", "FDR")
    print((ggplot(data, aes(`mRNA Expression`)) + geom_bar(aes(fill=`Risk Group`)) 
           + labs(title = "Primary Subtype Patients", x = paste(gene, "mRNA Expression"), y = "Count")
           + ylim(0,(nrow(data)*.60)) + stat_pvalue_manual(stats, label = "FDR", y.position = nrow(data)*.55)))
# Plot WT high-risk enrichment
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
      index2 <- which(gene == TCGA.mRNA$Hugo_Symbol)
      data[i,3] <- TCGA.mRNA[index2,index]
    }
    colnames(data) <- c("Score", "Risk Group", "Expression")
    # Order dataset by gene expression. Subset into high and low expressors.
    data <- data[order(data[,2], decreasing = TRUE),]
    # Stratify patients into n=2 quantiles.
    quantiles <- quantile(data$Expression, c(0.50))
    high.expressing <- data[data$Expression >= quantiles[1],]
    high.expressing[,4] <- "High Expressers"
    colnames(high.expressing)[4] <- "mRNA Expression"
    low.expressing <- data[data$Expression <= quantiles[1],]
    low.expressing[,4] <- "Low Expressers"
    colnames(low.expressing)[4] <- "mRNA Expression"
    data <- rbind(low.expressing, high.expressing)
    data$`mRNA Expression` <- factor(data$`mRNA Expression`, levels = c("Low Expressers", "High Expressers"))
    # Index high-risk enrichment from input table
    FDR <- input$`Wild-Type High Grade Enrichment pval`[FDR.index]
    stats <- t(c("Count", "Low Expressers", "High Expressers",
                 paste("p=", formatC(FDR, format = "e", digits = 2), sep = "")))
    stats <- as.data.frame(stats)
    colnames(stats) <- c(".y.", "group1", "group2", "p")
    print((ggplot(data, aes(`mRNA Expression`)) + geom_bar(aes(fill=`Risk Group`)) 
           + labs(title = "Primary WT Patients", x = paste(gene, "mRNA Expression"), y = "Count")
           + ylim(0,nrow(data)*.60) 
           + stat_pvalue_manual(stats, label = "p", y.position = nrow(data)*.55)))
  }
  dev.off()
}