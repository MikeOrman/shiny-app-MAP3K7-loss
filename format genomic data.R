#--------------------------LOAD cBIO CNA DATA-----------------------------------
# TCGA
TCGA <- read.table("TCGA_data.txt", check.names = FALSE, sep = "\t")
TCGA <- as.data.frame(TCGA)
# Taylor
Taylor <- read.table("Taylor_data.txt", check.names = FALSE, sep = "\t")
Taylor <- as.data.frame(Taylor)
Taylor <- na.omit(Taylor)
# Baca
Baca <- read.table("Baca_data.txt", check.names = FALSE, sep = "\t")
Baca <- as.data.frame(Baca)
# Abida
Abida <- read.table("Abida_data.txt", check.names = FALSE, sep = "\t")
Abida <- as.data.frame(Abida)
#--------------------------COMBINED PRIMARY CNA DATASETS------------------------
met <- Abida
# Combine primary sets
primary <- merge.data.frame(TCGA, Taylor)
primary <- primary[duplicated(primary$Hugo_Symbol)==FALSE,]
primary <- merge.data.frame(primary, Baca)
primary <- primary[duplicated(primary$Hugo_Symbol)==FALSE,]
# Filter genes missing in both primary and metastatic settings
genes <- intersect(primary$Hugo_Symbol, met$Hugo_Symbol)
# Create filtering index for primary dataset
primary.index <- c()
for (i in 1:nrow(primary)){
  if (sum(primary$Hugo_Symbol[i] == genes) == 1)
  {primary.index[i] = TRUE}
  if (sum(primary$Hugo_Symbol[i] == genes) == 0)
  {primary.index[i] = FALSE}
}
# Create filtering index for metastatic dataset
met.index <- c()
for (i in 1:nrow(met)){
  if (sum(met$Hugo_Symbol[i] == genes) == 1)
  {met.index[i] = TRUE}
  if (sum(met$Hugo_Symbol[i] == genes) == 0)
  {met.index[i] = FALSE}
}
# Filter primary and metastatic genes. The number of rows for primary and metastatic datasets should now be equal.
primary <- primary[primary.index,]
primary <- primary[duplicated(primary)==FALSE,]
met <- met[met.index,]
met <- met[duplicated(met)==FALSE,]
#--------------------------TRANSFORM CNA DATA TO ALTERED BINARY MATRICES--------
# Primary alteration matrix
primary.altered.data <- primary[,c(2:ncol(primary))]
primary.altered.data[primary[,c(2:ncol(primary))] < 0] <- 1
primary.altered.data[primary[,c(2:ncol(primary))] > 0] <- 1
primary.altered <- primary
primary.altered[,c(2:ncol(primary))] <- primary.altered.data
# Metastatic alteration matrix
met.altered.data <- met[,c(2:ncol(met))]
met.altered.data[met[,c(2:ncol(met))] < 0] <- 1
met.altered.data[met[,c(2:ncol(met))] > 0] <- 1
met.altered <- met
met.altered[,c(2:ncol(met))] <- met.altered.data
#--------------------------TRANSFORM CNA DATA TO GOF / LOF  BINARY MATRICES-----
# Primary LOF binary matrix
primary.LOF.data <- primary[,c(2:ncol(primary))]
primary.LOF.data[primary[,c(2:ncol(primary))] < 0] <- 1
primary.LOF.data[primary[,c(2:ncol(primary))] >= 0] <- 0
primary.LOF <- primary
primary.LOF[,c(2:ncol(primary))] <- primary.LOF.data
# Primary GOF binary matrix
primary.GOF.data <- primary[,c(2:ncol(primary))]
primary.GOF.data[primary[,c(2:ncol(primary))] > 0] <- 1
primary.GOF.data[primary[,c(2:ncol(primary))] <= 0] <- 0
primary.GOF <- primary
primary.GOF[,c(2:ncol(primary))] <- primary.GOF.data
# Met LOF binrary matrix
met.LOF.data <- met[,c(2:ncol(met))]
met.LOF.data[met[,c(2:ncol(met))] < 0] <- 1
met.LOF.data[met[,c(2:ncol(met))] >= 0] <- 0
met.LOF <- met
met.LOF[,c(2:ncol(met))] <- met.LOF.data
# Met GOF binary matrix
met.GOF.data <- met[,c(2:ncol(met))]
met.GOF.data[met[,c(2:ncol(met))] > 0] <- 1
met.GOF.data[met[,c(2:ncol(met))] <= 0] <- 0
met.GOF <- met
met.GOF[,c(2:ncol(met))] <- met.GOF.data
# TCGA LOF binary matrix
TCGA.LOF.data <- TCGA[,c(2:ncol(TCGA))]
TCGA.LOF.data[TCGA[,c(2:ncol(TCGA))] < 0] <- 1
TCGA.LOF.data[TCGA[,c(2:ncol(TCGA))] >= 0] <- 0
TCGA.LOF <- TCGA
TCGA.LOF[,c(2:ncol(TCGA))] <- TCGA.LOF.data
# TCGA GOF binary matrix
TCGA.GOF.data <- TCGA[,c(2:ncol(TCGA))]
TCGA.GOF.data[TCGA[,c(2:ncol(TCGA))] > 0] <- 1
TCGA.GOF.data[TCGA[,c(2:ncol(TCGA))] <= 0] <- 0
TCGA.GOF <- TCGA
TCGA.GOF[,c(2:ncol(TCGA))] <- TCGA.GOF.data
#--------------------------ADD DELETERIOUS MUTATIONS TO ALTERATION MATRIX-------
mutations <- read.table("Mutations.txt", sep = "\t", header = TRUE, check.names = FALSE)
# Add mutations to primary alteration matrix
# i loop: iterate through primary samples
for (i in 2:ncol(primary.altered)){
  sample <- colnames(primary.altered)[i]
  # Find if sample name exists in mutation data
  sample.index <- which(sample == mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.integer(sample.index))
    sample.mutations <- mutations[sample.index,]
  # k loop: iterate through the sample's mutations.
  for (k in 1:nrow(sample.mutations)){
    gene <- sample.mutations$`Hugo Symbol`[k]
    gene.index <- which(gene == primary.altered$Hugo_Symbol)
    # If the gene has a mutation present, then add its LOF prediction to the alteration count.
    if (is.integer(gene.index)){
      primary.altered[gene.index, i] <- primary.altered[gene.index, i] + sample.mutations$`LOF Pred`[k]
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
# Add mutations to met alteration matrix
# i loop: iterate through met samples
for (i in 2:ncol(met.altered)){
  sample <- colnames(met.altered)[i]
  # Find if sample name exists in mutation data
  sample.index <- which(sample == mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.integer(sample.index))
    sample.mutations <- mutations[sample.index,]
  # k loop: iterate through the sample's mutations.
  for (k in 1:nrow(sample.mutations)){
    gene <- sample.mutations$`Hugo Symbol`[k]
    gene.index <- which(gene == met.altered$Hugo_Symbol)
    # If the gene has a mutation present, then add its LOF prediction to the LOF alteration count.
    if (is.integer(gene.index)){
      met.altered[gene.index, i] <- met.altered[gene.index, i] + sample.mutations$`LOF Pred`[k]
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
#--------------------------ADD DELETERIOUS MUTATIONS TO TCGA ALTERATION MATRIX----
TCGA.mutation <- read.table("TCGA formatted mutations.txt", sep = "\t", header = TRUE, check.names = FALSE)
# Primary alteration matrix
TCGA.altered.data <- TCGA[,c(2:ncol(TCGA))]
TCGA.altered.data[TCGA[,c(2:ncol(TCGA))] < 0] <- 1
TCGA.altered.data[TCGA[,c(2:ncol(TCGA))] > 0] <- 1
TCGA.altered <- TCGA
TCGA.altered[,c(2:ncol(TCGA))] <- TCGA.altered.data
# Add TCGA mutations to TCGA alteration matrix
# i loop: iterate through TCGA primary samples
for (i in 2:ncol(TCGA.altered)){
  sample <- colnames(TCGA.altered)[i]
  # Find if sample name exists in mutation data
  sample.index <- which(sample == TCGA.mutation$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.integer(sample.index))
    sample.TCGA.mutation <- TCGA.mutation[sample.index,]
  # k loop: iterate through the sample's mutations.
  for (k in 1:nrow(sample.TCGA.mutation)){
    gene <- sample.TCGA.mutation$`Hugo Symbol`[k]
    gene.index <- which(gene == TCGA.altered$Hugo_Symbol)
    # If the gene has a mutation present, then add its LOF prediction to the alteration count.
    if (is.integer(gene.index)){
      TCGA.altered[gene.index, i] <- TCGA.altered[gene.index, i] + sample.TCGA.mutation$`LOF Pred`[k]
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
#--------------------------ADD DELETERIOUS MUTATIONS TO LOF MATRIX--------------
mutations <- read.table("Mutations.txt", sep = "\t", header = TRUE, check.names = FALSE)
# Add mutations to primary LOF alteration matrix
# i loop: iterate through primary samples
for (i in 2:ncol(primary.LOF)){
  sample <- colnames(primary.LOF)[i]
  # Find if sample name exists in mutation data
  sample.index <- which(sample == mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.integer(sample.index))
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- which(gene == primary.LOF$Hugo_Symbol)
      # If the gene has a mutation present, then add its LOF prediction to the LOF alteration count.
      if (is.integer(gene.index)){
        primary.LOF[gene.index, i] <- primary.LOF[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  # If gene does not exist in mutation data then move to next sample
}
# Add mutations to met LOF alteration matrix
# i loop: iterate through met samples
for (i in 2:ncol(met.LOF)){
  sample <- colnames(met.LOF)[i]
  # Find if sample name exists in mutation data
  sample.index <- which(sample == mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.integer(sample.index))
    sample.mutations <- mutations[sample.index,]
  # k loop: iterate through the sample's mutations.
  for (k in 1:nrow(sample.mutations)){
    gene <- sample.mutations$`Hugo Symbol`[k]
    gene.index <- which(gene == met.LOF$Hugo_Symbol)
    # If the gene has a mutation present, then add its LOF prediction to the LOF alteration count.
    if (is.integer(gene.index)){
      met.LOF[gene.index, i] <- met.LOF[gene.index, i] + sample.mutations$`LOF Pred`[k]
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
#--------------------------WRITE FORMATTED GOF / LOF ALTERATION TABLES----------
write.table(primary.altered, "primary alteration.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(primary.LOF, "primary LOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(primary.GOF, "primary GOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(met.altered, "met alteration.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(met.LOF, "met LOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(met.GOF, "met GOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(TCGA.altered, "TCGA alteration.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(TCGA.LOF, "TCGA LOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(TCGA.GOF, "TCGA GOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)