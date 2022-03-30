# Subset CNA dataset based on copy number loss, WT/gain genotypes
# x = c("Hugo_Symbol", "type of loss or gain")
subtype_subset <- function (dataset, x) {
  if (x[2] == "homdel") {
    samples = t(subset(dataset, Hugo_Symbol == x[1])[2:ncol(dataset)])
    index.loss = which(samples == -2) + 1
    subset = dataset[, c(1,index.loss)]
    return(subset)
  }
  if (x[2] == "hetloss") {
    samples = t(subset(dataset, Hugo_Symbol == x[1])[2:ncol(dataset)])
    index.loss = which(samples == -1) + 1
    subset = dataset[, c(1,index.loss)]
    return(subset)
  }
  if (x[2] == "loss"){
    samples = t(subset(dataset, Hugo_Symbol == x[1]))[2:ncol(dataset)]
    index.wtloss = which(samples < 0) + 1
    subset = dataset[, c(1,index.wtloss)]
    return(subset)
  }
  if (x[2] == "WT"){
    samples = t(subset(dataset, Hugo_Symbol == x[1]))[2:ncol(dataset)]
    index.wt = which(samples == 0) + 1
    subset = dataset[, c(1,index.wt)]
    return(subset)
  }
  if (x[2] == "gain"){
    samples = t(subset(dataset, Hugo_Symbol == x[1]))[2:ncol(dataset)]
    index.gain = which(samples > 0) + 1
    subset = dataset[, c(1,index.gain)]
    return(subset)
  }
  if (x[2] == "WT/gain"){
    samples = t(subset(dataset, Hugo_Symbol == x[1]))[2:ncol(dataset)]
    index.wtgain = which(samples >= 0) + 1
    subset = dataset[, c(1,index.wtgain)]
    return(subset)
  }
}