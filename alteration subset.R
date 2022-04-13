# Subsets a binary alteration matrix by alteration of a given gene
# Dataset = binary alteration matrix
# x = c("Hugo_Symbol", "unaltered or altered")
alteration_subset <- function (dataset, x) {
  if (x[2] == "altered") {
    samples = t(subset(dataset, Hugo_Symbol == x[1])[2:ncol(dataset)])
    index.loss = which(samples > 0) + 1
    subset = dataset[, c(1,index.loss)]
    return(subset)
  }
  if (x[2] == "unaltered") {
    samples = t(subset(dataset, Hugo_Symbol == x[1])[2:ncol(dataset)])
    index.loss = which(samples == 0) + 1
    subset = dataset[, c(1,index.loss)]
    return(subset)
  }
}