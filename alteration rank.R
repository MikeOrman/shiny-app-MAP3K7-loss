# Dataset = TRUE/FALSE GOF or LOF alteration matrix.
#   Col 1: "Hugo_Symbol"
#   Col 2 - End: Sample names
#   Rows 1 - End: Feature names (Hugo Symbols for genes)
alterationRank <- function(dataset){
  gene.names = dataset[,1]
  alteration.freq = data.frame()
  logical = (dataset > 0)
  for (i in 1:nrow(logical)) {
    alteration.freq[i,1] = gene.names[i]
    alteration.freq[i,2] = sum(logical[i,2:ncol(logical)], na.rm = TRUE)
    alteration.freq[i,3] = (ncol(logical)-1) - sum(logical[i,2:ncol(logical)], na.rm = TRUE)
    alteration.freq[i,4] = sum(logical[i,2:ncol(logical)], na.rm = TRUE)/(ncol(logical)-1)
  }
  ordered = order(alteration.freq[,4])
  x = c()
  y = c()
  z = c()
  a = c()
  b = c()
  for (i in 1:length(ordered)) {
    x[i] = i
    y[i] = alteration.freq[ordered[i], 4]
    z[i] = gene.names[ordered[i]]
    a[i] = alteration.freq[ordered[i], 2]
    b[i] = alteration.freq[ordered[i], 3]
  }
  output = data.frame(z, x, a, b, y)
  colnames(output) = c("Hugo_Symbol", "Gene Rank", "total altered", 
                       "total unaltered", "alteration freq")
  return(output)
}