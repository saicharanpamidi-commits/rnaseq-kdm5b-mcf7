countData <- read.table("results/featurecounts/gene_counts.txt", header=TRUE, skip=1, row.names=1)
print(dim(countData))
print(colSums(countData))
print(sum(rowSums(countData) == 0))
