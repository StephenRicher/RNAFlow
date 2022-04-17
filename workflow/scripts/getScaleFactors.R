#!/usr/bin/env Rscript

require(data.table)
require(edgeR)

## turn off exponential notation to avoid rounding errors
options(scipen=999) 

args = commandArgs(trailingOnly=TRUE)

i = 0
for (file in args) {
  sample = substring(file, 15, nchar(file)-18)
  ## Load count matrix, removing all columns that we are not interested in:
  sample.counts <- fread(file, data.table=FALSE, header=TRUE, skip=c(1))
  sample.counts <- sample.counts[,7:ncol(sample.counts)]
  if (i == 0) {
    raw.counts = as.data.frame(sample.counts)
    colnames(raw.counts) = c(sample)
  } else{
    raw.counts[sample] = sample.counts
  }
  i = i + 1
}

## effective normalization factors are the product 
## of TMM factor and library size for each sample:
norm.factor   <- calcNormFactors(object = raw.counts, method = c("TMM"))
lib.size      <- colSums(raw.counts)
final.factor  <- norm.factor * lib.size

## as one typically scales reads "per million" we divide this factor by 10^6
## we also need the reciprocal value when using bamCoverage later, you'll see why,
## see *comment below
perMillion.factor <- (final.factor / 1000000)^-1

write.table(x = data.frame(Sample     = names(perMillion.factor),
                           NormFactor = perMillion.factor),
            file = stdout(), sep="\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)
