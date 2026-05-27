#!/bin/Rscript

library(data.table)
library("scales")

data <- fread("qv_plot_input.tsv", header=T, sep="\t")
data$AsmGroup <- factor(data$AsmGroup, levels = c("HPRC_r2", "CuratedY"))

summary(data[data$AsmGroup=="HPRC_r2", ], digits = 3)
summary(data[data$AsmGroup=="CuratedY", ], digits = 3)

standard_deviation <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  sqrt(sum((x - mean_x)^2) / (n - 1))
}

print("Standard Deviation of HPRC_r2 QV:")
standard_deviation(data[data$AsmGroup=="HPRC_r2", ]$QV)

print("Standard Deviation of CuratedY QV:")
standard_deviation(data[data$AsmGroup=="CuratedY", ]$QV)
