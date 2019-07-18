
setwd("~/Projects/dke")
library(ggplot2)

df = read.csv('output.csv')
metadata =read.csv('metadata.csv')

data = merge(df, metadata, by.x="id", by.y="runid")


ggplot(data, aes(gen, gst, color=factor(n_alleles), group=id)) + facet_wrap(eff_pop ~ n_pops) + geom_line(alpha=0.5)

ggplot(data, aes(gen, jostd, color=factor(n_alleles), group=id)) + facet_wrap(eff_pop ~ n_pops) + geom_line(alpha=0.5)
 