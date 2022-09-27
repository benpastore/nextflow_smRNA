#!/usr/bin/env Rscript

options(warn = -1)
library(dplyr)
library(tidyr)

parser = argparse::ArgumentParser()
parser$add_argument("-i", "--input", type = "character", required = TRUE)
parser$add_argument("-c", "--condition", type = "character", required = TRUE)
parser$add_argument("-id",type = "character", required = FALSE)

args = parser$parse_args()

condition = read.delim(args$condition, sep = "\t", header = F)
names(condition) = c("sample", "condition")
condition$condition = as.character(condition$condition)
condition = condition %>% mutate(condition = ifelse(condition == "*", sample, condition))

df = read.delim(args$input, sep = "\t")

if (is.null(args$id)){
    id = c("gene_name", "seq_id", "locus_id", "class", "feature", "biotype")
} else {
    id = unlist(strsplit(args$id,","))
}

group_cols = c(id, "condition")

long = df %>%
    pivot_longer(!id, names_to = "sample", values_to = "count") %>%
    left_join(condition, by = "sample") %>%
    mutate(condition = ifelse(is.na(condition), sample, condition)) %>%
    group_by_at(group_cols) %>%
    summarise(count = mean(count)) %>%
    pivot_wider(names_from = condition, values_from = count)

outname = gsub(".tsv", ".averaged.tsv", args$input)
write.table(long, outname, sep = "\t", col.names = T, row.names = F, quote = F)