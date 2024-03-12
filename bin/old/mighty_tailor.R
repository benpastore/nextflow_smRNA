#!/usr/bin/env Rscript
options(warn = -1)
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE,repos='http://cran.us.r-project.org')
    if(!require(x,character.only = TRUE)) stop ("Failed to install the package. Please check the internet access or update your R if it is too old.")
  }
}


pkgTest('dplyr')
pkgTest('ggplot2')
pkgTest('lemon')
pkgTest('tidyverse')
pkgTest('scales')
pkgTest('RColorBrewer')
pkgTest('ggpubr')
pkgTest('ggrepel')
pkgTest('GGally')
pkgTest('ggpmisc')
pkgTest('rstatix')
pkgTest('argparse')
pkgTest('this.path')
pkgTest('stringr')

library(dplyr)
library(ggplot2)
library(lemon)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(GGally)
library(ggpmisc)
library(rstatix)
library(this.path)
library(argparse)
library(stringr)

source( paste0(this.dir(),"/mighty.R") )

parser = argparse::ArgumentParser()

parser$add_argument("-t", "--tailor", type = "character", required = TRUE)
parser$add_argument("-u", "--untailed", type = "character", required = TRUE)
parser$add_argument("-c", "--condition", type = "character", required = TRUE)
parser$add_argument("-s", "--samples", type = "character", required = TRUE)
parser$add_argument("-o", "--outname", type = "character", required = TRUE)
parser$add_argument("-f", "--feature", type = "character", required = TRUE)

args = parser$parse_args()

untail = args$untailed
untailed = read.delim(untail)

condition_file = args$condition
condition = read.delim(condition_file, header = F)
names(condition) = c("sample", "condition")
condition$sample = as.character(condition$sample)

tail = args$tailor
tailed = read.delim(tail)
tailed$tail = as.character(tailed$tail)

samps = unlist(strsplit(args$samples,","))

p = plot_tails(tailed %>% filter(feature == args$feature), untailed %>% filter(feature == args$feature), condition, samps)

ggsave(p, filename = paste0(args$outname), dpi = 300, height = 6, width = 6)







