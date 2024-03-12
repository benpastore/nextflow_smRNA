#!/usr/bin/env Rscript

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


source( paste0(this.dir(),"/mighty.R") )

parser = argparse::ArgumentParser()
parser$add_argument("-f", "--file", type = "character", required = TRUE)
parser$add_argument("-o", "--outname", type = "character", required = TRUE)
parser$add_argument("--sequence_column", type = "character", required = FALSE)
parser$add_argument("--strand_column", type = "character", required = FALSE)
parser$add_argument("--plot_column", type = "character", required = FALSE)
parser$add_argument("--abundance_column", type = "character", required = FALSE)

args = parser$parse_args()

df = read.delim(args$file, sep = "\t")

if (is.null(args$sequence_column)){  print("Hello") }

plots = plot_length_dist(df, args$sequence_column, args$strand_column, args$plot_column, args$abundance_column)
p = ggarrange(plotlist = plots, ncol = 2, nrow = 1)

ggsave(p, filename = paste0(args$outname), dpi = 300, height = 10, width = 10)

