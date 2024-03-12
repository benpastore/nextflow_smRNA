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
parser$add_argument("-i", "--input", type = "character", required = TRUE)
parser$add_argument("-o", "--outname", type = "character", required = TRUE)
parser$add_argument("-x",type = "character", required = TRUE)
parser$add_argument("-f",type = "character", required = FALSE)
parser$add_argument("-y",type = "character", required = TRUE)


args = parser$parse_args()

axmin=2^-5
axmax=2^15

data = read.delim(args$input)
if (!is.null(args$f)){
    p = plot_deseq_res(data %>% filter(feature == args$f) , args$x, args$y, axmin, axmax)
} else {
    p = plot_deseq_res(data, args$x, args$y, axmin, axmax)
}

ggsave(p, filename = args$outname, dpi = 300, height = 6, width = 6)




