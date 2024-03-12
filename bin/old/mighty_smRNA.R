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
pkgTest('gridExtra')


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
library(gridExtra)


source( paste0(this.dir(),"/mighty.R") )

parser = argparse::ArgumentParser()
parser$add_argument("-f", "--file", type = "character", required = TRUE)
parser$add_argument("-o", "--outname", type = "character", required = TRUE)
parser$add_argument("-x",type = "character", required = TRUE)
parser$add_argument("-y",type = "character", required = TRUE)


args = parser$parse_args()


axmin=2^-10
axmax=2^20

make_plots = function(data, x, y){
  
    # plot piRNAs, siRNAs, 26G RNAs, miRNA
    piRNA = xy_dge(data %>% filter(feature == "piRNA"), paste0(x), paste0(y), axmin, axmax) + ggtitle("piRNA") 
    siRNA = xy_dge(data %>% filter(feature == "siRNA"), paste0(x), paste0(y), axmin, axmax) + ggtitle("22G siRNA")
    miRNA = xy_dge(data %>% filter(feature == "miRNA"), paste0(x), paste0(y), axmin, axmax) + ggtitle("miRNA")
    transposon = xy_dge(data %>% filter(biotype == "transposable_element" | biotype == "transposon"), paste0(x), paste0(y), axmin, axmax) + ggtitle("Transposon")
  
    # plot piRNAs & siRNAs by gender
    piRNA_gender = xy_gender(data %>% filter(feature == "piRNA"), paste0(x), paste0(y), axmin, axmax) + ggtitle("piRNA")
    siRNA_gender = xy_gender(data %>% filter(feature == "siRNA"), paste0(x), paste0(y), axmin, axmax) + ggtitle("22G siRNA")

    # plot classes of small RNAs broken down into their subclasses
    pzm = plot_pzm(data, paste0(x), paste0(y), axmin, axmax)
    
    p = grid.arrange(piRNA, siRNA, miRNA, piRNA_gender, siRNA_gender, transposon, pzm, layout_matrix = rbind(c(1,2,3),c(4,5,6),c(7,7,7)))
    
  return(p)
}


data = read.delim(args$file)
   
p = make_plots(data, args$x, args$y)
    
ggsave(p, filename = args$outname, dpi = 300, height = 12, width = 14)





