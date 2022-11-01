
############# PACKAGE TEST ############# 
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE,repos='http://cran.us.r-project.org')
    if(!require(x,character.only = TRUE)) stop ("Failed to install the package. Please check the internet access or update your R if it is too old.")
  }
}

############# THEME ############# 
my_theme = function() {
    
    theme_classic() +
    theme(aspect.ratio = 1,
          axis.ticks = element_line(size = .5, color = "black"),
          axis.ticks.length=unit(-0.13, "cm"),
          text = element_text(size=10, color = "black"),
          plot.title = element_text(size = 11, hjust = 0.5, vjust = -0.5),
          axis.text = element_text(size = 11, color = "black"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0) ,color = "black"),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), color = "black"),
          panel.grid = element_blank(),
          axis.line = element_line(size = .5, colour = "black"),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 11),
          legend.key.size = unit(1,"line"),
          legend.justification = c(0, 1),
          legend.position = c(0, 1),
          legend.background = element_blank(),
          legend.key=element_blank())
}

my_theme_free_aspect = function() {
    
    theme_classic() +
    theme(
          axis.ticks = element_line(size = .5, color = "black"),
          axis.ticks.length=unit(-0.13, "cm"),
          text = element_text(size=10, color = "black"),
          plot.title = element_text(size = 10, hjust = 0.5, vjust = -0.5),
          axis.text = element_text(size = 10, color = "black"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0) ,color = "black"),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), color = "black"),
          panel.grid = element_blank(),
          axis.line = element_line(size = .5, colour = "black")
    ) 
}

############# GLOBAL PARAMETERS ############# 
BIGPOINT=1
MEDPOINT=0.6
SMALLPOINT=0.2
ALPHAPOINT=0.7
LINETYPE="solid"
LINECOLOR="grey70"
LINEALPHA=0.5
LINEWIDTH =.4
FC=1
PVALUE=0.05
RPM = 0


############# Helper Functions ###############
z_score = function(N, x, s){
  
  z = (N - x) / s
  return(z)
  
  
}

############# LENGTH DISTRIBUTION ############# 
plot_length_dist = function(dataframe, sequence_column=FALSE, strand_column=FALSE, plot_column=FALSE, abundance_column=FALSE){
  
  if (is.null(sequence_column)){
    dataframe['Sequence'] = dataframe[4] 
  } else {
    dataframe['Sequence'] = dataframe[paste0(sequence_column)]
  }
  
  if (is.null(strand_column)){
    dataframe['Strand'] = dataframe[6] 
  } else {
    dataframe['Strand'] = dataframe[paste0(strand_column)]
  }
  
  if (is.null(plot_column)){
    dataframe['smRNA_Type'] = "smRNA"
  } else {
    dataframe['smRNA_Type'] = dataframe[paste0(plot_column)]
  }
  
  if (is.null(abundance_column)){
    dataframe['Abundance'] = dataframe[5]
  } else {
    dataframe['Abundance'] = dataframe[paste0(gsub("-",".",abundance_column))]
  }
  
  dataframe$Sequence = as.character(dataframe$Sequence)

  dataframe = dataframe %>% 
    select(Sequence, smRNA_Type, Abundance) %>% 
    ungroup() %>% mutate(Length = nchar(Sequence)) %>% 
    mutate(First_NT = substr(Sequence, 1,1))

  grouped_data = dataframe %>% 
    group_by(smRNA_Type, Length, First_NT) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    ungroup() %>% 
    group_by(smRNA_Type) %>% 
    mutate(Percent = 100*Abundance/sum(Abundance)) %>% 
    ungroup()

  pie_data = grouped_data %>% 
    ungroup() %>% 
    group_by(smRNA_Type, First_NT) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    ungroup() %>% 
    mutate(Percent = 100*Abundance/sum(Abundance))
  
  cols = c("A" = "red", "G" = "green3", "C" = "violet", "T" = "turquoise")
  
  plots = list()
  p1 = ggplot(data = grouped_data, aes(x = Length, y = Abundance)) + 
    geom_bar(stat = "identity", width = 0.6, aes(fill = First_NT)) +
    my_theme() + 
    coord_cartesian() + 
    coord_capped_cart(bottom = "right", left = "top") + 
    facet_wrap(~smRNA_Type, scales = "free") + 
    scale_fill_manual(values = cols) + 
    scale_x_continuous(limits = c(15,40), breaks = c(seq(15,27,by=3),seq(30,40,5)))
  
  plots[[1]] = p1
  
  p2 = ggplot(data = pie_data, aes(x = smRNA_Type, y = Percent)) + 
    geom_bar(stat = "identity", width = 0.6, aes(fill = First_NT)) +
    my_theme() + 
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5)) + 
    coord_cartesian() + 
    coord_capped_cart(bottom = "both", left = "both") + 
    scale_fill_manual(values = cols) + 
    ylim(0,100)

  plots[[2]] = p2
  
  return(plots)
}

############# TAILING ANALYSIS ############# 
plot_tails = function(tailed, untailed, condition, samples) {

    freq = calculate_tail_frequency_total(tailed, untailed, condition, samples)
    print(freq)
    cols = c("A" = "red", "G" = "green3", "O" = "grey50", "C" = "violet", "T" = "turquoise", "TT" = "steelblue", "TTTn" = "blue")

    freq = freq %>% filter(condition %in% samples)
    freq$condition = factor(freq$condition, levels = samples)

    p = ggplot(data = freq, aes(x = condition, y = freq)) + 
      geom_bar(stat = "identity", aes(fill = tail), color = 'black') + 
      my_theme() + 
      theme(aspect.ratio = 0.5) + 
      theme(axis.text.x = element_text(angle = 60, vjust = 0.5)) + 
      coord_capped_cart(bottom = "both", left = "both") + 
      scale_fill_manual(values = cols)
    
    return(p)
 
}

calculate_tail_frequency = function(tailed, untailed, condition, samples) {
  
  # For each piRNA == sequence pair it is assigned a tail group either A, G, C, T, TT, TTT, or other
  # for each piRNA the abundance per tail group is calculated
  # for each piRNA -- tail group pair the frequency compared to the perfectly matched abundance + tailed abundance is calculated
  
  tailed$tail = as.character(tailed$tail)
  tailed = tailed %>%
    mutate(tail_group = ifelse(nchar(tail) == 1, tail,
                         ifelse(nchar(tail) >= 3 & !grepl("A|C|G", tail), "TTTn",
                                ifelse(nchar(tail) == 2 & !grepl("A|C|G", tail), "TT", "Other")))) %>%
    group_by(condition, sample, locus_id, tail_group) %>%  # get abundance per tail group
    summarise(tailed_rpm = sum(count_total_norm)) %>% 
    group_by(condition, sample, locus_id) %>% 
    mutate(total_tailed_rpm_per_piRNA = sum(tailed_rpm)) %>% 
    ungroup()

  keep = names(untailed)[(names(untailed) %in% condition$sample)]
  untailed = untailed[,c(keep, "locus_id")]
  
  untailed = untailed %>% 
    pivot_longer(!locus_id, names_to = 'sample', values_to = 'untailed_rpm') %>% 
    filter(untailed_rpm > 0) %>% 
    left_join(condition, 'sample') %>%
    group_by(sample) %>% 
    mutate(z_score_untailed_rpm = ( log2(untailed_rpm) - mean( log2(untailed_rpm) ) / sd(log2(untailed_rpm))) ) %>%
    mutate(percentile_untailed_rpm = pnorm(z_score_untailed_rpm)) %>% 
    filter(z_score_untailed_rpm >= 1)
  
  freq = tailed %>%
    inner_join(untailed, by = c('condition', 'sample', 'locus_id')) %>% 
    mutate(total_piRNA_abundance = untailed_rpm + total_tailed_rpm_per_piRNA) %>% 
    mutate(freq = 100*(tailed_rpm / total_piRNA_abundance)) %>% 
    select(sample, locus_id, tail_group, condition, tailed_rpm, total_tailed_rpm_per_piRNA, untailed_rpm, z_score_untailed_rpm, percentile_untailed_rpm, total_piRNA_abundance, freq) %>%
    filter(!is.na(freq))
  
  return(freq)

}

calculate_tail_frequency_total = function(tailed, untailed, condition, samples) {
  
  tailed$tail = as.character(tailed$tail)
  tailed = tailed %>%
    mutate(tail = ifelse(nchar(tail) == 1, tail,
                         ifelse(nchar(tail) >= 3 & !grepl("A|C|G", tail), "TTTn",
                                ifelse(nchar(tail) == 2 & !grepl("A|C|G", tail), "TT", "Other")))) %>%
    group_by(sample, condition, tail) %>%
    summarise(tailed_rpm = sum(count_total_norm))
  
  keep = names(untailed)[(names(untailed) %in% condition$sample)]
  untailed = untailed[,c(keep, "locus_id")]
  
  untailed = untailed %>% 
    pivot_longer(!locus_id, names_to = 'sample', values_to = 'untailed_rpm') %>% 
    filter(untailed_rpm > 0 ) %>% 
    left_join(condition, 'sample') %>%
    group_by(sample) %>% 
    mutate(z_score_untailed_rpm = ( log2(untailed_rpm) - mean( log2(untailed_rpm) ) / sd(log2(untailed_rpm))) ) %>%
    mutate(percentile_untailed_rpm = pnorm(z_score_untailed_rpm)) %>% 
    filter(z_score_untailed_rpm >= 1.5) %>% 
    group_by(condition, sample) %>% 
    summarise(untailed_rpm = sum(untailed_rpm)) 
  
  print(head(untailed))
  freq = tailed %>%
    left_join(untailed, by = c('condition', 'sample')) %>% 
    mutate(freq = 100*(tailed_rpm / (untailed_rpm + tailed_rpm))) %>%
    group_by(condition, tail) %>% 
    summarise(freq = mean(freq))
  
  return(freq)

}

############# HELPER FUNCTIONS ############# 
pick_pt_size = function(res){
    
    if (nrow(res)<1000){
    point_size = BIGPOINT
    } else if (nrow(res) < 3000){
    point_size = MEDPOINT
    } else {
    point_size = SMALLPOINT
    }
    
    return(point_size)
}
############# X vs Y plots ############# 
xy_plot_log2 = function(res, x, y, cols, labs, rho, m, M){

    p = ggplot(data = res %>% arrange(order), aes(x = X, y = Y, color = lab)) + 
    geom_point(aes(size = point_size), alpha = 1) + 
    scale_size_identity() + 
    scale_alpha_identity() + 
    my_theme() + 
    theme(legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(1,"line"),
          legend.justification = c(0, 1),
          legend.position = c(0, 1),
          legend.background = element_blank(),
          legend.key=element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 2, shape = 15))) + 
    coord_cartesian() + 
    #coord_capped_cart(bottom = "both", left = "both") + 
    scale_y_continuous(
      limits = c(m, M),
      trans = "log2",
      labels = trans_format("log2", math_format(2^.x)),
      breaks = trans_breaks("log2", function(x) 2^x)
    ) + 
    scale_x_continuous(
      limits = c(m, M),
      trans = "log2",
      labels = trans_format("log2", math_format(2^.x)),
      breaks = trans_breaks("log2", function(x) 2^x)
    ) +
    scale_color_manual(values = cols, labels = labs) + 
    labs(x = paste0(x), y = paste0(y), color = "") +
    #geom_abline(slope = 1, intercept = 0, lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) +
    geom_abline(slope = 1, intercept = 1, lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
    geom_abline(slope = 1, intercept = -1, lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) 
    
    if (rho != "None"){
        p = p + labs(subtitle = paste0("Pearson's rho: ", round(rho, 2)))
    } 
    
    return(p)
}

xy_plot = function(res, x, y, cols, labs, rho, point_size, m, M){

    p = ggplot(data = res %>% arrange(order), aes(x = X, y = Y, color = lab)) + 
    geom_point(size = point_size, alpha = 1) + 
    scale_size_identity() + 
    scale_alpha_identity() + 
    my_theme() + 
    theme(legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(1,"line"),
          legend.justification = c(0, 1),
          legend.position = c(0, 1),
          legend.background = element_blank(),
          legend.key=element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 2, shape = 15))) + 
    coord_cartesian() +
    scale_color_manual(values = cols, labels = labs) + 
    labs(x = paste0(x), y = paste0(y), color = "") +
    geom_abline(slope = 1, intercept = 1, lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH)
    
    if (rho != "None"){
        p = p + labs(subtitle = paste0("Pearson's rho: ", round(rho, 2)))
    } 
    
    return(p)
}

xy_dge = function(res, x, y, axmin, axmax) {
  
    point_size = pick_pt_size(res)
    
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    
    cor = cor.test(res$X, res$Y, method = "pearson")
    rho = cor(res$X, res$Y, method = "pearson")
    
    res = res %>% select(X, Y, log2FoldChange, pvalue)
    
    res = res %>% 
    mutate(lab = ifelse(log2FoldChange >= FC & pvalue < PVALUE & Y >= RPM, "Up",
                        ifelse(log2FoldChange <= -1*FC & pvalue < PVALUE & X >= RPM , "Down", "None"))) %>% 
    mutate(order = ifelse(lab == "None", 0, 1)) %>%
    mutate(size = ifelse(lab == "Up" | lab == "Down", 0.6, 0.4)) %>%
    mutate(alpha = ifelse(lab == "Up" | lab == "Down", 0.6, 0.4))
    
    res[which(res$X == 0), 'X'] = axmin
    res[which(res$Y == 0), 'Y'] = axmin
    
    up = res %>% filter(lab == "Up") 
    down = res %>% filter(lab == "Down") 
    unchanged = res %>% filter(lab == "None")
    
    print(up)
    
    cols = c("Up" = "springgreen2", 
           "Down" = "violetred1", 
           "None" = "grey80")
    
    labs = c(paste0("log2FC >= ",FC," & p < ",PVALUE," (n =",nrow(up),")"), 
           paste0("log2FC <= ",-1*FC," & p < ",PVALUE," (n =",nrow(down),")"), 
           "Other")
    
    res['point_size'] = point_size
    p = xy_plot_log2(res, x, y, cols, labs, rho, axmin, axmax)
    return(p)
}

xy_scat = function(res, x, y, axmin, axmax){

    point_size = pick_pt_size(res)
    
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    res = res %>% select(X, Y) %>% replace_na(list(X = 0, Y = 0))
    
    cor = cor.test(res$X, res$Y, method = "pearson")
    rho = cor(res$X, res$Y, method = "pearson")
    
    res['lfc'] = log2(res$Y / res$X)
    res[which(res$X == 0), 'X'] = axmin
    res[which(res$Y == 0), 'Y'] = axmin
    res['lab'] = "All"
    res['order'] = 1

    cols = c("All" = "blue")
    labs = c(paste0("All N = ",nrow(res)))
    
    res['point_size'] = point_size
    p = xy_plot_log2(res, x, y, cols, labs, rho, axmin, axmax)
    return(p)
    
}

xy_scat_basic = function(res, x, y, axmin, axmax){

    point_size = pick_pt_size(res)
    
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    
    res = res %>% select(X, Y) %>% replace_na(list(X = 0, Y = 0))
    
    cor = cor.test(res$X, res$Y, method = "pearson")
    rho = cor(res$X, res$Y, method = "pearson")
    
    res['lab'] = "All"
    res['order'] = 1

    cols = c("All" = "blue")
    labs = c(paste0("All = ",nrow(res)))
    
    res['point_size'] = point_size
    p = xy_plot(res, x, y, cols, labs, rho, axmin, axmax)
    return(p)
    
}

xy_feature = function(res, x, y, axmin, axmax) {
  
    fc = FC
    
    point_size = pick_pt_size(res)
        
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    
    res = res %>% select(X, Y, feature, class)
    
    res = res %>% 
    mutate(lab = ifelse(feature == "siRNA" & grepl("CSR", class), "siRNA_CSR", 
                        ifelse(feature == "siRNA" & grepl("PRG", class), "siRNA_PRG1",
                            ifelse(feature == "piRNA", "piRNA", 
                                ifelse(feature == "miRNA", "miRNA", 
                                    ifelse(feature == "risiRNA", "risiRNA", 
                                        ifelse(feature == "RepBase", "RepBase", "Other"))))))) %>% 
    mutate(order = ifelse(lab == "siRNA_CSR",3, 
                          ifelse(lab == "siRNA_PRG1",2,
                                ifelse(lab == "piRNA",1,
                                    ifelse(lab == "miRNA",4,
                                        ifelse(feature == "risiRNA", 5, 
                                            ifelse(lab == "RepBase",4, 0))))))) %>%
    mutate(point_size = ifelse(lab == "RepBase"  | lab == "miRNA" | lab == "risiRNA", 1, 0.2)) %>%
    arrange(order)
    
    res[which(res$X == 0), 'X'] = axmin
    res[which(res$Y == 0), 'Y'] = axmin
    
    cols = c("siRNA_CSR" = "royalblue1", 
           "siRNA_PRG1" = "orange", 
           "piRNA" = "green", 
           "miRNA" = "red",
           "risiRNA" = "magenta",
           "RepBase" = "orange", 
           "Other" = "grey80")
    
    n_csr =  res %>% filter(lab == "siRNA_CSR") %>% nrow()
    n_prg = res %>% filter(lab == "siRNA_PRG1") %>% nrow()
    n_piRNA = res %>%  filter(lab == "piRNA") %>% nrow()
    n_miRNA = res %>%  filter(lab == "miRNA") %>% nrow()
    n_repbase = res %>% filter(lab == "RepBase") %>% nrow()
    n_risi = res %>% filter(lab == "risiRNA") %>% nrow()
    n_other = res %>% filter(lab == "Other") %>% nrow()
    
    labs = c(paste0("CSR-1 22G"," (",n_prg,")"),
           paste0("PRG-1 22G"," (", n_csr,")"), 
           paste0("piRNA", " (", n_piRNA,")"), 
           paste0("miRNA", " (", n_miRNA,")"), 
           paste0("risiRNA", " (", n_risi,")"), 
           paste0("RepBase", " (", n_repbase,")"), 
           paste0("Other"," (",n_other,")"))
    
    p = xy_plot_log2(res, x, y, cols, labs, "None", axmin, axmax)
    return(p)
}

xy_highlight = function(res, x, y, feat, axmin, axmax) {
  
    fc = FC
    
    point_size = pick_pt_size(res)
        
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    
    res = res %>% select(X, Y, feature, class)
    
    res = res %>% 
    mutate(lab = ifelse(feature == feat, "Feature", "Other")) %>%
    mutate(order = ifelse(lab == "Feature",1,0)) %>%
    mutate(point_size = ifelse(lab == "Feature", point_size, 0.2)) %>%
    arrange(order)
    
    res[which(res$X == 0), 'X'] = axmin
    res[which(res$Y == 0), 'Y'] = axmin
    
    cols = c("Feature" = "blue", 
           "Other" = "grey80")
           
    n_feature =  res %>% filter(lab == "Feature") %>% nrow()
    n_other = res %>% filter(lab == "Other") %>% nrow()
    
    labs = c(paste0(feat," (",n_feature,")"),
           paste0("Other"," (", n_other,")"))
    
    p = xy_plot_log2(res, x, y, cols, labs, "None", axmin, axmax)
    return(p)
}

xy_siRNA_class = function(res, x, y, axmin, axmax) {
  
    fc = FC
    
    point_size = pick_pt_size(res)
    
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    
    res = res %>% select(X, Y, class)
    
    res = res %>% 
    mutate(lab = ifelse(grepl("WAGO1", class), "WAGO1", 
                        ifelse(grepl("PRG", class), "PRG1", 
                               ifelse(grepl("CSR", class), "CSR", "Other")))) %>% 
    mutate(order = ifelse(lab == "WAGO", 1, 
                          ifelse(lab == "PRG1", 3,
                                 ifelse(lab == "CSR", 2, 0)))) 
    
    res[which(res$X == 0), 'X'] = axmin
    res[which(res$Y == 0), 'Y'] = axmin
    
    cols = c("WAGO1" = "green", 
           "PRG1" = "red", 
           "CSR" = "blue", 
           "Other" = "grey80")
    
    n_wago =  res %>% filter(grepl("WAGO1", class)) %>% nrow()
    n_prg = res %>% filter(grepl("PRG", class)) %>% nrow()
    n_csr = res %>% filter(grepl("CSR", class)) %>% nrow()
    n_other = res %>% filter(!grepl("CSR", class)) %>% filter(!grepl("PRG", class)) %>% filter(!grepl("WAGO", class)) %>% nrow()
    
    labs = c(paste0("WAGO1"," (",n_wago,")"),
           paste0("PRG"," (", n_prg,")"), 
           paste0("CSR", " (", n_csr,")"), 
           paste0("Other"," (",n_other,")"))
           
    res['point_size'] = point_size
    p = xy_plot_log2(res, x, y, cols, labs, "None", axmin, axmax)
    return(p)
}

xy_26G_class = function(res, x, y, axmin, axmax) {
  
    fc = FC
    
    point_size = pick_pt_size(res)
    
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    
    res = res %>% select(X, Y, log2FoldChange, pvalue, class)
    
    res = res %>% 
    mutate(lab = ifelse(grepl("ALG", class), "ALG", 
                        ifelse(grepl("ERGO", class), "ERGO", "Other"))) %>% 
    mutate(order = ifelse(lab == "Other", 1, 
                          ifelse(lab == "ALG", 2,
                                 ifelse(lab == "ERGO", 3, 0))))
    
    res[which(res$X == 0), 'X'] = axmin
    res[which(res$Y == 0), 'Y'] = axmin
    
    cols = c("ALG" = "green3", 
           "ERGO" = "blue", 
           "Other" = "grey80")
    
    n_alg =  res %>% filter(grepl("ALG", class)) %>% nrow()
    n_ergo = res %>% filter(grepl("ERGO", class)) %>% nrow()
    n_other = res %>% filter(!grepl("ERGO", class)) %>% filter(!grepl("ALG", class)) %>% nrow()
    
    labs = c(paste0("ALG"," (",n_alg,")"),
           paste0("ERGO"," (", n_ergo,")"), 
           paste0("Other", " (", n_other,")"))
           
    res['point_size'] = point_size
    p = xy_plot_log2(res, x, y, cols, labs, "None", axmin, axmax)
    return(p)
    }
    
xy_gender = function(res, x, y, axmin, axmax) {
    
    fc = FC
    
    point_size = pick_pt_size(res)
    
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    
    res = res %>% select(X, Y, class)
    
    res = res %>% 
    mutate(lab = ifelse(grepl("Spermatogenic", class), "Sperm", 
                        ifelse(grepl("Oogenic", class), "Oocyte", 
                               ifelse(grepl("GenderNeutral", class), "Neutral", "Other")))) %>% 
    mutate(order = ifelse(lab == "Sperm", 4, 
                          ifelse(lab == "Oocyte", 3, 
                                 ifelse(lab == "Neutral", 2, 1))))
    
    res[which(res$X == 0), 'X'] = axmin
    res[which(res$Y == 0), 'Y'] = axmin
    
    cols = c("Sperm" = "blue", 
           "Oocyte" = "magenta", 
           "Neutral" = "green3",
           "Other" = "grey70")
    
    n_sperm =  res %>% filter(grepl("Spermatogenic", class)) %>% nrow()
    n_oocyte = res %>% filter(grepl("Oogenic", class)) %>% nrow()
    n_neutral = res %>% filter(grepl("GenderNeutral", class)) %>% nrow()
    n_other = res %>% filter(!grepl("Spermatogenic", class)) %>% filter(!grepl("Oogenic", class)) %>% filter(!grepl("neutral", class)) %>% nrow()
    
    labs = c(paste0("Sperm"," (",n_sperm,")"),
           paste0("Oocyte"," (", n_oocyte,")"), 
           paste0("Gender_Neutral"," (", n_neutral,")"), 
           paste0("Other", " (", n_other,")"))
          
    res['point_size'] = point_size
    p = xy_plot_log2(res, x, y, cols, labs, "None", axmin, axmax)
    return(p)
}

############## MA PLOTS #############
MA = function(res, x, y, cols, labs, point_size, m, M) {
    
    p = ggplot(data = res %>% arrange(order), aes(x = basemean, y = log2FoldChange, color = lab)) + 
    geom_point(size = point_size, alpha = 1) + 
    scale_size_identity() + 
    scale_alpha_identity() + 
    my_theme() + 
    theme(legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(1,"line"),
          legend.justification = c(0, 1),
          legend.position = c(0, 1),
          legend.background = element_blank(),
          legend.key=element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 2, shape = 15))) + 
    coord_cartesian() + 
    coord_capped_cart(bottom = "both", left = "both") + 
    scale_x_continuous(
      limits = c(m, M),
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      breaks = trans_breaks("log10", function(x) 10^x)
    ) +
    scale_color_manual(values = cols, labels = labs) + 
    labs(x = "Basemean RPM", y = paste0("log2FoldChange( ",x," / ", y ," )"), color = "") +
    geom_hline(yintercept = 1, lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
    geom_hline(yintercept = -1, lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH)
    
    return(p)
    
}

MA_plot_siRNA_class = function(res, x, y, axmin, axmax) {
     
     point_size = pick_pt_size(res)

    res = res %>% select(log2FoldChange, pvalue, basemean, class)
    
    res = res %>% 
    mutate(lab = ifelse(grepl("WAGO", class), "WAGO", 
                        ifelse(grepl("PRG", class), "PRG1", 
                               ifelse(grepl("CSR", class), "CSR", "Other")))) %>% 
    mutate(order = ifelse(lab == "WAGO", 1, 
                          ifelse(lab == "PRG1", 3,
                                 ifelse(lab == "CSR", 2, 0))))
                                 
    cols = c("WAGO" = "green3", 
       "PRG1" = "blue", 
       "CSR" = "magenta", 
       "Other" = "grey80")
    
    n_wago =  res %>% filter(grepl("WAGO", class)) %>% nrow()
    n_prg = res %>% filter(grepl("PRG", class)) %>% nrow()
    n_csr = res %>% filter(grepl("CSR", class)) %>% nrow()
    n_other = res %>% filter(!grepl("CSR", class)) %>% filter(!grepl("PRG", class)) %>% filter(!grepl("WAGO", class)) %>% nrow()
    
    labs = c(paste0("WAGO"," (",n_wago,")"),
       paste0("PRG"," (", n_prg,")"), 
       paste0("CSR", " (", n_csr,")"), 
       paste0("Other"," (",n_other,")"))
       
    p = MA(res, x, y, cols, labs, point_size, axmin, axmax)
    return(p)
}

MA_dge =  function(res, x, y, axmin, axmax) {
  
    point_size = pick_pt_size(res)

    res = res %>% select(log2FoldChange, pvalue, basemean)
    
    res = res %>% 
    mutate(lab = ifelse(log2FoldChange >= FC & pvalue < PVALUE, "Up",
                        ifelse(log2FoldChange <= -1*FC & pvalue < PVALUE, "Down", "None"))) %>% 
    mutate(order = ifelse(lab == "None", 0, 1)) %>%
    mutate(size = ifelse(lab == "Up" | lab == "Down", 0.6, 0.4)) %>%
    mutate(alpha = ifelse(lab == "Up" | lab == "Down", 0.6, 0.4))
    
    up = res %>% filter(lab == "Up") 
    down = res %>% filter(lab == "Down") 
    unchanged = res %>% filter(lab == "None")
    
    up_percent = round(100*( nrow(up) / nrow(res) ))
    down_percent = round(100*( nrow(down) / nrow(res) ))
    
    unchanged_percent = 100-(up_percent + down_percent)
    
    cols = c("Up" = "seagreen3", 
           "Down" = "violetred1", 
           "None" = "grey80")
    
    labs = c(paste0("log2FC >= ",FC," & p < ",PVALUE," (n =",nrow(up),")"), 
           paste0("log2FC <= ",-1*FC," & p < ",PVALUE," (n =",nrow(down),")"), 
           "Other")
    
    p = MA(res, x, y, cols, labs, point_size, axmin, axmax)
    return(p)
}


############## BOXPLOTS #############
plot_smRNA_boxplot = function(res, X, Y) {
  
  wago = res %>% filter(grepl("WAGO", class) & feature == "siRNA") %>% mutate(group = "WAGO") 
  prg1 = res %>% filter(grepl("PRG", class) & feature == "siRNA") %>% mutate(group = "PRG")
  csr = res %>% filter(grepl("CSR", class) & feature == "siRNA") %>% mutate(group = "CSR")
  alg_22G = res %>% filter(grepl("ALG", class) & feature == "siRNA") %>% mutate(group = "ALG_22G")
  alg_26G = res %>% filter(grepl("ALG", class) & feature == "siRNA26G") %>% mutate(group = "ALG_26G")
  ergo_22G = res %>% filter(grepl("ERGO", class) & feature == "siRNA") %>% mutate(group = "ERGO_22G")
  ergo_26G = res %>% filter(grepl("ERGO", class) & feature == "siRNA26G") %>% mutate(group = "ERGO_26G")
  
  dat = rbind(wago, prg1, csr, alg_22G, alg_26G, ergo_22G, ergo_26G)
  
  p = ggplot(data = dat, aes(x = group, y = log2FoldChange)) + 
    stat_boxplot(geom = "errorbar", width = 0.3) + 
    geom_boxplot(width = 0.4) + 
    stat_summary(fun.y=mean,col='red',geom='point') + 
    my_theme() + 
    theme(aspect.ratio = 1, axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)) + 
    coord_capped_cart(ylim = c(-10, 8), bottom = "both", left = "both") + 
    scale_y_continuous(breaks = seq(-10, 8, by = 2)) + 
    geom_hline(yintercept = 0, color = "skyblue", lwd = 0.7, linetype = "dashed")
  
  return(p)
}

plot_smRNA_boxplot_RPM = function(res, x, y, axmin, axmax) {
  
  res['X'] = res[ , paste0(x)]
  res['Y'] = res[ , paste0(y)]
  
  res = res %>% select(X, Y, class, feature)
  
  wago = res %>% filter(grepl("WAGO", class) & feature == "siRNA") %>% mutate(group = "WAGO") %>% pivot_longer(!c(class, feature, group), names_to = "sample", values_to = "count")
  prg1 = res %>% filter(grepl("PRG", class) & feature == "siRNA") %>% mutate(group = "PRG") %>% pivot_longer(!c(class, feature, group), names_to = "sample", values_to = "count")
  csr = res %>% filter(grepl("CSR", class) & feature == "siRNA") %>% mutate(group = "CSR")  %>% pivot_longer(!c(class, feature, group), names_to = "sample", values_to = "count")
  alg_22G = res %>% filter(grepl("ALG", class) & feature == "siRNA") %>% mutate(group = "ALG_22G") %>% pivot_longer(!c(class, feature, group), names_to = "sample", values_to = "count")
  alg_26G = res %>% filter(grepl("ALG", class) & feature == "siRNA26G") %>% mutate(group = "ALG_26G") %>% pivot_longer(!c(class, feature, group), names_to = "sample", values_to = "count")
  ergo_22G = res %>% filter(grepl("ERGO", class) & feature == "siRNA") %>% mutate(group = "ERGO_22G") %>% pivot_longer(!c(class, feature, group), names_to = "sample", values_to = "count")
  ergo_26G = res %>% filter(grepl("ERGO", class) & feature == "siRNA26G") %>% mutate(group = "ERGO_26G")  %>% pivot_longer(!c(class, feature, group), names_to = "sample", values_to = "count")
  
  dat = rbind(wago, prg1, csr, alg_22G, alg_26G, ergo_22G, ergo_26G)
  
  dat = dat %>% mutate(sample = ifelse(sample == "X", paste0("A_",x), paste0("B_",y)))
  
  stat.test <- dat %>%
    group_by(group) %>%
    wilcox_test(.,count ~ sample, alternative = "two.sided") %>%
    adjust_pvalue(method = "none") %>%
    add_significance() %>% 
    add_xy_position(x = "group", fun = "max") %>% 
    mutate(y.position = log2(y.position))
  
  dodge = 0.8
  p = ggplot(data = dat %>% arrange(sample), aes(x = group, y = count)) + 
    stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(dodge), aes(fill = sample)) + 
    geom_boxplot(width = 0.6, aes(fill = sample), outlier.size = 0.1, outlier.alpha = 0.5, outlier.color = "grey50", position = position_dodge(dodge)) + 
    stat_summary(fun.y=mean,col='magenta',geom='point', aes(fill = sample), position = position_dodge(dodge)) + 
    my_theme() + 
    theme(aspect.ratio = 1, axis.text.x = element_text(angle = 60, vjust = 0.5, size = 8)) + 
    coord_capped_cart(bottom = "both", left = "both") + 
    scale_y_continuous(
      limits = c(axmin, axmax),
      trans = "log2",
      labels = trans_format("log2", math_format(2^.x)),
      breaks = trans_breaks("log2", function(x) 2^x)) + 
    scale_fill_brewer(palette = "Set3") + 
    stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) + 
    labs(x = "", y = "RPM")
  # stat_compare_means(comparisons = list(c(paste0("A_",x), paste0("B_",y)))) can also do this
  
  return(p)
}


plot_gender_boxplot = function(res, x, y, axmin, axmax) {
  
  res['X'] = res[ , paste0(x)]
  res['Y'] = res[ , paste0(y)]
  
  res = res %>% select(X, Y, class, feature)
  
  sperm = res %>% filter(grepl("Spermatogenic", class)) %>% mutate(group = "sperm") %>% pivot_longer(!c(class, feature, group), names_to = "sample", values_to = "count")
  oocyte = res %>% filter(grepl("Oogenic", class)) %>% mutate(group = "oocyte") %>% pivot_longer(!c(class, feature, group), names_to = "sample", values_to = "count")
  neutral = res %>% filter(grepl("neutral", class)) %>% mutate(group = "neutral")  %>% pivot_longer(!c(class, feature, group), names_to = "sample", values_to = "count")
  other = res %>% filter(!grepl("neutral", class)) %>% filter(!grepl("Oogenic", class)) %>% filter(!grepl("Spermatogenic", class)) %>% mutate(group = "Other")  %>% pivot_longer(!c(class, feature, group), names_to = "sample", values_to = "count")
  
  dat = rbind(sperm, oocyte, neutral, other)
  
  dat = dat %>% mutate(sample = ifelse(sample == "X", paste0("A_",x), paste0("B_",y)))
  
  my_comparisons = list(c(paste0("A_",x), paste0("B_",y)))
  
  stat.test <- dat %>%
    group_by(group) %>%
    wilcox_test(.,count ~ sample, alternative = "two.sided") %>%
    adjust_pvalue(method = "none") %>%
    add_significance() %>% 
    add_xy_position(x = "group", fun = "max") %>% 
    mutate(y.position = log2(y.position))
  
  dodge = 0.8
  p = ggplot(data = dat %>% arrange(sample), aes(x = group, y = count)) + 
    stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(dodge), aes(fill = sample)) + 
    geom_boxplot(width = 0.6, aes(fill = sample), outlier.size = 0.1, outlier.alpha = 0.5, outlier.color = "grey50", position = position_dodge(dodge)) + 
    stat_summary(fun.y=mean,col='magenta',geom='point', aes(fill = sample), position = position_dodge(dodge)) + 
    my_theme() + 
    theme(aspect.ratio = 1, axis.text.x = element_text(angle = 60, vjust = 0.5, size = 8)) + 
    coord_capped_cart(bottom = "both", left = "both") + 
    scale_y_continuous(
      limits = c(axmin, axmax),
      trans = "log2",
      labels = trans_format("log2", math_format(2^.x)),
      breaks = trans_breaks("log2", function(x) 2^x)) + 
    scale_fill_brewer(palette = "Set3") + 
    #facet_wrap(~group) + 
    stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) + 
    labs(x = "", y = "RPM")
  
  # stat_compare_means(comparisons = list(c(paste0("A_",x), paste0("B_",y)))) can also do this
  
  return(p)
}

plot_pzm = function(res, x, y, axmin, axmax) {
  
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    
    res = res %>% select(X, Y, class, feature, locus_id)
    
    # P-granule = CSR-1, WAGO-1, PRG-1, ALG-3/4
    # M-granule = MUT-16 
    # Z-granule = WAGO-4
    # Others / Nuclear = ERGO-1, HRDE-1
    
    # First do PG 
    CSR1 = res %>% filter(grepl("CSR", class) & feature == "siRNA") %>% mutate(group = "CSR1") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    PRG1 = res %>% filter(grepl("PRG", class) & feature == "siRNA") %>% mutate(group = "PRG1") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    WAGO1 = res %>% filter(grepl("WAGO1", class) & feature == "siRNA") %>% mutate(group = "WAGO1") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    ALG34 = res %>% filter(grepl("ALG34", class) & feature == "siRNA") %>% mutate(group = "ALG34") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    ALG34_26G = res %>% filter(grepl("ALG34", class) & feature == "siRNA26G") %>% mutate(group = "ALG34_26G") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    ERGO = res %>% filter(grepl("ERGO", class) & feature == "siRNA") %>% mutate(group = "ERGO") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    ERGO_26G = res %>% filter(grepl("ERGO", class) & feature == "siRNA26G") %>% mutate(group = "ERGO_26G") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    piRNA = res %>% filter(feature == "piRNA") %>% mutate(group = "21URNA") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    P = rbind(CSR1, PRG1, WAGO1, ALG34, ALG34_26G, piRNA, ERGO, ERGO_26G)
    P$compartment = "P_granule"
    
    # M 
    MUT16 = res %>% filter(grepl("MUT16", class) & feature == "siRNA") %>% mutate(group = "MUT16") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    M = MUT16 
    M$compartment = "M_focus"
    
    # Z 
    WAGO4 = res %>% filter(grepl("WAGO4", class) & feature == "siRNA") %>% mutate(group = "WAGO4") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    Z = WAGO4
    Z$compartment = "Z_granule"

    # Others
    HRDE1 = res %>% filter(grepl("HRDE1|WAGO9", class) & feature == "siRNA") %>% mutate(group = "HRDE1") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    HISTONES = res %>% filter(grepl("his|htz", locus_id) & feature == "siRNA") %>% mutate(group = "Histone_22G") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    HISTONES_26G = res %>% filter(grepl("his|htz", locus_id) & feature == "siRNA26G") %>% mutate(group = "Histone_26G") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    RISIRNA = res %>% filter(feature == "risiRNA") %>% mutate(group = "risi_RNA") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    RISIRNA_26G = res %>% filter(feature == "risiRNA26G") %>% mutate(group = "risi_26GRNA") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    Others = rbind(HRDE1,HISTONES, HISTONES_26G, RISIRNA,RISIRNA_26G)
    Others$compartment = "Other/Nuclear"

    dat = rbind(P, M, Z, Others)
    
    dat = dat %>% mutate(sample = ifelse(sample == "X", paste0("A_",x), paste0("B_",y)))
    
    my_comparisons = list(c(paste0("A_",x), paste0("B_",y)))
    
    compartments = levels(factor(dat$compartment))
    k = 0
    res.stat.test = ""
    for (i in compartments){
        print(i)
        df = dat %>% filter(compartment == i)
        stat.test <- df %>%
            group_by(group) %>% 
            wilcox_test(count ~ sample, alternative = "two.sided") %>%
            adjust_pvalue(method = "none") %>%
            add_significance() %>%
            add_x_position(x = "group") %>%
            add_y_position(y.trans = function(x){log2(x)}, fun = 'max') %>%
            mutate(compartment = i)
        
        if (k == 0){
            res.stat.test = stat.test
        } else {
            res.stat.test = rbind(res.stat.test, stat.test)
        }
        k = k + 1
    }
    
    cols = c("turquoise", "red")
    labs = c(x, y)
    
    dodge = 0.8
    p = ggplot(data = dat, aes(x = group, y = count)) + 
        stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(dodge), aes(fill = sample)) + 
        geom_boxplot(width = 0.6, aes(fill = sample), outlier.size = 0.1, outlier.alpha = 0.5, outlier.color = "grey50", position = position_dodge(dodge)) + 
        stat_summary(fun=mean, color = "black", geom='point', size = 2, aes(fill = sample), position = position_dodge(dodge)) + 
        my_theme_free_aspect() +
        theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 8)) + 
        coord_capped_cart(left = "both") + 
        scale_y_continuous(
            limits = c(axmin,axmax),
            trans = "log2",
            labels = trans_format("log2", math_format(2^.x)),
            breaks = trans_breaks("log2", function(x) 2^x)) + 
        facet_grid(~compartment, scales = "free", space='free') + 
        stat_pvalue_manual(res.stat.test,label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.05, tip.length = 0.01, step.group.by = 'group') + 
        scale_fill_manual(values = cols, labels = labs)
    
    return(p)
}


plot_pzm_multi = function(res, axmin, axmax) {
  
    res = res %>% select(-seq_id, -gene_name, -biotype)
    # P-granule = CSR-1, WAGO-1, PRG-1, ALG-3/4
    # M-granule = MUT-16 
    # Z-granule = WAGO-4
    # Others / Nuclear = ERGO-1, HRDE-1
    
    # First do PG 
    CSR1 = res %>% filter(grepl("CSR", class) & feature == "siRNA") %>% mutate(group = "CSR1") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    PRG1 = res %>% filter(grepl("PRG", class) & feature == "siRNA") %>% mutate(group = "PRG1") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    WAGO1 = res %>% filter(grepl("WAGO1", class) & feature == "siRNA") %>% mutate(group = "WAGO1") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    ERGO = res %>% filter(grepl("ERGO", class) & feature == "siRNA") %>% mutate(group = "ERGO") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    ERGO_26G = res %>% filter(grepl("ERGO", class) & feature == "siRNA26G") %>% mutate(group = "ERGO_26G") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    piRNA = res %>% filter(feature == "piRNA") %>% mutate(group = "21URNA") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    ALG34 = res %>% filter(grepl("ALG34", class) & feature == "siRNA") %>% mutate(group = "ALG34") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    ALG34_26G = res %>% filter(grepl("ALG34", class) & feature == "siRNA26G") %>% mutate(group = "ALG34_26G") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    P = rbind(CSR1, PRG1, WAGO1, ALG34, ALG34_26G, ERGO, ERGO_26G, piRNA)
    P$compartment = "P_granule"
    
    # M 
    MUT16 = res %>% filter(grepl("MUT16", class) & feature == "siRNA") %>% mutate(group = "MUT16") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    M = MUT16 
    M$compartment = "M_focus"
    
    # Z 
    WAGO4 = res %>% filter(grepl("WAGO4", class) & feature == "siRNA") %>% mutate(group = "WAGO4") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    Z = WAGO4
    Z$compartment = "Z_granule"

    # Others
    HRDE1 = res %>% filter(grepl("HRDE1|WAGO9", class) & feature == "siRNA") %>% mutate(group = "HRDE1") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    HISTONES = res %>% filter(grepl("his|htz", locus_id) & feature == "siRNA") %>% mutate(group = "Histone_22G") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    HISTONES_26G = res %>% filter(grepl("his|htz", locus_id) & feature == "siRNA26G") %>% mutate(group = "Histone_26G") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    RISIRNA = res %>% filter(feature == "risiRNA") %>% mutate(group = "risi_RNA") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    RISIRNA_26G = res %>% filter(feature == "risiRNA26G") %>% mutate(group = "risi_26GRNA") %>% pivot_longer(!c(class, feature, group, locus_id), names_to = "sample", values_to = "count")
    Others = rbind(HRDE1,HISTONES, HISTONES_26G, RISIRNA, RISIRNA_26G)
    Others$compartment = "Other/Nuclear"

    dat = rbind(P, M, Z, Others)
    
    compartments = levels(factor(dat$compartment))
    k = 0
    res.stat.test = ""
    for (i in compartments){
        print(i)
        df = dat %>% filter(compartment == i)
        stat.test <- df %>%
            group_by(group) %>% 
            wilcox_test(count ~ sample, alternative = "two.sided") %>%
            adjust_pvalue(method = "none") %>%
            add_significance() %>%
            add_x_position(x = "group") %>%
            add_y_position(y.trans = function(x){log2(x)}, fun = 'max') %>%
            mutate(compartment = i)
        
        if (k == 0){
            res.stat.test = stat.test
        } else {
            res.stat.test = rbind(res.stat.test, stat.test)
        }
        k = k + 1
    }
    
    dodge = 0.8
    p = ggplot(data = dat, aes(x = group, y = count)) + 
        stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(dodge), aes(fill = sample)) + 
        geom_boxplot(width = 0.6, aes(fill = sample), outlier.size = 0.1, outlier.alpha = 0.5, outlier.color = "grey50", position = position_dodge(dodge)) + 
        stat_summary(fun=mean, aes(fill = sample), color = 'black', geom='errorbar', linetype = 2, position = position_dodge(dodge)) + 
        my_theme_free_aspect() +
        theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 8)) + 
        coord_capped_cart(left = "both") + 
        scale_y_continuous(
            trans = "log2",
            labels = trans_format("log2", math_format(2^.x)),
            breaks = trans_breaks("log2", function(x) 2^x)) + 
        facet_grid(~compartment, scales = "free", space='free') + 
        stat_pvalue_manual(res.stat.test,label = "p.adj.signif", hide.ns = TRUE, step.increase = 0.05, tip.length = 0.01, step.group.by = 'group') + 
        scale_fill_brewer(palette = "Set1") 

    
    return(p)
}

plot_boxplot_ylog2 = function(data, x, y) {
    
    data['sample'] = data[x]
    data['count'] = data[y]
    
    data = data %>% filter(count > 0)
    
    dodge = 0.8
    p = ggplot(data, aes(x = fct_reorder(sample, count, .fun = median, .desc = TRUE), y = count)) + 
        stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(dodge)) + 
        geom_boxplot(width = 0.4, outlier.size = 0.1, outlier.alpha = 0.5, outlier.color = "grey50", position = position_dodge(dodge)) + 
        stat_summary(fun=mean, color = 'red', geom='point', position = position_dodge(dodge)) + 
        my_theme() + 
        theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 8)) + 
        scale_y_continuous(
            trans = "log2",
            labels = trans_format("log2", math_format(2^.x)),
            breaks = trans_breaks("log2", function(x) 2^x)) 

    return(p)
}


############## MISC. plots #############
get_cartesian = function(a){
  
  return( t(combn(unique(a), 2)) )
  
}

plot_biorep_cor = function(df){
  
  plotlist = list()
  i = 1
  comps = get_cartesian(colnames(df))
  for (i in 1:nrow(comps)) {
    
    x = comps[i,][[1]]
    y = comps[i,][[2]]
    dsub = df %>% select(x, y) 
    
    c = round(cor(dsub[paste0(x)], dsub[paste0(y)]), 3)
    
    p = plot_scat(dsub, x, y, 2^-5, 2^15) 
    plotlist[[i]] = p
    i = i + 1
  }
  if (round(length(plotlist)/2) <= 4){
    z = ggarrange(plotlist=plotlist)
  } else {
    z = ggarrange(plotlist=plotlist, ncol = round(length(plotlist)/2), nrow = round(length(plotlist)/2))
  }
  
  return(z)
}









