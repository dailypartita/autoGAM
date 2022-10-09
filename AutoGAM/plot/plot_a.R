library(mgcv)
library(plyr)
library(dplyr)
library(knitr)
library(purrr)
library(scales)
library(readxl)
library(ggplot2)
library(stringi)
library(reshape2)
library(kableExtra)
P = rprojroot::find_rstudio_root_file
R.utils::sourceDirectory(P("R"))
gam_ALL_total = readRDS(P("temp/GAM0728/GAM_ALL.rds"))$model
gam_DNA_total = readRDS(P("temp/GAM0728/GAM_DNA.rds"))$model
gam_RNA_total = readRDS(P("temp/GAM0728/GAM_RNA.rds"))$model
gam_ALL = gam_ALL_total[[1]]
gam_DNA = gam_DNA_total[[1]]
gam_RNA = gam_RNA_total[[1]]
summary(gam_ALL)
summary(gam_DNA)
summary(gam_RNA)
CV = function(bgam){
  cv_gam(bgam) %>%
    mutate(p_value = round(p_value, digits=3)) %>%
    kable(align=c("l","r", "r", "r", "r"), format = "html") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
}
# CV(gam_ALL)
# CV(gam_DNA)
# CV(gam_RNA)
# gam_COV = readRDS(P("temp/GAM0408/GAM_virus_Coronaviridae_2022-04-08.rds"))$model[[1]]
# gam_COOSP = readRDS(P("temp/GAM0408/GAM_virus_COOSP_2022-04-06.rds"))$model[[1]]
# gam_ZOO = readRDS(P("temp/GAM0408/GAM_virus_ZOO_2022-04-08.rds"))$model[[1]]
# gam_HIP = readRDS(P("temp/GAM0408/GAM_virus_HIP_2022-04-15.rds"))$model[[1]]
gam_list = list(`Total virus (30.0%)` = gam_ALL,
                `DNA virus (49.0%)` = gam_DNA,
                `RNA virus (21.2%)` = gam_RNA)
rel_deviance = function(bgam){
  rd = get_relative_contribs(bgam) %>%
    mutate(rel_deviance_explained = signif(rel_deviance_explained*100, digits=3)) %>% 
    rename(Term = term, `Relative deviance explained` = rel_deviance_explained)
  return(rd)
}
res = lapply(gam_list, rel_deviance)
res2 = ldply(res, data.frame)
names(res2) = c("virus_type", "Factors", "Relative.deviance.explained")

if (T) {
  # res3 = res2 %>%
  #   ungroup() %>%
  #   left_join(new_name, by=c("Factors" = "variable")) %>%
  #   mutate(Factors = ifelse(is.na(best_name), Factors, best_name),
  #          ) %>%
  #   filter(!Factors %in% c("RNA",
  #                          "X15.1_LitterSize",
  #                          "Genus_Coleura",
  #                          "Genus_Hipposideros",
  #                          "Genus_Miniopterus",
  #                          "Clean Base")) %>%
  #   dplyr::select(virus_type, Factors, Relative.deviance.explained)
  # res4 = res3[order(-res3[,3]),]

  theme_mine <- function(base_size = 18, base_family = "Helvetica") {
    # Starts with theme_grey and then modify some parts
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(
        strip.background = element_blank(),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14,hjust=1),
        axis.ticks =  element_line(colour = "black"), 
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16,angle=90),
        panel.background = element_blank(), 
        panel.border =element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.margin = unit(1.0, "lines"), 
        plot.background = element_blank(), 
        plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)
      )
  }
  
  var_data = read.csv(P('plot/var_data_0811.csv'))
  
  var_order = var_data %>%
    group_by(term) %>%
    mutate(avg_exp=as.numeric(summary(rel_deviance_explained)[4])) %>%
    ungroup() %>%
    select(term, avg_exp) %>%
    distinct()
  var_order = var_order[order(var_order$avg_exp),]
  var_data$term = factor(var_data$term, levels=var_order$term)
  var_data$model_type = factor(var_data$model_type, levels=c('Total virus', 'DNA virus', 'RNA virus'))
  
  x11()
  ggplot(var_data, aes(x = rel_deviance_explained, y = term)) +
    geom_boxplot(size=1, aes(color = model_type)) +
    xlab("Relative deviance explained (%)") +
    ylab("Variable") +
    facet_wrap( ~ model_type) +
    theme_mine()
  ggsave("plot/gam_summary_0725_V1.pdf")
}






