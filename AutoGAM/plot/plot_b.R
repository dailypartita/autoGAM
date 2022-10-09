library(mgcv)
library(mgcViz)
library(knitr)
library(dplyr)
library(purrr)
library(purrrlyr)
library(caret)
library(MuMIn)
library(ggplot2)
library(ggpubr)
library(stringi)
library(parallel)
library(tidyverse)
library(DataCombine)
library(tidyr)
library(purrrlyr)
library(cowplot)
library(viridis)
library(svglite)
library(magrittr)
P = rprojroot::find_rstudio_root_file
R.utils::sourceDirectory(P("R"))
bgam = readRDS(P("temp/GAM0728/GAM_ALL.rds"))$model[[1]]
get_relative_contribs(bgam)
cross_validation(bgam)
##########################画图#######################
SHOW_DEV_EXPL = FALSE
partials_theme = theme(text = element_text(family="Helvetica", size=11),
                       panel.border=element_blank(),
                       panel.background=element_blank(),
                       axis.title.y = element_blank(),
                       panel.grid = element_blank(),
                       axis.ticks.x = element_line(size=0.3),
                       axis.ticks.y = element_blank(),
                       axis.text.y = element_text(color="black"),
                       axis.title.x = element_text(lineheight = 1.2),
                       legend.position="none"
)

theme_mine <- function(base_size = 18, base_family = "Helvetica") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 18),
      strip.text.y = element_text(size = 18),
      axis.text.x = element_text(size=14, colour = "black"),
      axis.text.y = element_text(size=14,hjust=1, colour = "black"),
      axis.ticks =  element_line(colour = "black"), 
      axis.title.x= element_text(size=16, colour = "black"),
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

blankPlot <- ggplot()+geom_blank(aes(1,1)) +
  cowplot::theme_nothing()
de_bgam =  get_relative_contribs(bgam) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%"))

binary_vars = c("Genus_")
preds <- predict(bgam, type="iterms", se.fit=TRUE)
intercept <- attributes(preds)$constant
preds = preds %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
model_data_vir = bgam$model
gterms = attributes(terms(bgam))
ilfun <- bgam$family$linkinv
lfun <- bgam$family$linkfun
binary_terms = which(stri_detect_regex(names(model_data_vir), paste0("(", paste0(binary_vars, collapse="|"), ")")))
smooth_data_vir = model_data_vir[, -c(binary_terms, gterms$response, gterms$offset)]
smooth_ranges = dmap(smooth_data_vir, ~seq(min(.), max(.), length.out = 100))
binary_data = data.frame(Genus_Coleura=model_data_vir[, binary_terms])
binary_ranges = setNames(as.data.frame(diag(length(binary_terms))), names(model_data_vir)[binary_terms])
offset_name = stri_replace_first_regex(names(model_data_vir)[gterms$offset], "offset\\(([^\\)]+)\\)", "$1")
smooth_ranges = cbind(smooth_ranges,
                      setNames(as.data.frame(lapply(c(names(binary_data), offset_name), function(x) rep(0, nrow(smooth_ranges)))), c(names(binary_data), offset_name)))
binary_ranges = cbind(binary_ranges,
                      setNames(as.data.frame(lapply(c(names(smooth_data_vir), offset_name), function(x) rep(0, nrow(binary_ranges)))), c(names(smooth_data_vir), offset_name)))
smooth_preds <- predict(bgam, newdata=smooth_ranges, type="iterms", se.fit=TRUE)  %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
binary_preds <- predict(bgam, newdata=binary_ranges, type="iterms", se.fit=TRUE) %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
partials <- as.data.frame(lapply(1:ncol(preds$fit), function(cl) {
  x = lfun(model_data_vir[[gterms$response]]) - rowSums(preds$fit[,-cl]) - intercept
  x
}))
names(partials) = names(preds$fit)
relc = get_relative_contribs(bgam)
relc_order = relc[order(-relc$rel_deviance_explained),] 
relc_order = relc_order %>%
  filter(!grepl('Genus_', term))
smooth_data_vir = dplyr::select(smooth_data_vir, relc_order$term)
smooth_titles = paste(relc_order$term , ' (', round(relc_order$rel_deviance_explained *100, 2), '%)', sep='')
names(smooth_titles) = names(smooth_data_vir)

smooth_plots_vir = map(names(smooth_data_vir), function(smooth_term_vir) {
  pl =  ggplot() +
    geom_hline(yintercept = 0, size=0.1, col="grey50") +
    geom_point(mapping = aes(x=model_data_vir[[smooth_term_vir]], y = (partials[[smooth_term_vir]])),
               shape=21, fill="#0d17d6", col="black", alpha=0.25, size=1, stroke=0.1) +
    geom_ribbon(mapping = aes(x = smooth_ranges[[smooth_term_vir]],
                              ymin = (smooth_preds$fit[[smooth_term_vir]] - 2 * smooth_preds$se.fit[[smooth_term_vir]]),
                              ymax = (smooth_preds$fit[[smooth_term_vir]] + 2 * smooth_preds$se.fit[[smooth_term_vir]])),
                alpha = 0.5, fill="#20d226") +
    geom_line(mapping = aes(x = smooth_ranges[[smooth_term_vir]], y = (smooth_preds$fit[[smooth_term_vir]])), size=0.3) +
    xlab(smooth_titles[[smooth_term_vir]]) +
    scale_y_continuous(limits=c(-6,7), oob=scales::rescale_none) +
    theme_mine()
  return(pl)
})
smooth_plots_vir[[1]] = smooth_plots_vir[[1]] + ylab("strength of effect on\nviruses per host")
smooth_plots_vir[[2]] = smooth_plots_vir[[2]] + theme(axis.title.y= element_blank())
smooth_plots_vir[[3]] = smooth_plots_vir[[3]] + theme(axis.title.y= element_blank())
smooth_plots_vir[[4]] = smooth_plots_vir[[4]] + theme(axis.title.y= element_blank())

saveRDS(smooth_plots_vir, P("plot/smooth_plots_vir_0728.rds"))
# smooth_plots_vir = readRDS(P("plot/smooth_plots_vir.rds"))
p_smooth = c(smooth_plots_vir[1], smooth_plots_vir[2], smooth_plots_vir[3], smooth_plots_vir[4])

x11()
plot_grid(plotlist = p_smooth,
          nrow=1,
          align="h",
          rel_widths = c(1, 1, 1, 1))





# bin_vir_data = binary_preds %>% map(function(x) {
#   x = x[, stri_detect_regex(names(x), paste0("(", paste0(binary_vars, collapse="|"), ")"))]
#   n = names(model_data_vir[binary_terms])
#   tibble(response=x, variable=n)
# })
# bin_vir_data$fit$se = bin_vir_data$se.fit$response
# bin_vir_data = bin_vir_data$fit
# bin_vir_data$response = bin_vir_data$response
# bin_vir_data$labels = stri_replace_first_regex(bin_vir_data$variable, "hOrder", "")
# bin_vir_data$signif = summary(bgam)$s.table[stri_detect_regex(rownames(summary(bgam)$s.table), paste0("(", paste0(binary_vars, collapse="|"), ")")), "p-value"] < 0.05
# bin_vir_data = bin_vir_data %>%
#   arrange(desc(signif), response) %>%
#   mutate(no = 1:nrow(bin_vir_data))
# 
# bin_vir_partials = lapply(binary_terms, function(x) {
#   vals = partials[as.logical(model_data_vir[[x]]), x-1]
#   variable = names(model_data_vir)[x]
#   data_frame(variable=variable, partial=vals, no=bin_vir_data$no[bin_vir_data$variable == variable])
# }) %>% bind_rows
# 
# bin_vir_data = bin_vir_partials %>%
#   group_by(variable) %>%
#   summarize(minval = min(partial)) %>%
#   inner_join(bin_vir_data, by="variable") %>%
#   mutate(minval = pmin(minval, response - 2*se)) %>%
#   left_join(de_bgam, by=c('variable' = 'term'))
# 
# bin_plot_vir = ggplot() +
#   geom_hline(yintercept = 0, size=0.1, col="grey50") +
#   geom_point(data=bin_vir_partials, mapping=aes(x=no, y=(partial)), position=position_jitter(width=0.25),
#              shape=21, fill="#DA006C", col="black", alpha=0.25, size=1, stroke=0.1) +
#   geom_rect(data = bin_vir_data, mapping=aes(xmin = no - 0.35, xmax  = no + 0.35, ymin=(response-2*se), ymax=(response+2*se), fill=signif), alpha = 0.5) +
#   geom_segment(data = bin_vir_data, mapping=aes(x=no - 0.35, xend = no + 0.35, y=(response), yend=(response)), col="black", size=0.3) +
#   
#   geom_text(data = bin_vir_data, mapping=aes(x=no, y=(minval - 0.4), label = stri_trans_totitle(labels)),
#             color="black", family="Lato", size=2, angle =90, hjust=1, vjust =0.5) +
#   scale_fill_manual(breaks = c(TRUE, FALSE), values=c(viridis(5)[4]), "grey") +
#   scale_x_continuous(breaks = bin_vir_data$no, labels = stri_trans_totitle(bin_vir_data$labels)) +
#   scale_y_continuous(limits=c(-3.8,2.2), breaks=seq(-2,2, by=1), oob=scales::rescale_none, name="") +
#   theme_bw() + theme_mine()





