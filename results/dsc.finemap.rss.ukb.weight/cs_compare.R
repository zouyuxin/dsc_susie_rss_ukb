library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
## Functions
plot_panel = function(dat, quantity, legend = TRUE) {
  p = ggplot(dat, aes_string(x="z_ld_weight", y=quantity[1], fill="Method")) + 
    scale_color_manual("Method", values = c("SuSiE" = "#348ABD", "SuSiE-RSS" = "#A60628", "SuSiE-RSS-lambda" = "#31a354")) + 
    geom_point(aes(colour=Method), position=position_dodge(.25), size=2.5)
  if (quantity[1] == 'power') p = p + geom_errorbar(aes(ymin=power-power_se, ymax=power+power_se, color=Method), width=.2, position=position_dodge(.25))
  if (quantity[1] == 'coverage') p = p + geom_errorbar(aes(ymin=coverage-coverage_se, ymax=coverage+coverage_se, color=Method), width=.2, position=position_dodge(.25)) + geom_hline(yintercept = 0.95, colour = 'gray') 
  p = p + labs(x = "z LD weight", y = "") + theme_cowplot() + background_grid(major = "x", minor = "none") + 
    ggtitle(quantity[2]) + theme(axis.title.x = element_text(size = 8), axis.text.x = element_text(size=5))
  if (!legend) p = p + theme(legend.position="none")
  return(p)
}
## parameters
input = 'susie_rss_ukb_weight_query.rds'
ldmethod = c('in_sample_Z', 'ref_sample_Z')
all.comb = expand.grid(ldmethod)
colnames(all.comb) = c('ldmethod')

dat = readRDS(input)

for(case in 1:nrow(all.comb)){
  ldmethod = all.comb[case, 'ldmethod']
  output = paste0('susie_rss_ukb_weight_cs_compare/susie_rss_ukb_weight_ld',ldmethod,'_cs')
  susie_out = dat$susie %>% filter(method_susie == 'susie')
  susierss_out = dat$susie %>% filter(method_susie == 'susie_rss', ld_method == ldmethod)
  
  data_sets = unique(susie_out$dataset)
  n_signals = unique(susie_out$n_signal)
  zldweight = unique(susierss_out$z_ld_weight)
  n_r = 1
  n_experiments = n_r * length(data_sets)
  result = NULL
  
  for (s in zldweight) {
    # I cannot find a good median tracker so do it stupid way: save all and take median later
    susie_sizes = susierss_sizes = vector()
    # do the same for mean ...
    susie_avg_ld = susierss_avg_ld = vector()
    susie_valid = susierss_valid = 0
    susie_total = susierss_total = 0
    susie_expected = susierss_expected = 0
    susie_fail = susierss_fail = 0
    # for debug
    susie_overlap = susierss_overlap = 0
    for (d in data_sets) {
      # should all be one row
      susie_row = susie_out %>% filter(dataset == d)
      susierss_row = susierss_out %>% filter(z_ld_weight == s, dataset == d)
      #
      susie_fail = susie_fail + sum(susie_row$converged[[1]]!=T, na.rm=T) + sum(is.na(susie_row$converged[[1]]))
      if(any(susie_row$converged[[1]])){
        converged = which(susie_row$converged[[1]] == TRUE)
        susie_total = susie_total + sum(unlist(susie_row$total)[converged])
        susie_valid = susie_valid + sum(unlist(susie_row$valid)[converged])
        susie_sizes = c(susie_sizes, unlist(susie_row$size)[converged])
        susie_avg_ld = c(susie_avg_ld, unlist(susie_row$avgr2)[converged])
        susie_overlap = susie_overlap + sum(unlist(susie_row$overlap)[converged])
        susie_expected = susie_expected + 2 * length(unlist(susie_row$total)[converged])
      }
      
      susierss_fail = susierss_fail + sum(susierss_row$converged[[1]]!=T, na.rm=T) + sum(is.na(susierss_row$converged[[1]]))
      if(any(susierss_row$converged[[1]])){
        converged = which(susierss_row$converged[[1]] == TRUE)
        susierss_total = susierss_total + sum(unlist(susierss_row$total)[converged])
        susierss_valid = susierss_valid + sum(unlist(susierss_row$valid)[converged])
        susierss_sizes = c(susierss_sizes, unlist(susierss_row$size)[converged])
        susierss_avg_ld = c(susierss_avg_ld, unlist(susierss_row$avgr2)[converged])
        susierss_overlap = susierss_overlap + sum(unlist(susierss_row$overlap)[converged])
        susierss_expected = susierss_expected + 2 * length(unlist(susierss_row$total)[converged])
      }
    }
    rates = c(s, susie_total, susierss_total, 
              susie_expected, susierss_expected, 
              susie_valid/susie_expected, susierss_valid/susierss_expected, 
              susie_valid/susie_total, susierss_valid/susierss_total, 
              median(susie_sizes), median(susierss_sizes), 
              mean(susie_avg_ld), mean(susierss_avg_ld,na.rm=T), 
              susie_overlap, susierss_overlap, 
              susie_fail, susierss_fail)
    if (is.null(result)) {
      result = rates
    } else {
      result = rbind(result, rates)
    }
  }
  colnames(result) = c('z_ld_weight', 'susie_discoveries', 'rss_discoveries', 
                       'susie_expected', 'rss_expected', 
                       'susie_power', 'rss_power',
                       'susie_coverage', 'rss_coverage',
                       'susie_median_size', 'rss_median_size',
                       'susie_avg_ld', 'rss_avg_ld',
                       'susie_overlap', 'rss_overlap',
                       'susie_fail', 'rss_fail')
  rownames(result) = as.character(result[,1])
  result = data.frame(result)
  result$susie_power_se = sqrt(result$susie_power * (1-result$susie_power) / result$susie_expected)
  result$rss_power_se = sqrt(result$rss_power * (1-result$rss_power) / result$rss_expected)
  
  result$susie_coverage_se = sqrt(result$susie_coverage * (1-result$susie_coverage) / result$susie_discoveries)
  result$rss_coverage_se = sqrt(result$rss_coverage * (1-result$rss_coverage) / result$rss_discoveries)
  saveRDS(result, paste0(output, '.rds'))
  
  susierss = cbind(result[, c("z_ld_weight", "rss_expected", "rss_power", "rss_coverage", 
                              "rss_power_se", "rss_coverage_se", "rss_median_size", "rss_avg_ld")], 
                   'SuSiE-RSS')
  colnames(susierss) = c("z_ld_weight", "expected", "power", "coverage", "power_se", "coverage_se", 
                         "median_size", "avg_ld", "Method")
  susie = cbind(result[, c("z_ld_weight", "susie_expected", "susie_power", "susie_coverage", 
                           "susie_power_se", "susie_coverage_se", "susie_median_size", "susie_avg_ld")], 
                   'SuSiE')
  colnames(susie) = c("z_ld_weight", "expected", "power", "coverage", "power_se", "coverage_se", 
                      "median_size", "avg_ld", "Method")
  tmp = rbind(susie, susierss)
  tmp$z_ld_weight = as.factor(tmp$z_ld_weight)
  p1 = plot_panel(tmp, c('coverage', 'coverage'), legend=F)
  p2 = plot_panel(tmp, c('power', 'power'), legend=F)
  p3 = plot_panel(tmp, c('median_size', 'median number of variables'), legend=F)
  p4 = plot_panel(tmp, c('avg_ld', 'average r2'), legend=T)
  pdf(paste0(output, '_plots.pdf'), width=14, height=3)
  grid.arrange(p1,p2,p3,p4, ncol=4, widths=c(3,3,3,5))
  dev.off()
  system(paste0("convert -flatten -density 120 ", paste0(output, '_plots.pdf'), " ", paste0(output, '_plots.png')))
}
