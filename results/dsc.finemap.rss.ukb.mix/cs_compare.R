library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
## Functions
plot_panel = function(dat, quantity, legend = TRUE) {
  p = ggplot(dat, aes_string(x="n_signal", y=quantity[1], fill="Method")) + 
    scale_color_manual("Method", values = c("lm in sample" = "#A60628",
                                            "mix in sample" = "#00CC00",
                                            "lm ref sample" = "#00FFFF",
                                            "mix ref sample" = "#FF00FF")) + 
    geom_point(aes(colour=Method), position=position_dodge(.25), size=2.5)
  if (quantity[1] == 'power') p = p + geom_errorbar(aes(ymin=power-power_se, ymax=power+power_se, color=Method), width=.2, position=position_dodge(.25))
  if (quantity[1] == 'coverage') p = p + geom_errorbar(aes(ymin=coverage-coverage_se, ymax=coverage+coverage_se, color=Method), width=.2, position=position_dodge(.25)) + geom_hline(yintercept = 0.95, colour = 'gray') 
  p = p + labs(x = "number of effect variables", y = "") + theme_cowplot() + background_grid(major = "x", minor = "none") + 
    ggtitle(quantity[2]) + theme(axis.title.x = element_text(size = 8))
  if (!legend) p = p + theme(legend.position="none")
  return(p)
}
## parameters
input = 'susie_rss_ukb_mix_query.rds'
add_z = c(FALSE, TRUE)
estimate_resid = c(FALSE, TRUE)
all.comb = expand.grid(add_z, estimate_resid)
colnames(all.comb) = c('add_z', 'estimate_resid')

dat = readRDS(input)

for(case in 1:nrow(all.comb)){
  add_z = all.comb[case, 'add_z']
  estimate_resid = all.comb[case, 'estimate_resid']
  output = paste0('susie_rss_ukb_mix_cs_compare/susie_rss_ukb_mix_pip_AZ',add_z,'_ER',
                  estimate_resid, '_cs')
  
  susierss_lm_insample = dat$susierss_lm_insample %>% filter(addz == add_z, estimate_residual_variance == estimate_resid)
  susierss_mix_insample = dat$susierss_mix_insample %>% filter(addz == add_z, estimate_residual_variance == estimate_resid)
  susierss_lm_refsample = dat$susierss_lm_refsample %>% filter(addz == add_z, estimate_residual_variance == estimate_resid)
  susierss_mix_refsample = dat$susierss_mix_refsample %>% filter(addz == add_z, estimate_residual_variance == estimate_resid)  

  data_sets = unique(susierss_mix_insample$dataset)
  n_signals = unique(susierss_mix_insample$n_signal)
  n_r = 1
  n_experiments = n_r * length(data_sets)
  result = NULL
  
  for (s in n_signals) {
    # I cannot find a good median tracker so do it stupid way: save all and take median later
    susierss_lm_insample_sizes = susierss_mix_insample_sizes = susierss_lm_refsample_sizes = susierss_mix_refsample_sizes = vector()
    # do the same for mean ...
    susierss_lm_insample_avg_ld = susierss_mix_insample_avg_ld =susierss_lm_refsample_avg_ld =susierss_mix_refsample_avg_ld = vector()
    susierss_lm_insample_valid = susierss_mix_insample_valid = susierss_lm_refsample_valid = susierss_mix_refsample_valid = 0
    susierss_lm_insample_total = susierss_mix_insample_total = susierss_lm_refsample_total = susierss_mix_refsample_total = 0
    susierss_lm_insample_expected = susierss_mix_insample_expected = susierss_lm_refsample_expected = susierss_mix_refsample_expected = 0
    susierss_lm_insample_fail = susierss_mix_insample_fail = susierss_lm_refsample_fail = susierss_mix_refsample_fail = 0
    # for debug
    susierss_lm_insample_overlap = susierss_mix_insample_overlap = susierss_lm_refsample_overlap = susierss_mix_refsample_overlap = 0
    for (d in data_sets) {
      # should all be one row
      susierss_lm_insample_row = susierss_lm_insample %>% filter(n_signal == s, dataset == d)
      susierss_mix_insample_row = susierss_mix_insample %>% filter(n_signal == s, dataset == d)
      susierss_lm_refsample_row = susierss_lm_refsample %>% filter(n_signal == s, dataset == d)
      susierss_mix_refsample_row = susierss_mix_refsample %>% filter(n_signal == s, dataset == d)
 
      susierss_lm_insample_fail = susierss_lm_insample_fail + sum(susierss_lm_insample_row$converged[[1]]!=T, na.rm=T) + sum(is.na(susierss_lm_insample_row$converged[[1]]))
      if(any(!is.na(susierss_lm_insample_row$converged[[1]]) & any(susierss_lm_insample_row$converged[[1]]))){
        converged = which(susierss_lm_insample_row$converged[[1]] == TRUE)
        susierss_lm_insample_total = susierss_lm_insample_total + sum(unlist(susierss_lm_insample_row$total)[converged])
        susierss_lm_insample_valid = susierss_lm_insample_valid + sum(unlist(susierss_lm_insample_row$valid)[converged])
        susierss_lm_insample_sizes = c(susierss_lm_insample_sizes, unlist(susierss_lm_insample_row$size)[converged])
        susierss_lm_insample_avg_ld = c(susierss_lm_insample_avg_ld, unlist(susierss_lm_insample_row$avgr2)[converged])
        susierss_lm_insample_overlap = susierss_lm_insample_overlap + sum(unlist(susierss_lm_insample_row$overlap)[converged])
        susierss_lm_insample_expected = susierss_lm_insample_expected + s * length(unlist(susierss_lm_insample_row$total)[converged])
      }
      
      susierss_mix_insample_fail = susierss_mix_insample_fail + sum(susierss_mix_insample_row$converged[[1]]!=T, na.rm=T) + sum(is.na(susierss_mix_insample_row$converged[[1]]))
      if(any(!is.na(susierss_mix_insample_row$converged[[1]]) & any(susierss_mix_insample_row$converged[[1]]))){
        converged = which(susierss_mix_insample_row$converged[[1]] == TRUE)
        susierss_mix_insample_total = susierss_mix_insample_total + sum(unlist(susierss_mix_insample_row$total)[converged])
        susierss_mix_insample_valid = susierss_mix_insample_valid + sum(unlist(susierss_mix_insample_row$valid)[converged])
        susierss_mix_insample_sizes = c(susierss_mix_insample_sizes, unlist(susierss_mix_insample_row$size)[converged])
        susierss_mix_insample_avg_ld = c(susierss_mix_insample_avg_ld, unlist(susierss_mix_insample_row$avgr2)[converged])
        susierss_mix_insample_overlap = susierss_mix_insample_overlap + sum(unlist(susierss_mix_insample_row$overlap)[converged])
        susierss_mix_insample_expected = susierss_mix_insample_expected + s * length(unlist(susierss_mix_insample_row$total)[converged])
      }

      susierss_lm_refsample_fail = susierss_lm_refsample_fail + sum(susierss_lm_refsample_row$converged[[1]]!=T, na.rm=T) + sum(is.na(susierss_lm_refsample_row$converged[[1]]))
      if(any(!is.na(susierss_lm_refsample_row$converged[[1]]) & any(susierss_lm_refsample_row$converged[[1]]))){
        converged = which(susierss_lm_refsample_row$converged[[1]] == TRUE)
        susierss_lm_refsample_total = susierss_lm_refsample_total + sum(unlist(susierss_lm_refsample_row$total)[converged])
        susierss_lm_refsample_valid = susierss_lm_refsample_valid + sum(unlist(susierss_lm_refsample_row$valid)[converged])
        susierss_lm_refsample_sizes = c(susierss_lm_refsample_sizes, unlist(susierss_lm_refsample_row$size)[converged])
        susierss_lm_refsample_avg_ld = c(susierss_lm_refsample_avg_ld, unlist(susierss_lm_refsample_row$avgr2)[converged])
        susierss_lm_refsample_overlap = susierss_lm_refsample_overlap + sum(unlist(susierss_lm_refsample_row$overlap)[converged])
        susierss_lm_refsample_expected = susierss_lm_refsample_expected + s * length(unlist(susierss_lm_refsample_row$total)[converged])
      }
      
      susierss_mix_refsample_fail = susierss_mix_refsample_fail + sum(susierss_mix_refsample_row$converged[[1]]!=T, na.rm=T) + sum(is.na(susierss_mix_refsample_row$converged[[1]]))
      if(any(!is.na(susierss_mix_refsample_row$converged[[1]]) & any(susierss_mix_refsample_row$converged[[1]]))){
        converged = which(susierss_mix_refsample_row$converged[[1]] == TRUE)
        susierss_mix_refsample_total = susierss_mix_refsample_total + sum(unlist(susierss_mix_refsample_row$total)[converged])
        susierss_mix_refsample_valid = susierss_mix_refsample_valid + sum(unlist(susierss_mix_refsample_row$valid)[converged])
        susierss_mix_refsample_sizes = c(susierss_mix_refsample_sizes, unlist(susierss_mix_refsample_row$size)[converged])
        susierss_mix_refsample_avg_ld = c(susierss_mix_refsample_avg_ld, unlist(susierss_mix_refsample_row$avgr2)[converged])
        susierss_mix_refsample_overlap = susierss_mix_refsample_overlap + sum(unlist(susierss_mix_refsample_row$overlap)[converged])
        susierss_mix_refsample_expected = susierss_mix_refsample_expected + s * length(unlist(susierss_mix_refsample_row$total)[converged])
      }
    }
    rates = c(s, susierss_lm_insample_total, susierss_mix_insample_total, susierss_lm_refsample_total, susierss_mix_refsample_total,
              susierss_lm_insample_expected,
              susierss_lm_insample_valid/susierss_lm_insample_expected, susierss_mix_insample_valid/susierss_mix_insample_expected,
              susierss_lm_refsample_valid/susierss_lm_refsample_expected, susierss_mix_refsample_valid/susierss_mix_refsample_expected,
              susierss_lm_insample_valid/susierss_lm_insample_total, susierss_mix_insample_valid/susierss_mix_insample_total,
              susierss_lm_refsample_valid/susierss_lm_refsample_total, susierss_mix_refsample_valid/susierss_mix_refsample_total, 
              median(susierss_lm_insample_sizes),median(susierss_mix_insample_sizes),median(susierss_lm_refsample_sizes),median(susierss_mix_refsample_sizes),
              mean(susierss_lm_insample_avg_ld,na.rm=T), mean(susierss_mix_insample_avg_ld,na.rm=T), mean(susierss_lm_refsample_avg_ld,na.rm=T), mean(susierss_mix_refsample_avg_ld,na.rm=T),
              susierss_lm_insample_overlap, susierss_mix_insample_overlap, susierss_lm_refsample_overlap, susierss_mix_refsample_overlap,
              susierss_lm_insample_fail, susierss_mix_insample_fail, susierss_lm_refsample_fail, susierss_mix_refsample_fail)
    if (is.null(result)) {
      result = rates
    } else {
      result = rbind(result, rates)
    }
  }
  result = matrix(result, length(n_signals), 30)
  colnames(result) = c('n_signal', 'susierss_lm_insample_discoveries', 'susierss_mix_insample_discoveries', 'susierss_lm_refsample_discoveries', 'susierss_mix_refsample_discoveries',
                       'expected_discoveries',
                       'susierss_lm_insample_power', 'susierss_mix_insample_power', 'susierss_lm_refsample_power','susierss_mix_refsample_power',
                       'susierss_lm_insample_coverage','susierss_mix_insample_coverage','susierss_lm_refsample_coverage','susierss_mix_refsample_coverage',
                       'susierss_lm_insample_median_size','susierss_mix_insample_median_size','susierss_lm_refsample_median_size','susierss_mix_refsample_median_size',
                       'susierss_lm_insample_avg_ld', 'susierss_mix_insample_avg_ld', 'susierss_lm_refsample_avg_ld', 'susierss_mix_refsample_avg_ld',
                       'susierss_lm_insample_overlap', 'susierss_mix_insample_overlap', 'susierss_lm_refsample_overlap', 'susierss_mix_refsample_overlap',
                       'susierss_lm_insample_fail', 'susierss_mix_insample_fail', 'susierss_lm_refsample_fail', 'susierss_mix_refsample_fail')
  rownames(result) = as.character(result[,1])
  result = data.frame(result)
  result$susierss_lm_insample_power_se = sqrt(result$susierss_lm_insample_power * (1-result$susierss_lm_insample_power) / result$expected_discoveries)
  result$susierss_mix_insample_power_se = sqrt(result$susierss_mix_insample_power * (1-result$susierss_mix_insample_power) / result$expected_discoveries)
  result$susierss_lm_refsample_power_se = sqrt(result$susierss_lm_refsample_power * (1-result$susierss_lm_refsample_power) / result$expected_discoveries)
  result$susierss_mix_refsample_power_se = sqrt(result$susierss_mix_refsample_power * (1-result$susierss_mix_refsample_power) / result$expected_discoveries)

  result$susierss_lm_insample_coverage_se = sqrt(result$susierss_lm_insample_coverage * (1-result$susierss_lm_insample_coverage) / result$susierss_lm_insample_discoveries)
  result$susierss_mix_insample_coverage_se = sqrt(result$susierss_mix_insample_coverage * (1-result$susierss_mix_insample_coverage) / result$susierss_mix_insample_discoveries)
  result$susierss_lm_refsample_coverage_se = sqrt(result$susierss_lm_refsample_coverage * (1-result$susierss_lm_refsample_coverage) / result$susierss_lm_refsample_discoveries)
  result$susierss_mix_refsample_coverage_se = sqrt(result$susierss_mix_refsample_coverage * (1-result$susierss_mix_refsample_coverage) / result$susierss_mix_refsample_discoveries)
  saveRDS(result, paste0(output, '.rds'))
  
  susierss_lm_insample = cbind(result[, c("n_signal", "expected_discoveries", "susierss_lm_insample_power", "susierss_lm_insample_coverage", 
                              "susierss_lm_insample_power_se", "susierss_lm_insample_coverage_se", "susierss_lm_insample_median_size", "susierss_lm_insample_avg_ld")], 
                   'lm in sample')
  colnames(susierss_lm_insample) = c("n_signal", "expected_discoveries", "power", "coverage", "power_se", "coverage_se", 
                         "median_size", "avg_ld", "Method")

  susierss_mix_insample = cbind(result[, c("n_signal", "expected_discoveries", "susierss_mix_insample_power", "susierss_mix_insample_coverage",
                              "susierss_mix_insample_power_se", "susierss_mix_insample_coverage_se", "susierss_mix_insample_median_size", "susierss_mix_insample_avg_ld")],
                   'mix in sample')
  colnames(susierss_mix_insample) = c("n_signal", "expected_discoveries", "power", "coverage", "power_se", "coverage_se",
                         "median_size", "avg_ld", "Method")

  susierss_lm_refsample = cbind(result[, c("n_signal", "expected_discoveries", "susierss_lm_refsample_power", "susierss_lm_refsample_coverage",
                              "susierss_lm_refsample_power_se", "susierss_lm_refsample_coverage_se", "susierss_lm_refsample_median_size", "susierss_lm_refsample_avg_ld")],
                   'lm ref sample')
  colnames(susierss_lm_refsample) = c("n_signal", "expected_discoveries", "power", "coverage", "power_se", "coverage_se",
                         "median_size", "avg_ld", "Method")
  
  susierss_mix_refsample = cbind(result[, c("n_signal", "expected_discoveries", "susierss_mix_refsample_power", "susierss_mix_refsample_coverage",
                              "susierss_mix_refsample_power_se", "susierss_mix_refsample_coverage_se", "susierss_mix_refsample_median_size", "susierss_mix_refsample_avg_ld")],
                   'mix ref sample')
  colnames(susierss_mix_refsample) = c("n_signal", "expected_discoveries", "power", "coverage", "power_se", "coverage_se",
                         "median_size", "avg_ld", "Method")

  tmp = rbind(susierss_lm_insample, susierss_mix_insample, susierss_lm_refsample, susierss_mix_refsample)
  tmp$n_signal = as.factor(tmp$n_signal)
  p1 = plot_panel(tmp, c('coverage', 'coverage'), legend=F)
  p2 = plot_panel(tmp, c('power', 'power'), legend=F)
  p3 = plot_panel(tmp, c('median_size', 'median number of variables'), legend=F)
  p4 = plot_panel(tmp, c('avg_ld', 'average r2'), legend=T)
  pdf(paste0(output, '_plots.pdf'), width=14, height=3)
  grid.arrange(p1,p2,p3,p4, ncol=4, widths=c(3,3,3,5))
  dev.off()
  system(paste0("convert -flatten -density 120 ", paste0(output, '_plots.pdf'), " ", paste0(output, '_plots.png')))
}
