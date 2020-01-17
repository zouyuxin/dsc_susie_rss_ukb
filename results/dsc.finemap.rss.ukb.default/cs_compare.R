library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
## Functions
plot_panel = function(dat, quantity, legend = TRUE) {
  p = ggplot(dat, aes_string(x="n_signal", y=quantity[1], fill="Method")) + 
    scale_color_manual("Method", values = c("SuSiE" = "#348ABD", "SuSiE-RSS" = "#A60628", "FINEMAPv3" = "#31a354")) + 
    geom_point(aes(colour=Method), position=position_dodge(.25), size=2.5)
  if (quantity[1] == 'power') p = p + geom_errorbar(aes(ymin=power-power_se, ymax=power+power_se, color=Method), width=.2, position=position_dodge(.25))
  if (quantity[1] == 'coverage') p = p + geom_errorbar(aes(ymin=coverage-coverage_se, ymax=coverage+coverage_se, color=Method), width=.2, position=position_dodge(.25)) + geom_hline(yintercept = 0.95, colour = 'gray') 
  p = p + labs(x = "number of effect variables", y = "") + theme_cowplot() + background_grid(major = "x", minor = "none") + 
    ggtitle(quantity[2]) + theme(axis.title.x = element_text(size = 8))
  if (!legend) p = p + theme(legend.position="none")
  return(p)
}
## parameters
input = 'susie_rss_ukb_default_query.rds'
addz = c(FALSE, TRUE)
ldmethod = c('in_sample', 'ref_sample')
estimate_resid = c(FALSE, TRUE)
lambda = c(0, 1e-04, 0.1)
all.comb = expand.grid(addz, ldmethod, estimate_resid, lambda)
colnames(all.comb) = c('addz', 'ldmethod', 'estimate_resid', 'lambda')
all.comb = all.comb %>% filter(!(ldmethod == 'in_sample' & addz == TRUE))

dat = readRDS(input)

for(case in 1:nrow(all.comb)){
  addz = all.comb[case, 'addz']
  ldmethod = all.comb[case, 'ldmethod']
  estimate_resid = all.comb[case, 'estimate_resid']
  lambda = all.comb[case, 'lambda']
  output = paste0('susie_rss_ukb_default_cs_compare/susie_rss_ukb_default_pip_AZ',addz,'_ld',ldmethod,'_ER',
                  estimate_resid, '_lamb', lambda,'_cs')
  susie_out = dat$susie %>% filter(method_susie == 'susie', 
                                   estimate_residual_variance == estimate_resid)
  susierss_out = dat$susie %>% filter(method_susie != 'susie', 
                                      ld_method == ldmethod, add_z == addz,
                                      estimate_residual_variance == estimate_resid,
                                      lamb == lambda)
  finemapv3_out = dat$finemapv3 %>% filter(add_z == addz, ld_method == ldmethod)
  data_sets = unique(susierss_out$dataset)
  n_signals = unique(susierss_out$n_signal)
  n_r = 1
  n_experiments = n_r * length(data_sets)
  result = NULL
  
  for (s in n_signals) {
    # I cannot find a good median tracker so do it stupid way: save all and take median later
    susierss_sizes = susie_sizes = finemapv3_sizes = vector()
    # do the same for mean ...
    susierss_avg_ld = susie_avg_ld = vector()
    susierss_valid = susie_valid = finemapv3_valid = 0
    susierss_total = susie_total = finemapv3_total = 0
    susierss_expected = susie_expected = finemapv3_expected = 0
    susierss_fail = susie_fail = finemapv3_fail = 0
    # for debug
    susierss_overlap = susie_overlap = 0
    for (d in data_sets) {
      # should all be one row
      susie_row = susie_out %>% filter(n_signal == s, dataset == d)
      susierss_row = susierss_out %>% filter(n_signal == s, dataset == d, L == 10)
      
      susierss_fail = susierss_fail + sum(susierss_row$converged[[1]]!=T, na.rm=T) + sum(is.na(susierss_row$converged[[1]]))
      if(any(susierss_row$converged[[1]])){
        converged = which(susierss_row$converged[[1]] == TRUE)
        susierss_total = susierss_total + sum(unlist(susierss_row$total)[converged])
        susierss_valid = susierss_valid + sum(unlist(susierss_row$valid)[converged])
        susierss_sizes = c(susierss_sizes, unlist(susierss_row$size)[converged])
        susierss_avg_ld = c(susierss_avg_ld, unlist(susierss_row$avgr2)[converged])
        susierss_overlap = susierss_overlap + sum(unlist(susierss_row$overlap)[converged])
        susierss_expected = susierss_expected + s * length(unlist(susierss_row$total)[converged])
      }
      #
      susie_fail = susie_fail + sum(susie_row$converged[[1]]!=T, na.rm=T) + sum(is.na(susie_row$converged[[1]]))
      if(any(susie_row$converged[[1]])){
        converged = which(susie_row$converged[[1]] == TRUE)
        susie_total = susie_total + sum(unlist(susie_row$total)[converged])
        susie_valid = susie_valid + sum(unlist(susie_row$valid)[converged])
        susie_sizes = c(susie_sizes, unlist(susie_row$size)[converged])
        susie_avg_ld = c(susie_avg_ld, unlist(susie_row$avgr2)[converged])
        susie_overlap = susie_overlap + sum(unlist(susie_row$overlap)[converged])
        susie_expected = susie_expected + s * length(unlist(susie_row$total)[converged])
      }
      if(s <= 3){
        finemapv3_row = finemapv3_out %>% filter(n_signal == s, dataset == d, args == paste('--n-causal-snps', s))
        finemapv3_fail = finemapv3_fail + sum(is.na(finemapv3_row$total))
        if(any(!is.na(finemapv3_row$total))){
          converged = which(!is.na(finemapv3_row$total))
          finemapv3_total = finemapv3_total + sum(unlist(finemapv3_row$total)[converged])
          finemapv3_valid = finemapv3_valid + sum(unlist(finemapv3_row$valid)[converged])
          finemapv3_sizes = c(finemapv3_sizes, unlist(finemapv3_row$size)[converged])
          finemapv3_expected = finemapv3_expected + s * length(unlist(finemapv3_row$total)[converged])
        }
      }
    }
    rates = c(s, susie_total, susierss_total, finemapv3_total, susierss_expected,
              susie_valid/susie_expected, susierss_valid/susierss_expected, finemapv3_valid/finemapv3_expected,
              susie_valid/susie_total, susierss_valid/susierss_total, finemapv3_valid/finemapv3_total,
              median(susie_sizes), median(susierss_sizes), median(finemapv3_sizes),
              mean(susie_avg_ld), mean(susierss_avg_ld,na.rm=T), 
              susie_overlap, susierss_overlap, 
              susie_fail, susierss_fail, finemapv3_fail)
    if (is.null(result)) {
      result = rates
    } else {
      result = rbind(result, rates)
    }
  }
  colnames(result) = c('n_signal', 'susie_discoveries', 'susierss_discoveries', 'finemapv3_discoveries',
                       'expected_discoveries',
                       'susie_power', 'susierss_power', 'finemapv3_power',
                       'susie_coverage', 'susierss_coverage', 'finemapv3_coverage',
                       'susie_median_size', 'susierss_median_size', 'finemapv3_median_size',
                       'susie_avg_ld', 'susierss_avg_ld', 
                       'susie_overlap', 'susierss_overlap', 
                       'susie_fail', 'susierss_fail', 'finemapv3_fail')
  rownames(result) = as.character(result[,1])
  result = data.frame(result)
  result$susie_power_se = sqrt(result$susie_power * (1-result$susie_power) / result$expected_discoveries)
  result$susierss_power_se = sqrt(result$susierss_power * (1-result$susierss_power) / result$expected_discoveries)
  result$finemapv3_power_se = sqrt(result$finemapv3_power * (1-result$finemapv3_power) / result$expected_discoveries)
  result$finemapv3_power_se[is.nan(result$finemapv3_power_se)] = NA
  result$finemapv3_power[is.nan(result$finemapv3_power)] = NA
  
  result$susie_coverage_se = sqrt(result$susie_coverage * (1-result$susie_coverage) / result$susie_discoveries)
  result$susierss_coverage_se = sqrt(result$susierss_coverage * (1-result$susierss_coverage) / result$susierss_discoveries)
  result$finemapv3_coverage_se = sqrt(result$finemapv3_coverage * (1-result$finemapv3_coverage) / result$finemapv3_discoveries)
  result$finemapv3_coverage_se[is.nan(result$finemapv3_coverage_se)] = NA
  result$finemapv3_coverage[is.nan(result$finemapv3_coverage)] = NA
  saveRDS(result, paste0(output, '.rds'))
  
  susierss = cbind(result[, c("n_signal", "expected_discoveries", "susierss_power", "susierss_coverage", 
                              "susierss_power_se", "susierss_coverage_se", "susierss_median_size", "susierss_avg_ld")], 
                   'SuSiE-RSS')
  colnames(susierss) = c("n_signal", "expected_discoveries", "power", "coverage", "power_se", "coverage_se", 
                         "median_size", "avg_ld", "Method")
  susie = cbind(result[, c("n_signal", "expected_discoveries", "susie_power", "susie_coverage", 
                           "susie_power_se", "susie_coverage_se", "susie_median_size", "susie_avg_ld")], 
                   'SuSiE')
  colnames(susie) = c("n_signal", "expected_discoveries", "power", "coverage", "power_se", "coverage_se", 
                      "median_size", "avg_ld", "Method")
  finemapv3 = cbind(result[, c("n_signal", "expected_discoveries", "finemapv3_power", "finemapv3_coverage", 
                               "finemapv3_power_se", "finemapv3_coverage_se", "finemapv3_median_size")], 
                    'finemapv3_avg_ld','FINEMAPv3')
  finemapv3[,8] = NA
  colnames(finemapv3) = c("n_signal", "expected_discoveries", "power", "coverage", "power_se", "coverage_se", 
                          "median_size", "avg_ld","Method")
  tmp = rbind(susie, susierss, finemapv3)
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
