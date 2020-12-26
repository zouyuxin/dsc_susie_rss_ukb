library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
## Functions
plot_panel = function(dat, quantity, legend = TRUE) {
  p = ggplot(dat, aes_string(x="n_signal", y=quantity[1], fill="Method")) + 
    scale_color_manual("Method", values = c("SuSiE" = "#348ABD", "SuSiE-RSS" = "#A60628")) + 
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
add_z = c(FALSE, TRUE)
ldmethod = c('in_sample', 'ref_sample')
lambda = c(0, 1e-04, 0.01)
all.comb = expand.grid(add_z, ldmethod, lambda)
colnames(all.comb) = c('add_z', 'ldmethod', 'lambda')
all.comb = all.comb %>% filter(!(ldmethod == 'in_sample' & add_z == TRUE))

dat = readRDS(input)
dat$susie$addz[dat$susie$z_ld_weight == 0] = FALSE
dat$susie$addz[dat$susie$z_ld_weight == 0.002] = TRUE

for(case in 1:nrow(all.comb)){
  add_z = all.comb[case, 'add_z']
  ldmethod = all.comb[case, 'ldmethod']
  lambda = all.comb[case, 'lambda']
  output = paste0('susierss_ukb_default_cs_compare/susierss_ukb_default_pip_AZ',add_z,'_ld',ldmethod,
                  '_lamb', lambda,'_cs')
  susie_out = dat$susie %>% filter(method_susie_full == 'susie')
  if(lambda==0){
    susierss_out = dat$susie %>% filter(method_susie_full == 'susie_rss', 
                                        ld_method == ldmethod, addz == add_z)
  }else{
    susierss_out = dat$susie %>% filter(method_susie_full == 'susie_rss_lambda', 
                                        ld_method == ldmethod, addz == add_z,
                                        lamb == lambda)
  }
  data_sets = unique(susierss_out$dataset)
  n_signals = unique(susierss_out$n_signal)
  n_r = 1
  n_experiments = n_r * length(data_sets)
  result = NULL
  
  for (s in n_signals) {
    # I cannot find a good median tracker so do it stupid way: save all and take median later
    susierss_sizes = susie_sizes = vector()
    # do the same for mean ...
    susierss_avg_ld = susie_avg_ld = vector()
    susierss_valid = susie_valid = 0
    susierss_total = susie_total = 0
    susierss_expected = susie_expected = 0
    susierss_fail = susie_fail = 0
    # for debug
    susierss_overlap = susie_overlap = 0
    for (d in data_sets) {
      # should all be one row
      susie_row = susie_out %>% filter(n_signal == s, dataset == d)
      susierss_row = susierss_out %>% filter(n_signal == s, dataset == d)
      
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
    }
    rates = c(s, susie_total, susierss_total, susierss_expected,
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
  colnames(result) = c('n_signal', 'susie_discoveries', 'susierss_discoveries',
                       'expected_discoveries',
                       'susie_power', 'susierss_power', 
                       'susie_coverage', 'susierss_coverage', 
                       'susie_median_size', 'susierss_median_size', 
                       'susie_avg_ld', 'susierss_avg_ld', 
                       'susie_overlap', 'susierss_overlap', 
                       'susie_fail', 'susierss_fail')
  rownames(result) = as.character(result[,1])
  result = data.frame(result)
  result$susie_power_se = sqrt(result$susie_power * (1-result$susie_power) / result$expected_discoveries)
  result$susierss_power_se = sqrt(result$susierss_power * (1-result$susierss_power) / result$expected_discoveries)
  
  result$susie_coverage_se = sqrt(result$susie_coverage * (1-result$susie_coverage) / result$susie_discoveries)
  result$susierss_coverage_se = sqrt(result$susierss_coverage * (1-result$susierss_coverage) / result$susierss_discoveries)
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

  tmp = rbind(susie, susierss)
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
