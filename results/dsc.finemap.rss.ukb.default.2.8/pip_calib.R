library(ggplot2)
library(cowplot)
library(dplyr)
## Functions
get_cali = function(alist, col) {
  res = alist[[1]][[col]]
  if(length(alist) > 1){
    for (i in 2:length(alist)) {
      if (!is.null(alist[[i]][[col]])) res = res + alist[[i]][[col]]
    }
  }
  res[,c(1,2)] = res[,c(1,2)] / res[,3]
  return(res[-1,])
}
dot_plot = function(dataframe) {
  ggplot(dataframe, aes(x=mean_pip, y=observed_freq)) +
    geom_errorbar(aes(ymin=observed_freq-se, ymax=observed_freq+se), colour="gray", size = 0.2, width=.01) +
    geom_point(size=1.5, shape=21, fill="#002b36") + # 21 is filled circle
    xlab("Mean PIP") +
    ylab("Observed frequency") +
    coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
    geom_abline(slope=1,intercept=0,colour='red', size=0.2) +
    ggtitle(name) +
    expand_limits(y=0) +                        # Expand y range
    theme_cowplot()
}

## parameters
initial = c(NA, 'oracle', 'lasso')
add_z = c(FALSE, TRUE)
ldmethod = c('in_sample', 'ref_sample')
estimate_resid = c(FALSE, TRUE)
lambda = c(0, 1e-04, 0.1)
all.comb = expand.grid(initial, add_z, ldmethod, estimate_resid, lambda)
colnames(all.comb) = c('initial', 'add_z', 'ldmethod', 'estimate_resid', 'lambda')
all.comb = all.comb %>% filter(!(ldmethod == 'in_sample' & add_z == TRUE))

bin_size = 20
pip_cutoff = 0

for (case in 1:nrow(all.comb)){
  initial = all.comb[case, 'initial']
  add_z = all.comb[case, 'add_z']
  ldmethod = all.comb[case, 'ldmethod']
  estimate_resid = all.comb[case, 'estimate_resid']
  lambda = all.comb[case, 'lambda']
  input = paste0('susie_rss_ukb_default_pip_extraction/susie_rss_ukb_default_pip_init',initial,'_AZ',add_z,'_ld',ldmethod,'_ER',
                 estimate_resid, '_lamb', lambda,'.rds')
  output_name = paste0('susie_rss_ukb_default_pip_calibration/susie_rss_ukb_default_pip_init',initial,'_AZ',add_z,'_ld',ldmethod,'_ER',
                       estimate_resid, '_lamb', lambda,'_pipcali')
  output = paste0(output_name,'.rds')
  
  dat = readRDS(input)
  pip_cali = list()
  bins = cbind(seq(1:bin_size)/bin_size-1/bin_size, seq(1:bin_size)/bin_size)
  for (s in 1:length(dat)) {
    res = dat[[as.character(s)]]
    pip_cali[[as.character(s)]] = list()
    for (name in rev(colnames(res))[-1]) {
      for (i in 1:nrow(bins)) {
        tmp = res[which(res[[name]] > bins[i,1] & res[[name]] < bins[i,2]),]
        if (is.null(pip_cali[[as.character(s)]][[name]])) pip_cali[[as.character(s)]][[name]] = c(sum(tmp[[name]]), sum(tmp$is_signal), length(tmp$is_signal))
        else pip_cali[[as.character(s)]][[name]] = rbind(pip_cali[[as.character(s)]][[name]], c(sum(tmp[[name]]), sum(tmp$is_signal), length(tmp$is_signal)))
      }
      pip_cali[[as.character(s)]][[name]][which(is.na(pip_cali[[as.character(s)]][[name]]))] = 0
    }
  }
  saveRDS(list("SuSiE"=get_cali(pip_cali, 'susie'),
               "SuSiE-RSS"=get_cali(pip_cali, 'susierss'),
               "CAVIAR"=get_cali(pip_cali, 'caviar'),
               "FINEMAP"=get_cali(pip_cali, 'finemap'),
               "FINEMAPv3"=get_cali(pip_cali, 'finemapv3')),
          output)
  
  dat = readRDS(output)
  idx = 0
  for (name in names(dat)) {
    idx = idx + 1
    dat[[name]][,3] = sqrt(dat[[name]][,2] * (1 - dat[[name]][,2]) / dat[[name]][,3]) * 2
    dat[[name]] = as.data.frame(dat[[name]])
    colnames(dat[[name]]) = c("mean_pip", "observed_freq", "se")
    pdf(paste0(output_name, '_' , idx, '.pdf'), width=3, height=3, pointsize=16)
    print(dot_plot(dat[[name]]))
    dev.off()
    system(paste0("convert -flatten -density 120 ", output_name, '_' , idx, '.pdf', " ", output_name, '_' , idx, '.png'))
  }
  files = paste0(output_name, '_', seq(1:idx), '.png')
  cmd = paste('convert +append', paste(files, collapse=" "), paste0(output_name, '.png'))
  system(cmd)
  # system(paste('rm -f', paste(files, collapse=" ")))
  
}
