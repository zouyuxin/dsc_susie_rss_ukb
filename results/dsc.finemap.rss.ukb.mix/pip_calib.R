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
add_z = c(FALSE, TRUE)
estimate_resid = c(FALSE, TRUE)
all.comb = expand.grid(add_z, estimate_resid)
colnames(all.comb) = c('add_z', 'estimate_resid')

bin_size = 20
pip_cutoff = 0

for (case in 1:nrow(all.comb)){
  add_z = all.comb[case, 'add_z']
  estimate_resid = all.comb[case, 'estimate_resid']
  input = paste0('susie_rss_ukb_mix_pip_extraction/susie_rss_ukb_mix_pip_AZ',add_z,'_ER',
                 estimate_resid, '.rds')
  output_name = paste0('susie_rss_ukb_mix_pip_calibration/susie_rss_ukb_mix_pip_AZ',add_z,'_ER',
                       estimate_resid,'_pipcali')
  output = paste0(output_name,'.rds')
  
  dat = readRDS(input)
  pip_cali = list()
  bins = cbind(seq(1:bin_size)/bin_size-1/bin_size, seq(1:bin_size)/bin_size)
  for (s in 1:length(dat)) {
    res = dat[[s]]
    pip_cali[[s]] = list()
    for (name in rev(colnames(res))[-1]) {
      for (i in 1:nrow(bins)) {
        tmp = res[which(res[[name]] > bins[i,1] & res[[name]] < bins[i,2]),]
        if (is.null(pip_cali[[s]][[name]])) pip_cali[[s]][[name]] = c(sum(tmp[[name]]), sum(tmp$is_signal), length(tmp$is_signal))
        else pip_cali[[s]][[name]] = rbind(pip_cali[[s]][[name]], c(sum(tmp[[name]]), sum(tmp$is_signal), length(tmp$is_signal)))
      }
      pip_cali[[s]][[name]][which(is.na(pip_cali[[s]][[name]]))] = 0
    }
  }
  saveRDS(list("lm in sample"=get_cali(pip_cali, 'susierss_lm_insample'),
               "mix in sample"=get_cali(pip_cali, 'susierss_mix_insample'),
               "lm ref sample"=get_cali(pip_cali, 'susierss_lm_refsample'), 
               "mix ref sample"=get_cali(pip_cali, 'susierss_mix_refsample')),
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
