library(dplyr)
## Functions
merge_img = function(prefix, n) {
  files = paste0(prefix, '_', seq(1:n), '.png')
  cmd = paste('convert +append', paste(files, collapse=" "), paste0(prefix, '.png'))
  system(cmd)
  files = paste0(prefix, '_', seq(1:n), '.md')
  cmd = paste('cat', paste(files, collapse=" "), '>', paste0(prefix, '.md'))
  system(cmd)
  system(paste('rm -f', paste(files, collapse=" ")))
}

merge_lists = function(lists) {
  lists = do.call(rbind,lapply(cbind(lists), unlist))
  names = colnames(lists)
  lists = lapply(1:ncol(lists), function(i) lists[,i])
  names(lists) = names
  return(lists)
}

plot_pip = function(x, n_causal, s, output_prefix, xname, yname, xlab, ylab, pip_cutoff = -1, pip_diff_categories = c(0.1, 0.15, 0.2)) {
  x = x[x[[xname]] > pip_cutoff & x[[yname]] > pip_cutoff,]
  if (pip_cutoff >=0) {
    xlab = paste0(xlab, ">", pip_cutoff)
    ylab = paste0(ylab, ">", pip_cutoff)
  }
  colors = sapply(1:length(x$is_signal), function(i) ifelse(x$is_signal[i],'#800000','#D3D3D3'))
  # compute difference in pip and proportions
  pip_diff = abs(x[[xname]] - x[[yname]])
  pip_diff_cnts = sapply(pip_diff_categories, function(x) length(which(pip_diff >= x)))
  pip_diff_text = vector()
  for (i in 1:length(pip_diff_categories)) {
    tmp = paste0('- ', pip_diff_cnts[i], '/', length(pip_diff), ' (', formatC(pip_diff_cnts[i] / length(pip_diff) * 100, format = "e", digits = 2), '%) differ by ', pip_diff_categories[i])
    pip_diff_text = c(pip_diff_text, tmp)
  }
  pip_diff_text = paste(pip_diff_text, collapse = '\n')
  corr_text = paste('- correlation', round(cor(x)[1,2],2))
  header = paste('#', xname, 'vs', yname, s, n_causal, 'causal\n')
  pdf(paste0(output_prefix, '_', n_causal, '.pdf'), width=5, height=5, pointsize=10)
  plot(x[[xname]][colors == '#D3D3D3'], x[[yname]][colors == '#D3D3D3'], xlab = xlab, ylab = ylab,
       main = paste0(n_causal, s , ifelse(n_causal>1, ' effect variables', ' effect variable')),
       col = '#D3D3D3', pch = 20, cex = 1.2, bty='l', xlim=c(0,1), ylim=c(0,1))
  points(x[[xname]][colors == '#800000'], x[[yname]][colors == '#800000'], col = '#800000', pch=20, cex=1.2)
  abline(0,1,col=2)
  #abline(h=0.95, col='gray')
  #abline(v=0.95, col='gray')
  dev.off()
  system(paste0("convert -flatten -density 120 ", paste0(output_prefix, '_', n_causal, '.pdf'), " ", paste0(output_prefix, '_', n_causal, '.png')))
  write(paste(header, corr_text, pip_diff_text, '\n', sep = '\n'), paste0(output_prefix, '_', n_causal, '.md'))
}

## parameters
initial = 'lasso'
add_z = c(FALSE, TRUE)
ldmethod = c('in_sample', 'ref_sample')
estimate_resid = c(FALSE, TRUE)
lambda = 0.1
all.comb = expand.grid(initial, add_z, ldmethod, estimate_resid, lambda)
colnames(all.comb) = c('initial','add_z', 'ldmethod', 'estimate_resid', 'lambda')
all.comb = all.comb %>% filter(!(ldmethod == 'in_sample' & add_z == TRUE))

for (case in 1:nrow(all.comb)){
  initial = all.comb[case, 'initial']
  add_z = all.comb[case, 'add_z']
  ldmethod = all.comb[case, 'ldmethod']
  estimate_resid = all.comb[case, 'estimate_resid']
  lambda = all.comb[case, 'lambda']
  input = paste0('susie_rss_ukb_default_pip_extraction/susie_rss_ukb_default_pip_init',initial,'_AZ',add_z,'_ld',ldmethod,'_ER',
                 estimate_resid, '_lamb', lambda,'.rds')
  output = paste0('susie_rss_ukb_default_pip_comparison/susie_rss_ukb_default_pip_init',initial,'_AZ',add_z,'_ld',ldmethod,'_ER',
                  estimate_resid, '_lamb', lambda)
  result = readRDS(input)
  
  result_1 = result[[1]]
  result_2 = result[[2]]
  result_3 = result[[3]]
  result_5 = do.call(rbind, lapply(3:5, function(i) result[[i]][,c('susie','susierss', 'is_signal')]))
  
  # susierss vs susie
  plot_pip(result_1, 1, '', paste0(output, '_rssVSsusie'), 'susierss', 'susie', 'PIP SuSiE-RSS', 'PIP SuSiE')
  plot_pip(result_2, 2, '', paste0(output, '_rssVSsusie'), 'susierss', 'susie', 'PIP SuSiE-RSS', 'PIP SuSiE')
  plot_pip(result_5, 3, ' ~ 5', paste0(output, '_rssVSsusie'), 'susierss', 'susie', 'PIP SuSiE-RSS', 'PIP SuSiE')
  merge_img(paste0(output, '_rssVSsusie'), 3)
  
  # susierss vs caviar
  plot_pip(result_1, 1, '', paste0(output, '_rssVScaviar'), 'susierss', 'caviar', 'PIP SuSiE-RSS', 'PIP CAVIAR')
  plot_pip(result_2, 2, '', paste0(output, '_rssVScaviar'), 'susierss', 'caviar', 'PIP SuSiE-RSS', 'PIP CAVIAR')
  plot_pip(result_3, 3, '', paste0(output, '_rssVScaviar'), 'susierss', 'caviar', 'PIP SuSiE-RSS', 'PIP CAVIAR')
  merge_img(paste0(output, '_rssVScaviar'), 3)
  
  # susierss vs finemap
  plot_pip(result_1, 1, '', paste0(output, '_rssVSfinemap'), 'susierss', 'finemap', 'PIP SuSiE-RSS', 'PIP FINEMAP')
  plot_pip(result_2, 2, '', paste0(output, '_rssVSfinemap'), 'susierss', 'finemap', 'PIP SuSiE-RSS', 'PIP FINEMAP')
  plot_pip(result_3, 3, '', paste0(output, '_rssVSfinemap'), 'susierss', 'finemap', 'PIP SuSiE-RSS', 'PIP FINEMAP')
  merge_img(paste0(output, '_rssVSfinemap'), 3)
  
  # susierss vs finemapv3
  plot_pip(result_1, 1, '', paste0(output, '_rssVSfinemapv3'), 'susierss', 'finemapv3', 'PIP SuSiE-RSS', 'PIP FINEMAPv3')
  plot_pip(result_2, 2, '', paste0(output, '_rssVSfinemapv3'), 'susierss', 'finemapv3', 'PIP SuSiE-RSS', 'PIP FINEMAPv3')
  plot_pip(result_3, 3, '', paste0(output, '_rssVSfinemapv3'), 'susierss', 'finemapv3', 'PIP SuSiE-RSS', 'PIP FINEMAPv3')
  merge_img(paste0(output, '_rssVSfinemapv3'), 3)
  
  # susie vs caviar
  plot_pip(result_1, 1, '', paste0(output, '_susieVScaviar'), 'susie', 'caviar', 'PIP SuSiE', 'PIP CAVIAR')
  plot_pip(result_2, 2, '', paste0(output, '_susieVScaviar'), 'susie', 'caviar', 'PIP SuSiE', 'PIP CAVIAR')
  plot_pip(result_3, 3, '', paste0(output, '_susieVScaviar'), 'susie', 'caviar', 'PIP SuSiE', 'PIP CAVIAR')
  merge_img(paste0(output, '_susieVScaviar'), 3)
  
  # susie vs finemap
  plot_pip(result_1, 1, '', paste0(output, '_susieVSfinemap'), 'susie', 'finemap', 'PIP SuSiE', 'PIP FINEMAP')
  plot_pip(result_2, 2, '', paste0(output, '_susieVSfinemap'), 'susie', 'finemap', 'PIP SuSiE', 'PIP FINEMAP')
  plot_pip(result_3, 3, '', paste0(output, '_susieVSfinemap'), 'susie', 'finemap', 'PIP SuSiE', 'PIP FINEMAP')
  merge_img(paste0(output, '_susieVSfinemap'), 3)
  
  # susie vs finemapv3
  plot_pip(result_1, 1, '', paste0(output, '_susieVSfinemapv3'), 'susie', 'finemapv3', 'PIP SuSiE', 'PIP FINEMAPv3')
  plot_pip(result_2, 2, '', paste0(output, '_susieVSfinemapv3'), 'susie', 'finemapv3', 'PIP SuSiE', 'PIP FINEMAPv3')
  plot_pip(result_3, 3, '', paste0(output, '_susieVSfinemapv3'), 'susie', 'finemapv3', 'PIP SuSiE', 'PIP FINEMAPv3')
  merge_img(paste0(output, '_susieVSfinemapv3'), 3)
  
}
