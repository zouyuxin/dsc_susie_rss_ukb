input = 'susie_rss_ukb_mix_query.rds'

library(dplyr)
dat = readRDS(input)

# parameters
add_z = c(FALSE, TRUE)
estimate_resid = c(FALSE, TRUE)
all.comb = expand.grid(add_z, estimate_resid)
colnames(all.comb) = c('add_z', 'estimate_resid')

for(case in 1:nrow(all.comb)){
  add_z = all.comb[case, 'add_z']
  estimate_resid = all.comb[case, 'estimate_resid']
  output = paste0('susie_rss_ukb_mix_pip_extraction/susie_rss_ukb_mix_pip_AZ',add_z,'_ER',
                  estimate_resid,'.rds')
  
  susierss_lm_insample = dat$susierss_lm_insample %>% filter(addz == add_z, estimate_residual_variance == estimate_resid)
  susierss_mix_insample = dat$susierss_mix_insample %>% filter(addz == add_z, estimate_residual_variance == estimate_resid)
  susierss_lm_refsample = dat$susierss_lm_refsample %>% filter(addz == add_z, estimate_residual_variance == estimate_resid)
  susierss_mix_refsample = dat$susierss_mix_refsample %>% filter(addz == add_z, estimate_residual_variance == estimate_resid)

  data_sets = unique(susierss_mix_insample$dataset)
  n_signals = unique(susierss_mix_insample$n_signal)
  result = list()
  for (s in n_signals) {
    result[[as.character(s)]] = NULL
    print(paste('==============', s, '=============='))
    for (d in data_sets) {
      # should all be one row
      # if(s > 3){
      susierss_lm_insample_row = susierss_lm_insample %>% filter(n_signal == s, dataset == d)
      susierss_mix_insample_row = susierss_mix_insample %>% filter(n_signal == s, dataset == d)
      susierss_lm_refsample_row = susierss_lm_refsample %>% filter(n_signal == s, dataset == d)
      susierss_mix_refsample_row = susierss_mix_refsample %>% filter(n_signal == s, dataset == d)

      # }else{
      #   susiebhat_row = susiebhat_out %>% filter(n_signal == s, dataset == d, L == s)
      #   susierss_row = susierss_out %>% filter(n_signal == s, dataset == d, L == s)
      # }
      if (nrow(susierss_lm_insample_row) != 1 || nrow(susierss_mix_insample_row) != 1 || 
          nrow(susierss_lm_refsample_row) != 1 || nrow(susierss_mix_refsample_row) != 1) {
        stop("Rows are not unique for susie & susierss")
      }
      converged = as.logical(susierss_lm_insample_row$converged[[1]] * susierss_mix_insample_row$converged[[1]] * 
                               susierss_lm_refsample_row$converged[[1]] * susierss_mix_refsample_row$converged[[1]])
      # a slightly awkward yet convenient syntax thanks to tibble format
      if(any(!is.na(converged)) & any(converged)){
        truth = susierss_mix_insample_row$meta[[1]]$true_coef
        for (r in which(converged == TRUE)) {
          susierss_lm_insample_pip = susierss_lm_insample_row$pip[[1]][,r]
          susierss_mix_insample_pip = susierss_mix_insample_row$pip[[1]][,r]
          susierss_lm_refsample_pip = susierss_lm_refsample_row$pip[[1]][,r]
          susierss_mix_refsample_pip = susierss_mix_refsample_row$pip[[1]][,r]

          pip = cbind(susierss_lm_insample_pip, susierss_mix_insample_pip, 
                      susierss_lm_refsample_pip, susierss_mix_refsample_pip, truth[,r] != 0)
          # remove all zero PIP / table
          pip = pip[rowSums(pip) > 0, ]
          if (is.null(result[[as.character(s)]])) {
             result[[as.character(s)]] = pip
          } else {
             result[[as.character(s)]] = rbind(result[[as.character(s)]], pip)
          }
        }
      }
    }
    result[[as.character(s)]] = data.frame(result[[as.character(s)]])
    colnames(result[[as.character(s)]]) = c('susierss_lm_insample', 'susierss_mix_insample', 
                                            'susierss_lm_refsample', 'susierss_mix_refsample', 'is_signal')
  }
  saveRDS(result, output)
}
