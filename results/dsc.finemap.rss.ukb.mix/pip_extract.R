input = 'susie_rss_ukb_mix_query.rds'

library(dplyr)
dat = readRDS(input)

# parameters
estimate_resid = c(FALSE, TRUE)
all.comb = expand.grid(estimate_resid)
colnames(all.comb) = c('estimate_resid')

for(case in 1:nrow(all.comb)){
  estimate_resid = all.comb[case, 'estimate_resid']
  output = paste0('susie_rss_ukb_mix_pip_extraction/susie_rss_ukb_mix_pip_ER',
                  estimate_resid,'.rds')
  
  susierss_lm_insample = dat$susierss_lm_insample %>% filter(estimate_residual_variance == estimate_resid)
  susierss_mix_insample = dat$susierss_mix_insample %>% filter(estimate_residual_variance == estimate_resid)
  susierss_lm_insamplez = dat$susierss_lm_insamplez %>% filter(estimate_residual_variance == estimate_resid)
  susierss_mix_insamplez = dat$susierss_mix_insamplez %>% filter(estimate_residual_variance == estimate_resid)

  data_sets = unique(susierss_mix_insample$dataset)
  result = list()
  for(d in data_sets){
    susierss_lm_insample_row = susierss_lm_insample %>% filter(dataset == d)
    susierss_mix_insample_row = susierss_mix_insample %>% filter(dataset == d)
    susierss_lm_insamplez_row = susierss_lm_insamplez %>% filter(dataset == d)
    susierss_mix_insamplez_row = susierss_mix_insamplez %>% filter(dataset == d)
    
    if (nrow(susierss_lm_insample_row) != 1 || nrow(susierss_mix_insample_row) != 1 || 
        nrow(susierss_lm_insamplez_row) != 1 || nrow(susierss_mix_insamplez_row) != 1) {
      stop("Rows are not unique for susie & susierss")
    }
    
    converged = as.logical(susierss_lm_insample_row$converged[[1]] * susierss_mix_insample_row$converged[[1]] * 
                             susierss_lm_insamplez_row$converged[[1]] * susierss_mix_insamplez_row$converged[[1]])
    if(any(!is.na(converged)) & any(converged)){
      truth = susierss_mix_insample_row$meta[[1]]$true_coef
      for (r in which(converged == TRUE)) {
        susierss_lm_insample_pip = susierss_lm_insample_row$pip[[1]][,r]
        susierss_mix_insample_pip = susierss_mix_insample_row$pip[[1]][,r]
        susierss_lm_insamplez_pip = susierss_lm_insamplez_row$pip[[1]][,r]
        susierss_mix_insamplez_pip = susierss_mix_insamplez_row$pip[[1]][,r]
        
        pip = cbind(susierss_lm_insample_pip, susierss_mix_insample_pip, 
                    susierss_lm_insamplez_pip, susierss_mix_insamplez_pip, truth[,r] != 0)
        # remove all zero PIP / table
        pip = pip[rowSums(pip) > 0, ]
        if (is.null(result[['2']])) {
          result[['2']] = pip
        } else {
          result[['2']] = rbind(result[['2']], pip)
        }
      }
    }
  }
  result[['2']] = data.frame(result[['2']])
  colnames(result[['2']]) = c('susierss_lm_insample', 'susierss_mix_insample', 
                              'susierss_lm_insamplez', 'susierss_mix_insamplez', 'is_signal')
  saveRDS(result, output)
}
