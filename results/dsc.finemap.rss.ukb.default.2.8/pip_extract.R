input = 'susie_rss_ukb_default_query.rds'

library(dplyr)
dat = readRDS(input)

# parameters
add_z = c(FALSE, TRUE)
ldmethod = c('in_sample', 'ref_sample')
estimate_resid = c(FALSE, TRUE)
lambda = c(0, 1e-04, 0.1)
initial = c(NA, 'oracle', 'lasso')
all.comb = expand.grid(initial, add_z, ldmethod, estimate_resid, lambda)
colnames(all.comb) = c('initial','add_z', 'ldmethod', 'estimate_resid', 'lambda')
all.comb = all.comb %>% filter(!(ldmethod == 'in_sample' & add_z == TRUE))

for(case in 1:nrow(all.comb)){
  initial = all.comb[case, 'initial']
  add_z = all.comb[case, 'add_z']
  ldmethod = all.comb[case, 'ldmethod']
  estimate_resid = all.comb[case, 'estimate_resid']
  lambda = all.comb[case, 'lambda']
  output = paste0('susie_rss_ukb_default_pip_extraction/susie_rss_ukb_default_pip_init',initial,'_AZ',add_z,'_ld',ldmethod,'_ER',
                  estimate_resid, '_lamb', lambda,'.rds')
  
  caviar_out = dat$caviar %>% filter(addz == add_z, ld_method == ldmethod)
  finemap_out = dat$finemap %>% filter(addz == add_z, ld_method == ldmethod)
  finemapv3_out = dat$finemapv3 %>% filter(addz == add_z, ld_method == ldmethod)
  susie_out = dat$susie %>% filter(method_susie == 'susie',
                                   estimate_residual_variance == estimate_resid)
  if(is.na(initial)){
    susie_out = susie_out[is.na(susie_out$init),]
  }else{
    susie_out = susie_out[which(susie_out$init == initial),]
  }
  susierss_out = dat$susie %>% filter(method_susie != 'susie',
                                      addz == add_z,
                                      estimate_residual_variance == estimate_resid,
                                      ld_method == ldmethod, lamb == lambda)
  if(is.na(initial)){
    susierss_out = susierss_out[is.na(susierss_out$init),]
  }else{
    susierss_out = susierss_out[which(susierss_out$init == initial),]
  }
  # and for the rest of the PIP analysis we pool across all data-set
  # but we evaluate for each number of signals settings
  # PVE is fixed to 0.2 here using `lm_pve02` module
  data_sets = unique(susierss_out$dataset)
  n_signals = unique(susierss_out$n_signal)
  result = list()
  for (s in n_signals) {
    result[[as.character(s)]] = NULL
    if (s > 3) {
      has_caviar = FALSE
    } else {
      has_caviar = TRUE
    }
    print(paste('==============', s, '=============='))
    for (d in data_sets) {
      # should all be one row
      # if(s > 3){
      susie_row = susie_out %>% filter(n_signal == s, dataset == d)
      susierss_row = susierss_out %>% filter(n_signal == s, dataset == d)
      # }else{
      #   susiebhat_row = susiebhat_out %>% filter(n_signal == s, dataset == d, L == s)
      #   susierss_row = susierss_out %>% filter(n_signal == s, dataset == d, L == s)
      # }
      if (nrow(susie_row) != 1 || nrow(susierss_row) != 1) stop("Rows are not unique for susie & susierss")
      converged = as.logical(susie_row$converged[[1]] * susierss_row$converged[[1]])
      if(any(!is.na(converged)) & any(converged)){
        if (has_caviar) {
          caviar_row = caviar_out %>% filter(n_signal == s, dataset == d, args == paste('-g 0.001 -c', s))
          finemap_row = finemap_out %>% filter(n_signal == s, dataset == d, args == paste('--n-causal-max', s))
          finemapv3_row = finemapv3_out %>% filter(n_signal == s, dataset == d, args == paste('--n-causal-snps', s))
          if (nrow(caviar_row) != 1 || nrow(finemap_row) != 1 || nrow(finemapv3_row) != 1) stop("Rows are not unique for caviar & finemap & finemapv3")
        }
        # a slightly awkward yet convenient syntax thanks to tibble format
        truth = susie_row$meta[[1]]$true_coef
        for (r in which(converged == TRUE)) {
          susie_pip = susie_row$pip[[1]][,r]
          susierss_pip = susierss_row$pip[[1]][,r]
          if (has_caviar) {
            caviar_pip = caviar_row$pip[[1]][,r]
            finemap_pip = finemap_row$pip[[1]][,r]
            finemapv3_pip = finemapv3_row$pip[[1]][,r]
            pip = cbind(susie_pip, susierss_pip, caviar_pip, finemap_pip, finemapv3_pip, truth[,r] != 0)
          } else {
            pip = cbind(susie_pip, susierss_pip, truth[,r] != 0)
          }
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
    if (has_caviar) {
      colnames(result[[as.character(s)]]) = c('susie', 'susierss', 'caviar', 'finemap', 'finemapv3', 'is_signal')
    } else {
      colnames(result[[as.character(s)]]) = c('susie', 'susierss', 'is_signal')
    }
  }
  saveRDS(result, output)
}
