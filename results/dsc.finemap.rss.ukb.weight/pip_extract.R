input = 'susie_rss_ukb_weight_query_weights.rds'

library(dplyr)
dat = readRDS(input)

# parameters
ldmethod = c('in_sample_Z', 'ref_sample_Z')
all.comb = expand.grid(ldmethod)

colnames(all.comb) = c('ldmethod')

for(case in 1:nrow(all.comb)){
  ldmethod = all.comb[case, 'ldmethod']
  output = paste0('susie_rss_ukb_weight_pip_extraction/susie_rss_ukb_weight_pip_ld',ldmethod,
                  '.rds')
  susie_out = dat$susie
  susierss_z0_out = dat$susierss_z0 %>% filter(ld_method == ldmethod)
  susierss_z001_out = dat$susierss_z001 %>% filter(ld_method == ldmethod)
  susierss_z002_out = dat$susierss_z002 %>% filter(ld_method == ldmethod)
  susierss_z005_out = dat$susierss_z005 %>% filter(ld_method == ldmethod)
  susierss_z01_out = dat$susierss_z01 %>% filter(ld_method == ldmethod)
  # susierss_z02_out = dat$susierss_z02 %>% filter(ld_method == ldmethod)
  # susierss_z05_out = dat$susierss_z05 %>% filter(ld_method == ldmethod)
  
  data_sets = unique(susie_out$dataset)
  result = NULL
  for (d in data_sets) {
    # should all be one row
    susie_row = susie_out %>% filter(dataset == d)
    susierss_z0_row = susierss_z0_out %>% filter(dataset == d)
    susierss_z001_row = susierss_z001_out %>% filter(dataset == d)
    susierss_z002_row = susierss_z002_out %>% filter(dataset == d)
    susierss_z005_row = susierss_z005_out %>% filter(dataset == d)
    susierss_z01_row = susierss_z01_out %>% filter(dataset == d)
    # susierss_z02_row = susierss_z02_out %>% filter(dataset == d)
    # susierss_z05_row = susierss_z05_out %>% filter(dataset == d)

    if (nrow(susie_row) != 1 || nrow(susierss_z0_row) != 1 || nrow(susierss_z001_row) != 1 || 
        nrow(susierss_z002_row) != 1 || nrow(susierss_z005_row) != 1 ||
        nrow(susierss_z01_row) != 1 ) #|| nrow(susierss_z02_row) != 1 || nrow(susierss_z05_row) != 1) 
      stop("Rows are not unique for susie & susierss")
    converged = as.logical(susie_row$converged[[1]] * susierss_z0_row$converged[[1]] * 
                             susierss_z001_row$converged[[1]] * susierss_z002_row$converged[[1]] * 
                             susierss_z005_row$converged[[1]] * susierss_z01_row$converged[[1]]) 
    # * susierss_z02_row$converged[[1]] * susierss_z05_row$converged[[1]])
    if(any(!is.na(converged)) & any(converged)){
      # a slightly awkward yet convenient syntax thanks to tibble format
      truth = susie_row$meta[[1]]$true_coef
      for (r in which(converged == TRUE)) {
        susie_pip = susie_row$pip[[1]][,r]
        susierss_z0_pip = susierss_z0_row$pip[[1]][,r]
        susierss_z001_pip = susierss_z001_row$pip[[1]][,r]
        susierss_z002_pip = susierss_z002_row$pip[[1]][,r]
        susierss_z005_pip = susierss_z005_row$pip[[1]][,r]
        susierss_z01_pip = susierss_z01_row$pip[[1]][,r]
        # susierss_z02_pip = susierss_z02_row$pip[[1]][,r]
        # susierss_z05_pip = susierss_z05_row$pip[[1]][,r]
          
        pip = cbind(susie_pip, susierss_z0_pip, susierss_z001_pip, 
                    susierss_z002_pip, susierss_z005_pip, susierss_z01_pip, truth[,r] != 0)
                    # susierss_z02_pip, susierss_z05_pip, truth[,r] != 0)
          
        # remove all zero PIP / table
        pip = pip[rowSums(pip) > 0, ]
        if (is.null(result)) {
          result = pip
        } else {
          result = rbind(result, pip)
        }
      }
    }
  }
  result = data.frame(result)
  colnames(result) = c('susie', 'RSS_z0', 'RSS_z001', 'RSS_z002', 'RSS_z005', 
                       'RSS_z01', 'is_signal') # 'RSS_z02', 'RSS_z05', 
  saveRDS(result, output)
}
  
