dat = readRDS('susie_rss_ukb_default_query.rds')
initial = NA
add_z = FALSE
lambda = 0
estimate_resid = TRUE
ldmethod = 'in_sample'

s = 3

caviar_out = dat$caviar %>% filter(addz == add_z, ld_method == ldmethod, 
                                   n_signal == s, args == paste('-g 0.001 -c', s))
finemap_out = dat$finemap %>% filter(addz == add_z, ld_method == ldmethod, 
                                     n_signal == s, args == paste('--n-causal-max', s))
finemapv3_out = dat$finemapv3 %>% filter(addz == add_z, ld_method == ldmethod, 
                                         n_signal == s, args == paste('--n-causal-snps', s))
susie_out = dat$susie %>% filter(method_susie == 'susie', estimate_residual_variance == estimate_resid, n_signal == s)
susie_out = susie_out[is.na(susie_out$init),]
susierss_out = dat$susie %>% filter(method_susie != 'susie', addz == add_z,
                                    estimate_residual_variance == estimate_resid,
                                    ld_method == ldmethod, lamb == lambda,  n_signal == s)
susierss_out = susierss_out[is.na(susierss_out$init),]
summary(susie_out$DSC_TIME)
summary(susierss_out$DSC_TIME)
summary(caviar_out$DSC_TIME)
summary(finemap_out$DSC_TIME)
summary(finemapv3_out$DSC_TIME)




