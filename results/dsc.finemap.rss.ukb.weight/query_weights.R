out = readRDS('susie_rss_ukb_weight_query.rds')

library(dplyr)

susie_out = out$susie %>% filter(method_susie == 'susie')
susierss_z0 = out$susie %>% filter(method_susie == 'susie_rss', 
                                   z_ld_weight == 0)
susierss_z001 = out$susie %>% filter(method_susie == 'susie_rss', 
                                     z_ld_weight == 0.001)
susierss_z002 = out$susie %>% filter(method_susie == 'susie_rss', 
                                     z_ld_weight == 0.002)
susierss_z005 = out$susie %>% filter(method_susie == 'susie_rss', 
                                     z_ld_weight == 0.005)
susierss_z01 = out$susie %>% filter(method_susie == 'susie_rss', 
                                     z_ld_weight == 0.01)
susierss_z02 = out$susie %>% filter(method_susie == 'susie_rss', 
                                     z_ld_weight == 0.02)

saveRDS(list(susie = susie_out, susierss_z0 = susierss_z0,
             susierss_z001 = susierss_z001, susierss_z002 = susierss_z002, susierss_z005 = susierss_z005,
             susierss_z01 = susierss_z01, susierss_z02 = susierss_z02), 
        'susie_rss_ukb_weight_query_weights.rds')

