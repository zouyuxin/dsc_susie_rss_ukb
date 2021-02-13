output_dir = '../../output/rss_compare_ukb_mix'
output = 'susie_rss_ukb_mix_query.rds'

library(tibble)
library(dplyr)
susie_out = dscrutils::dscquery(output_dir,
                                targets = c("small_data.dataset",
                                            "sim_gaussian_simple.meta",
                                            "get_sumstats_simple.method",
                                            'susie_rss_mix.z_ld_weight',
                                            "susie_rss_mix.estimate_residual_variance",
                                            "susie_rss_mix.DSC_TIME", "score_susie",
                                            "score_susie.total", "score_susie.valid",
                                            "score_susie.size", "score_susie.purity", 
                                            "score_susie.avgr2", "score_susie.niter",
                                            "score_susie.top", "score_susie.converged","score_susie.objective",
                                            "score_susie.overlap", "score_susie.signal_pip", "score_susie.pip"),
                                    module.output.files = "score_susie", ignore.missing.files = TRUE)


rename_cols = function(dat) {
  for (item in names(dat)) {
    tmp = strsplit(colnames(dat[[item]]), "[.]")
    colnames(dat[[item]]) = unlist(lapply(1:length(tmp), function(i) ifelse(length(tmp[[i]])>1, tmp[[i]][2], tmp[[i]][1])))
  }
  return(dat)
}
# remove module names from column names; this is Okay because every column field here are unique
res = rename_cols(list(susie=as_tibble(susie_out)))
susierss_lm_in_out = res$susie %>% filter(method == 'lm',z_ld_weight==0)
susierss_lm_inz_out = res$susie %>% filter(method == 'lm',z_ld_weight==0.0001)
susierss_mix_in_out = res$susie %>% filter(method == 'mixed', z_ld_weight==0)
susierss_mix_inz_out = res$susie %>% filter(method == 'mixed',z_ld_weight==0.0001)

saveRDS(list(susierss_lm_insample = susierss_lm_in_out,
             susierss_lm_insamplez = susierss_lm_inz_out,
             susierss_mix_insamplez = susierss_mix_inz_out,
             susierss_mix_insample = susierss_mix_in_out), output)

