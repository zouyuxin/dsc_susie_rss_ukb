output_dir = '../../output/rss_compare_ukb_default20201221/'
output = 'susie_rss_ukb_default_query.rds'

library(tibble)
susie_out = dscrutils::dscquery(output_dir,
                                targets = c("small_data.dataset","small_data.r_ref_2dist", "small_data.r_ref_Mdist",
                                            "sim_gaussian", "sim_gaussian.meta", "sim_gaussian.n_signal",
                                            "method_susie_full","method_susie_full.lamb",
                                            "method_susie_full.ld_method", "method_susie_full.addz",
                                            "method_susie_full.z_ld_weight",
                                            "method_susie_full.estimate_residual_variance",
                                            "method_susie_full.DSC_TIME",
                                            "score_susie.total", "score_susie.valid",
                                            "score_susie.size", "score_susie.purity", 
                                            "score_susie.avgr2", "score_susie.niter",
                                            "score_susie.top", "score_susie.converged","score_susie.objective",
                                            "score_susie.overlap", "score_susie.signal_pip", "score_susie.pip"),
                                    module.output.files = "sim_gaussian",
                                groups = c('method_susie_ldweight:'))

caviar_out = dscrutils::dscquery(output_dir,
                                 targets = c("small_data.dataset",
                                             "sim_gaussian","sim_gaussian.meta", "sim_gaussian.n_signal",
                                             "caviar.ld_method", "caviar.addz",
                                             "caviar.args", "caviar.DSC_TIME",
                                             "score_caviar.total", "score_caviar.valid", "score_caviar.size",
                                             "score_caviar.signal_pip", "score_caviar.pip"),
                                 module.output.files = "sim_gaussian",
                                 groups = c('method_susie_ldweight:'))

finemap_out = dscrutils::dscquery(output_dir,
                                  targets = c("small_data.dataset",
                                              "sim_gaussian", "sim_gaussian.meta", "sim_gaussian.n_signal",
                                              "finemap.args", "finemap.ld_method",
                                              "finemap.addz", "finemap.DSC_TIME",
                                              "score_finemap.total", "score_finemap.valid", "score_finemap.size",
                                              "score_finemap.signal_pip", "score_finemap.pip"),
                                  module.output.files ="sim_gaussian",
                                  groups = c('method_susie_ldweight:'))

rename_cols = function(dat) {
  for (item in names(dat)) {
    tmp = strsplit(colnames(dat[[item]]), "[.]")
    colnames(dat[[item]]) = unlist(lapply(1:length(tmp), function(i) ifelse(length(tmp[[i]])>1, tmp[[i]][2], tmp[[i]][1])))
  }
  return(dat)
}
# remove module names from column names; this is Okay because every column field here are unique
res = rename_cols(list(susie=as_tibble(susie_out),caviar=as_tibble(caviar_out), 
                       finemap=as_tibble(finemap_out)))
saveRDS(res, output)
