output_dir = '../../output/rss_compare_ukb'
output = 'susie_rss_ukb_default_query.rds'

library(tibble)
susie_out = dscrutils::dscquery(output_dir,
                                targets = c("small_data.dataset","small_data.r_Fdist", "small_data.r_Mdist",
                                            "sim_gaussian", "sim_gaussian.meta", "sim_gaussian.n_signal",
                                            "method_susie","method_susie.init","method_susie.lamb",
                                            "method_susie.ld_method", "method_susie.addz",
                                            "method_susie.estimate_residual_variance",
                                            "method_susie.DSC_TIME",
                                            "score_susie.total", "score_susie.valid",
                                            "score_susie.size", "score_susie.purity", 
                                            "score_susie.avgr2", "score_susie.niter",
                                            "score_susie.top", "score_susie.converged","score_susie.objective",
                                            "score_susie.overlap", "score_susie.signal_pip", "score_susie.pip"),
                                    module.output.files = "sim_gaussian")

caviar_out = dscrutils::dscquery(output_dir,
                                 targets = c("small_data.dataset", "small_data.r_Fdist", "small_data.r_Mdist",
                                             "sim_gaussian","sim_gaussian.meta", "sim_gaussian.n_signal",
                                             "caviar.ld_method", "caviar.addz",
                                             "caviar.args", "caviar.DSC_TIME",
                                             "score_caviar.total", "score_caviar.valid", "score_caviar.size",
                                             "score_caviar.signal_pip", "score_caviar.pip"),
                                 module.output.files = "sim_gaussian")

finemap_out = dscrutils::dscquery(output_dir,
                                  targets = c("small_data.dataset","small_data.r_Fdist", "small_data.r_Mdist",
                                              "sim_gaussian", "sim_gaussian.meta", "sim_gaussian.n_signal",
                                              "finemap.args", "finemap.ld_method",
                                              "finemap.addz", "finemap.DSC_TIME",
                                              "score_finemap.total", "score_finemap.valid", "score_finemap.size",
                                              "score_finemap.signal_pip", "score_finemap.pip"),
                                  module.output.files ="sim_gaussian")

finemapv3_out = dscrutils::dscquery(output_dir,
                                  targets = c("small_data.dataset","small_data.r_Fdist", "small_data.r_Mdist",
                                              "sim_gaussian", "sim_gaussian.meta", "sim_gaussian.n_signal",
                                              "finemapv3.args", "finemapv3.ld_method",
                                              "finemapv3.addz", "finemapv3.DSC_TIME",
                                              "score_finemapv3.total", "score_finemapv3.valid", 
                                              "score_finemapv3.size",
                                              "score_finemapv3.signal_pip", "score_finemapv3.pip"),
                                  module.output.files ="sim_gaussian")

rename_cols = function(dat) {
  for (item in names(dat)) {
    tmp = strsplit(colnames(dat[[item]]), "[.]")
    colnames(dat[[item]]) = unlist(lapply(1:length(tmp), function(i) ifelse(length(tmp[[i]])>1, tmp[[i]][2], tmp[[i]][1])))
  }
  return(dat)
}
# remove module names from column names; this is Okay because every column field here are unique
res = rename_cols(list(susie=as_tibble(susie_out),caviar=as_tibble(caviar_out), 
                       finemap=as_tibble(finemap_out), finemapv3=as_tibble(finemapv3_out)))
saveRDS(res, output)
