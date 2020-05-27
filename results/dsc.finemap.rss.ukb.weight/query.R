output_dir = '../../output/rss_compare_ukb_weight'
output = 'susie_rss_ukb_weight_query.rds'

library(tibble)
out = dscrutils::dscquery(output_dir,
                                targets = c("small_data.dataset","small_data.r_ref_2dist", "small_data.r_ref_Mdist",
                                            "sim_gaussian_s2", "sim_gaussian_s2.meta", "sim_gaussian_s2.n_signal",
                                            "method_susie",
                                            "method_susie.ld_method", "method_susie.z_ld_weight",
                                            "method_susie.DSC_TIME",
                                            "score_susie.total", "score_susie.valid",
                                            "score_susie.size", "score_susie.purity", 
                                            "score_susie.avgr2", "score_susie.niter",
                                            "score_susie.top", "score_susie.converged","score_susie.objective",
                                            "score_susie.overlap", "score_susie.signal_pip", "score_susie.pip"),
                                    module.output.files = "sim_gaussian_s2")

rename_cols = function(dat) {
  for (item in names(dat)) {
    tmp = strsplit(colnames(dat[[item]]), "[.]")
    colnames(dat[[item]]) = unlist(lapply(1:length(tmp), function(i) ifelse(length(tmp[[i]])>1, tmp[[i]][2], tmp[[i]][1])))
  }
  return(dat)
}

# remove module names from column names; this is Okay because every column field here are unique
res = rename_cols(list(susie=as_tibble(out)))
saveRDS(res, output)
