# workhorse(s) for finemapping

# Module input
# ============
# $X, $Y: full data; or
# $sumstats: summary statistics; or / and
# $ld: LD information

# Module output
# =============
# $fitted: for diagnostics
# $posterior: for inference

caviar: fit_caviar.R + add_z.R + R(posterior = finemap_mcaviar(z,ld_file, args, prefix=cache))
  @CONF: R_libs = (dplyr, magrittr, data.table)
  sumstats: $sumstats
  ld: $ld
  N_ref: $N_ref
  (addz, ld_method): (FALSE, "in_sample"),(FALSE, "ref_sample"),(TRUE, "ref_sample")
  ld_ref_z_file: file(ref.z.ld)
  args: "-g 0.001 -c 1", "-g 0.001 -c 2", "-g 0.001 -c 3"
  cache: file(CAVIAR)
  $posterior: posterior

caviar_simple(caviar):
  args: "-g 0.001 -c 3"

finemap(caviar): fit_finemap.R + add_z.R + R(posterior = finemap_mvar(z,ld_file, N_in, k, args, prefix=cache))
  k: NULL
  N_in: $N_sample
  args: "--n-causal-max 1", "--n-causal-max 2", "--n-causal-max 3"
  cache: file(FM)

finemap_simple(finemap):
  args: "--n-causal-max 3"

finemapv3(caviar): fit_finemap_v3.R + add_z.R + R(posterior = finemap_mvar_v1.3.1(sumstats$bhat, sumstats$shat,
                                        maf[[ld_method]], ld_file, N_in, k, method, args, prefix=cache))
  k: NULL
  maf: $maf
  N_in: $N_sample
  method: 'sss'
  args: "--n-causal-snps 1", "--n-causal-snps 2", "--n-causal-snps 3"
  cache: file(FM)

finemapv3_simple(finemapv3):
  args: "--n-causal-snps 3"

dap_z: fit_dap.py + Python(z = sumstats['bhat']/sumstats['shat'];
                           numpy.nan_to_num(z, copy=False);
                           posterior = dap_batch_z(z, ld[ld_method], cache, args))
  sumstats: $sumstats
  ld: $ld
  ld_method: "in_sample", "ref_sample"
  args: "-ld_control 0.20 --all"
  cache: file(DAP)
  $posterior: posterior

susie: initialize.R + fit_susie.R
  # Prior variance of nonzero effects.
  @CONF: R_libs = susieR
  maxI: 1000
  maxL: 10
  null_weight: 0
  prior_var: 0
  X: $X_sample
  Zpve: $Zpve
  Z: $PC_sample
  Y: $Y
  meta: $meta
  estimate_residual_variance: TRUE
  init: NA
  $posterior: posterior
  $fitted: fitted

susie_init(susie):
  init: NA, 'oracle', 'lasso'

susie_auto: fit_susie.R
  @CONF: R_libs = susieR
  X: $X_sample
  Y: $Y
  prior_var: "auto"
  $posterior: posterior
  $fitted: fitted

susie01(susie):
  null_weight: 0
  maxL: 1

susie10(susie):
  null_weight: 0
  maxL: 15
  prior_var: 0, 0.1

#------------------------------
# SuSiE with summary statistics
#------------------------------

susie_rss: initialize.R + R(if(is.na(init)){
                          s_init = NULL
                        }else if(init == 'oracle'){
                          s_init = init_susie_rss_true($(meta)$true_coef, n)
                        }else if(init == 'lasso'){
                          s_init = init_rss_lasso(z,r,L)
                        }) + fit_susie_rss.R
  @CONF: R_libs = (susieR, data.table)
  sumstats: $sumstats
  ld: $ld
  L: 10
  z_ld_weight: 0, 0.002
  n: $N_sample
  estimate_residual_variance: TRUE
  ld_method: "in_sample", "ref_sample"
  init: NA
  $fitted: res$fitted
  $posterior: res$posterior

susie_rss_large(susie_rss):
  L: 15

susie_rss_simple(susie_rss):
  L: 10
  
susie_rss_zldweight(susie_rss):
  z_ld_weight: 0, 0.001, 0.002, 0.005, 0.01

susie_rss_init(susie_rss):
  init: NA, 'oracle', 'lasso'

susie_rss_lambda: add_z_susierss.R + initialize.R + R(if(is.na(init)){
                          s_init = NULL
                        }else if(init == 'oracle'){
                          s_init = init_susie_rss_true($(meta)$true_coef, n)
                        }else if(init == 'lasso'){
                          s_init = init_rss_lasso(z,r,L)
                        }) + fit_susie_rss_lambda.R
  @CONF: R_libs = (susieR, data.table)
  sumstats: $sumstats
  ld: $ld
  L: 10
  n: $N_sample
  N_ref: $N_ref
  estimate_residual_variance: TRUE
  lamb: 0.0001, 0.01
  addz: FALSE, TRUE
  ld_method: "in_sample", "ref_sample"
  init: NA
  $fitted: res$fitted
  $posterior: res$posterior

susie_rss_lambda_init(susie_rss_lambda):
  init: NA, 'oracle', 'lasso'
