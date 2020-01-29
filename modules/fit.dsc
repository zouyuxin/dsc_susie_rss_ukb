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

add_z: add_z.R
  sumstats: $sumstats
  ld: $ld
  N_ref: $N_ref
  (addz, ld_method): (FALSE, "in_sample"), (FALSE, "ref_sample"), (TRUE, "ref_sample")
  ld_ref_z_file: file(ref.z.ld)
  $z: z
  $ld_file: ld_file
  $ldmethod: ld_method

caviar: fit_caviar.R + R(posterior = finemap_mcaviar(z,ld_file, args, prefix=cache))
  @CONF: R_libs = (dplyr, magrittr, data.table)
  z: $z
  ld_file: $ld_file
  args: "-g 0.001 -c 1", "-g 0.001 -c 2", "-g 0.001 -c 3"
  cache: file(CAVIAR)
  $posterior: posterior

caviar_simple(caviar):
  args: "-g 0.001 -c 3"

finemap(caviar): fit_finemap.R + R(posterior = finemap_mvar(z,ld_file, N_in, k, args, prefix=cache))
  k: NULL
  N_in: $N_sample
  args: "--n-causal-max 1", "--n-causal-max 2", "--n-causal-max 3"
  cache: file(FM)

finemap_simple(finemap):
  args: "--n-causal-max 3"

finemapv3(caviar): fit_finemap_v3.R + R(posterior = finemap_mvar_v1.3.1(sumstats$bhat, sumstats$shat,
                                        maf[[$(ldmethod)]], ld_file, N_in, k, method, args, prefix=cache))
  k: NULL
  sumstats: $sumstats
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

susie: fit_susie.R
  # Prior variance of nonzero effects.
  @CONF: R_libs = susieR
  maxI: 200
  maxL: 10
  null_weight: 0
  prior_var: 0
  X: $X_sample
  Y: $Y
  estimate_residual_variance: TRUE, FALSE
  $posterior: posterior
  $fitted: fitted

susie_auto: fit_susie.R
  @CONF: R_libs = susieR
  X: $X_sample
  Y: $Y
  prior_var: "auto"
  $posterior: posterior
  $fitted: fitted

init_lasso: initialize.R + R(s_init=init_lasso(X,Y,L))
  @CONF: R_libs = susieR
  X: $X_sample
  Y: $Y
  L: 10
  $s_init: s_init

init_oracle: initialize.R + R(s_init=init_susie_true($(meta)$true_coef))
  @CONF: R_libs = susieR
  $s_init: s_init

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

add_z_susierss: add_z_susierss.R
  sumstats: $sumstats
  ld: $ld
  (addz, ld_method): (FALSE, "in_sample"), (FALSE, "ref_sample"), (TRUE, "ref_sample")
  N_ref: $N_ref
  $z: z
  $r: r


susie_rss: fit_susie_rss.R
  @CONF: R_libs = (susieR, data.table)
  z: $z
  r: $r
  s_init: NA
  L: 10
  lamb: 0, 1e-4, 0.1
  estimate_residual_variance: TRUE, FALSE
  $fitted: res$fitted
  $posterior: res$posterior

susie_rss_large(susie_rss):
  L: 15

susie_rss_simple(susie_rss):
  L: 10
  lamb: 0

init_susie_rss_oracle: initialize.R + R(s_init = init_susie_rss_true($(meta)$true_coef, n))
  @CONF: R_libs = susieR
  n: $N_sample
  $s_init: s_init

init_susie_rss_lasso: initialize.R + R(s_init=init_rss_lasso(z,r,L))
  @CONF: R_libs = (susieR, glmnet)
  z: $z
  r: $r
  L: 10
  $s_init: s_init

