sim_gaussian: simulate.R + \
                R(neale_effects = readRDS(effects);
                set.seed(seed);
                res=sim_gaussian_multiple(X, Z, pve, Z_pve, n_signal, neale_effects$beta, n_traits, file_name, sample_file))
  @CONF: R_libs = susieR
  seed: $seed
  X: $X_sample
  Z: $PC_sample
  pve: 0.02
  Z_pve: ${Z_pve}
  n_signal: 1,2,3,4,5
  effects: ${neale_effect_size}
  n_traits: 1
  file_name: file(pheno)
  sample_file: $sample_file
  $Y: res$Y
  $meta: res$meta
  $pheno_file: file_name

sim_gaussian_s2(sim_gaussian):
  X: $X_sample_batch
  n_signal: 2

sim_gaussian_null(sim_gaussian):
  pve: 0
  n_signal: 0

sim_gaussian_large(sim_gaussian):
  n_signal: 200
