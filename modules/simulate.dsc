sim_gaussian: simulate.R + \
                R(neale_effects = readRDS(effects);
                set.seed(seed);
                res=sim_gaussian_multiple(X, Z, pve, Z_pve, n_signal, neale_effects$beta, n_traits, file_name, sample_file))
  @CONF: R_libs = susieR
  seed: $seed
  X: $X_sample
  Z: $Z_sample
  pve: 0.01
  Z_pve: 0.01
  n_signal: 3
  effects: ${neale_effect_size}
  n_traits: 1
  file_name: file(pheno)
  sample_file: $sample_file
  $Y: res$Y
  $meta: res$meta
  $pheno_file: file_name

sim_gaussian_null(sim_gaussian):
  pve: 0
  n_signal: 0

sim_gaussian_large(sim_gaussian):
  n_signal: 200
