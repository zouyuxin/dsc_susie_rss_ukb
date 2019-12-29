sim_gaussian: simulate.R + \
                R(neale_effects = readRDS(effects);
                res=sim_gaussian_multiple(X, pve, n_signal, neale_effects$beta, n_traits))
  @CONF: R_libs = susieR
  X: $X_sample
  pve: 0.02
  n_signal: 1,2,3,4,5
  effects: ${neale_effect_size}
  n_traits: 1
  $Y: res$Y
  $meta: res$meta

sim_gaussian_null(sim_gaussian):
  pve: 0
  n_signal: 0

sim_gaussian_large(sim_gaussian):
  n_signal: 200
