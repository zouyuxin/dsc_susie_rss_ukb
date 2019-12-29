sim_gaussian: simulate.R + \
                R(print(effects);
                t_neale = readRDS(effects);
                effects = t_neale$tstat / sqrt(t_neale$n_sample);
                res=sim_gaussian_multiple(X, pve, n_signal, effects, n_traits))
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
