sim_gaussian: simulate.R + \
                R(set.seed(seed);
                res=sim_gaussian_multiple(X, Z, pve, Z_pve, n_signal, n_traits, file_name, sample_file))
  @CONF: R_libs = susieR
  seed: $seed
  X: $X_sample
  Z: $PC_sample
  pve: 0.04
  Z_pve: 0
  n_signal: 1,2,3,4,5
  n_traits: 1
  file_name: file(pheno)
  sample_file: $sample_file
  $Y: res$Y
  $meta: res$meta
  $pheno_file: file_name
  $Zpve: Z_pve

sim_gaussian_s2(sim_gaussian):
  n_signal: 2

sim_gaussian_Z_s2(sim_gaussian):
  Z_pve: 0.01
  n_signal: 2
  
sim_gaussian_null(sim_gaussian):
  pve: 0
  n_signal: 0

sim_gaussian_large(sim_gaussian):
  n_signal: 200
