# Modules to provide data
# Real or simulated

# Module output
# =============
# $data: full data
# $sumstats: summary statistics

full_data: sim_utils.R + data_prepare.R
  tag: "full"
  dataset: Shell{head -200 ${data_file}}
  subset: NULL
  maf_thresh: 0
  GWASsample: ${GWAS_sample}
  REFsample: ${REF_sample}
  in_sample_id_file: file(.txt)
  snps_id_file: file(.txt)
  ld_sample_file: file(sample.ld)
  ld_sample_batch_file: file(sample.batch.ld)
  ld_sample_Z_file: file(sample.Z.ld)
  ld_ref_file: file(ref.ld)
  ld_ref_batch_file: file(ref.batch.ld)
  ld_ref_Z_file: file(ref.Z.ld)
  $X_sample: X.sample
  $X_sample_batch: X.sample.batch
  $PC_sample: PC.sample
  $X_sample_resid: X.sample.resid
  $N_sample: nrow(X.sample)
  $N_ref: nrow(X.ref)
  $maf: list(in_sample=maf.sample, ref_sample=maf.ref)
  $ld: list(in_sample=ld_sample_file, in_sample_batch = ld_sample_batch_file, in_sample_Z=ld_sample_Z_file, ref_sample=ld_ref_file, ref_sample_batch=ld_ref_batch_file, ref_sample_Z=ld_ref_Z_file)
  $r_Z_2dist: r.Z.2dist
  $r_ref_2dist: r.ref.2dist
  $r_Z_Mdist: r.Z.Mdist
  $r_ref_Mdist: r.ref.Mdist
  $seed: seed
  $X_file: dataset
  $sample_file: in_sample_id_file
  $snp_file: snps_id_file

lite_data(full_data):
  tag: "2k"
  subset: 2000

small_data(full_data):
  tag: "1k"
  maf_thresh: 0
  subset: 1000

tiny_data(full_data):
  tag: "300"
  maf_thresh: 0
  subset: 300

get_sumstats: regression.R
  @CONF: R_libs = (abind, data.table)
  method: 'lm', 'mixed'
  X: $X_sample_batch
  Z: $PC_sample
  Z_pve: ${Z_pve}
  Y: $Y
  X_file: $X_file
  Y_file: $pheno_file
  sample_file: $sample_file
  snp_file: $snp_file
  n_trait: ncol(Y)
  $sumstats: res

get_sumstats_lm(get_sumstats):
  method: 'lm'

summarize_ld: lib_regression_simulator.py + \
                regression_simulator.py + \
                Python(res = summarize_LD(X, ld_file['in_sample'], ld_plot))
  X: $X_sample
  ld_file: $ld
  $ld_plot: file(png)
  $top_idx: res
