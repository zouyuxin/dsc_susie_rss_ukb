# Modules to provide data
# Real or simulated

# Module output
# =============
# $data: full data
# $sumstats: summary statistics

prepare_data: data_removeZ.R
  dataset: Shell{cat ${data_file}}
  XtX_full_file: file(full.XTX)
  $X: res$X
  $pos: res$pos
  $XtX_full: XtX_full_file
  
full_data: sim_utils.R + data_sim.R
  tag: "full"
  X: $X
  pos: $pos
  XtX_full: $XtX_full
  subset: NULL
  GWASsample: ${GWAS_sample}
  REFsample: ${REF_sample}
  ld_sample_file: file(sample.ld)
  ld_ref_file: file(ref.ld)
  $X_sample: X.sample
  $X_ref: X.ref
  $N_sample: nrow(X.sample)
  $N_ref: nrow(X.ref)
  $ld: list(in_sample=ld_sample_file, ref_sample=ld_ref_file)

lite_data(full_data):
  tag: "2k"
  subset: 2000

small_data(full_data):
  tag: "1k"
  subset: 1000

tiny_data(full_data):
  tag: "300"
  subset: 300

get_sumstats: regression.R + R(res = mm_regression(as.matrix(X), as.matrix(Y)))
  @CONF: R_libs = abind
  X: $X_sample
  Y: $Y
  $sumstats: res

summarize_ld: lib_regression_simulator.py + \
                regression_simulator.py + \
                Python(res = summarize_LD(X, ld_file['in_sample'], ld_plot))
  X: $X_sample
  ld_file: $ld
  $ld_plot: file(png)
  $top_idx: res
