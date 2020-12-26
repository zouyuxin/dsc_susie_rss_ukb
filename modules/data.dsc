# Modules to provide data
# Real or simulated

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
  ld_ref_file: file(ref.ld)
  $X_sample: X.sample
  $PC_sample: pcs
  $N_sample: nrow(X.sample)
  $N_ref: nrow(X.ref)
  $maf: list(in_sample=maf.sample, ref_sample=maf.ref)
  $ld: list(in_sample=ld_sample_file, ref_sample=ld_ref_file)
  $r_ref_2dist: r.ref.2dist
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

summarize_ld: lib_regression_simulator.py + \
                regression_simulator.py + \
                Python(res = summarize_LD(X, ld_file['in_sample'], ld_plot))
  X: $X_sample
  ld_file: $ld
  $ld_plot: file(png)
  $top_idx: res
