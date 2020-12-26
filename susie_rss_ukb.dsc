#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/get_sumstats
%include modules/fit
%include modules/evaluate

DSC:
  define:
    method_susie_full: susie, susie_rss, susie_rss_lambda
    method_susie_ldweight: susie, susie_rss_zldweight
  run:
    default: small_data * sim_gaussian * get_sumstats_lm * (method_susie_full * score_susie, finemap * score_finemap, caviar * score_caviar)
    mix_effect: small_data * sim_gaussian_Z_s2 * get_sumstats * susie_rss * score_susie
    choose_weight: small_data * sim_gaussian_s2 * get_sumstats_lm * method_susie_ldweight * score_susie
  exec_path: code
  global:
    data_file: data/ukb_genotypes_files.txt
    GWAS_sample: 10000
    REF_sample: 500
  output: output/rss_compare_ukb
