#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/get_sumstats
%include modules/fit
%include modules/evaluate

DSC:
  define:
    method_susie_init: susie_init, susie_rss_lambda_init
    method_susie: susie, susie_rss_zldweight
  run:
    default: small_data * sim_gaussian * get_sumstats_lm * (method_susie_init * score_susie, finemap * score_finemap, finemapv3 * score_finemapv3, caviar * score_caviar) # with Z_pve = 0
    mix_effect: small_data * sim_gaussian_simple * get_sumstats_simple * susie_rss_mix * score_susie
    choose_weight: small_data * sim_gaussian_s2 * get_sumstats_lm * method_susie * score_susie # Z_pve 0.01
  exec_path: code
  global:
    data_file: data/ukb_genotypes_files.txt
    GWAS_sample: 10000
    REF_sample: 500
    n_dataset: 200
    neale_effect_size: "data/UKB_neale_effectsize.rds"
  output: output/rss_compare_ukb
