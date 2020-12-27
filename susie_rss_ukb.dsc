#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/fit
%include modules/evaluate

DSC:
  define:
    method_susie: susie, susie_rss
  run:
    default: small_data * sim_gaussian * get_sumstats_lm * (method_susie * score_susie, finemap * score_finemap, finemapv3 * score_finemapv3, caviar * score_caviar)
    checkLD: small_data * sim_gaussian * get_sumstats_lm * method_susie * score_susie
    large_region: full_data * sim_gaussian * get_sumstats_lm * (susie * score_susie, method_susie_rss_simple * score_susie, method_finemap_simple * score_finemap, method_finemap_v3_simple * score_finemapv3, method_caviar_simple * score_caviar)
    mix_effect: small_data * sim_gaussian * get_sumstats * susie_rss * score_susie
    choose_weight: small_data * sim_gaussian_s2 * get_sumstats_lm * method_susie * score_susie
  exec_path: code
  global:
    data_file: data/ukb_genotypes_files.txt
    GWAS_sample: 10000
    REF_sample: 500
    Z_pve: 0.01
    neale_effect_size: "data/UKB_neale_effectsize.rds"
  output: output/rss_compare_ukb
