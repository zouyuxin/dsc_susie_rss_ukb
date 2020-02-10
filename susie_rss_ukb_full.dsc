#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/fit
%include modules/evaluate

DSC:
  define:
    method_susie: susie, susie_rss
    method_susie_simple: susie_simple, susie_rss_simple
  run:
    default: small_data * sim_gaussian * get_sumstats * (method_susie * score_susie, finemap * score_finemap, finemapv3 * score_finemapv3, caviar * score_caviar)
    large_region: full_data * sim_gaussian * get_sumstats * (method_susie_simple * score_susie, finemap_simple * score_finemap, finemapv3_simple * score_finemapv3, caviar_simple * score_caviar)
    mix_effect: prepare_data * small_data
  exec_path: code
  global:
    data_file: data/ukb_genotypes_files_top20.txt
    GWAS_sample: 274549
    REF_sample: 500
    neale_effect_size: "data/UKB_neale_effectsize.rds"
  output: output/rss_compare_ukb_full
