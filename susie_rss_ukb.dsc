#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/fit
%include modules/evaluate

DSC:
  define:
    method_susie: susie, add_z_susierss*susie_rss
    method_susie_oracle: init_oracle*susie, add_z_susierss*init_susie_rss_oracle*susie_rss
    method_susie_lasso: init_lasso*susie, add_z_susierss*init_susie_rss_lasso*susie_rss
    method_susie_rss_simple: susie_rss_simple
    method_finemap: add_z * finemap
    method_finemap_simple: finemap_simple
    method_finemap_v3: add_z * finemapv3
    method_finemap_v3_simple: finemapv3_simple
    method_caviar: add_z * caviar
    method_caviar_simple: caviar_simple
  run:
    default: small_data * sim_gaussian * get_sumstats * (method_susie * score_susie, method_finemap * score_finemap, method_finemap_v3 * score_finemapv3, method_caviar * score_caviar)
    large_region: full_data * sim_gaussian * get_sumstats * (susie * score_susie, method_susie_rss_simple * score_susie, method_finemap_simple * score_finemap, method_finemap_v3_simple * score_finemapv3, method_caviar_simple * score_caviar)
    mix_effect: prepare_data * small_data
  exec_path: code
  global:
    data_file: data/ukb_genotypes_files.txt
    GWAS_sample: 10000
    REF_sample: 500
    neale_effect_size: "data/UKB_neale_effectsize.rds"
  output: output/rss_compare_ukb
