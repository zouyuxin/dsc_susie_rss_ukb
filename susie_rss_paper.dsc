#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/fit
%include modules/evaluate

DSC:
  define:
    method_susie: susie, susie_rss, susie_rss_add_z
    method_finemap: finemap, finemap_add_z
    method_finemap_v3: finemapv3, finemapv3_add_z
    method_caviar: caviar, caviar_add_z
  run:
    default: small_data * sim_gaussian * get_sumstats * (method_susie * score_susie, method_finemap * score_finemap, method_finemap_v3 * score_finemapv3, method_caviar * score_caviar)
    large_region: prepare_data * full_data * sim_gaussian * get_sumstats * (method_susie * score_susie, method_finemap * score_finemap, method_caviar * score_caviar)
    mix_effect: prepare_data * small_data
  exec_path: code
  global:
    data_file: data/ukb_genotypes_files.txt
    GWAS_sample: 10000
    REF_sample: 500
    neale_effect_size: "data/UKB_neale_effectsize.rds"
  output: output/rss_compare_ukb
