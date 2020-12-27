# Modules to compute summary statistics

get_sumstats: regression.R
  @CONF: R_libs = (abind, data.table)
  method: 'lm', 'mixed'
  X: $X_sample_resid
  Z: $PC_sample
  Z_pve: $Zpve
  Y: $Y
  X_file: $X_file
  Y_file: $pheno_file
  sample_file: $sample_file
  snp_file: $snp_file
  n_trait: ncol(Y)
  $sumstats: res

get_sumstats_lm(get_sumstats):
  method: 'lm'

get_sumstats_pca(get_sumstats): R(Z = scale(Z, center=T, scale=T);
                                qrZ <- qr(Z);
                                Z = Z[, qrZ$pivot[1:qrZ$rank]];
                                qrZ.R = qr.R(qrZ)[1:qrZ$rank, 1:qrZ$rank];
                                qrZ.Q = qr.Q(qrZ)[, 1:qrZ$rank];
                                W = crossprod(qrZ.Q, X);
                                SZX = backsolve(qrZ.R, W);
                                X = X - Z %*% SZX;
                                r = cor(X);
                                write.table(r, ldfile, quote=F, col.names=F, row.names=F);
                                ld$in_sample_pca = ldfile) + regression.R
  X: $X_sample
  Z: $Z
  Z_pve: $Zpve
  ld: $ld
  ldfile: file(pca.ld)
  $ld: ld
  
  