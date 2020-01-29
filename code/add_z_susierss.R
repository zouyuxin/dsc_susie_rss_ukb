library(data.table);
z = sumstats$bhat/sumstats$shat;
r = as.matrix(fread(ld[[ld_method]]));
if (addz == TRUE && ld_method == 'ref_sample') {
  if (is.null(N_ref)) stop("Cannot use add_z out sample LD when N_ref is not available (NULL)")
  r = cov2cor(r*(N_ref-1) + tcrossprod(z));
  r = (r + t(r))/2;
}
