library(data.table);
z = sumstats$bhat / sumstats$shat;
if(add_z){
  r = as.matrix(fread(ld[[ld_method]]));
  if(ld_method == 'ref_sample'){
    if (is.null(N_ref)) stop("Cannot use add_z out sample LD when N_ref is not available (NULL)")
    r = cov2cor(r*(N_ref-1) + tcrossprod(z));
    r = (r + t(r))/2;
    write.table(r,ld_ref_z_file,quote=F,col.names=F,row.names=F);
    ld_file = ld_ref_z_file;
  }else{
    r = cov2cor(r*(N_sample-1) + tcrossprod(z));
    r = (r + t(r))/2;
    write.table(r,ld_sample_z_file,quote=F,col.names=F,row.names=F);
    ld_file = ld_sample_z_file;
  }
} else { ld_file = ld[[ld_method]] } 
