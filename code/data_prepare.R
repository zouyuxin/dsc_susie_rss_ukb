library(data.table)
library(Matrix)
library(readr)
library(digest)
library(rsvd)
# get unique seed for this genotype
str_to_int = function(x){
  h = digest(x, algo = "xxhash32")
  xx = strsplit(tolower(h), "")[[1L]]
  pos = match(xx, c(0L:9L, letters[1L:6L]))
  sum((pos - 1L) * 4^(rev(seq_along(xx) - 1)))
}
seed = str_to_int(dataset)

# read genotypes
geno <- fread(paste0(dataset, '.raw.gz'),sep = "\t", 
              header = TRUE, stringsAsFactors = FALSE)
class(geno) <- "data.frame"
# Extract the genotypes.
X <- as(as.matrix(geno[-(1:6)]),'dgCMatrix')
ind <- geno[(1:2)]
# Get subset of X (individuals)
set.seed(1)
n = nrow(X)
in_sample = sample(1:n, GWASsample)
X.sample = X[in_sample,]
if(REFsample > 0){
  if(GWASsample == n){ # random choose individuals
    ref_sample = sample(1:n, REFsample)
  }else{ # choose individuals different from samples
    ref_sample = sample(setdiff(1:n, in_sample), REFsample)
  }
  X.ref = X[ref_sample,]
}else{
  X.ref = NA
}
rm(X)

# Remove invariant SNPs
sample.idx = apply(X.sample, 2, var, na.rm=TRUE) != 0
if (all(!is.na(X.ref))) {
  ref.idx = apply(X.ref, 2, var, na.rm=TRUE) != 0
} else {
  ref.idx = 1
}
choose.idx = which(sample.idx * ref.idx == 1)
X.sample = X.sample[, choose.idx]
if (all(!is.na(X.ref))) X.ref = X.ref[, choose.idx]

# Compute MAF
maf.sample = apply(X.sample, 2, function(x) sum(x)/(2*length(x)))
maf.sample = pmin(maf.sample, 1-maf.sample)
if (all(!is.na(X.ref))) {
    maf.ref = apply(X.ref, 2, function(x) sum(x)/(2*length(x)))
    maf.ref = pmin(maf.ref, 1-maf.ref)
} else {
    maf.ref = NA
}
# Filter MAF (UKB genotypes already have MAF > 0.01)
if (maf_thresh > 0) {
    sample.idx = maf.sample > maf_thresh
    if (all(!is.na(X.ref))) {
        ref.idx = maf.ref > maf_thresh
    } else {
        ref.idx = 1
    }
    overlap.idx = which(sample.idx*ref.idx == 1)
    X.sample = X.sample[, overlap.idx]
    maf.sample = maf.sample[overlap.idx]
    if (all(!is.na(X.ref))) {
        X.ref = X.ref[, overlap.idx]
        maf.ref = maf.ref[overlap.idx]
    }
}else{
  overlap.idx = 1:length(choose.idx)
}

# subset SNPs
snps = fread(paste0(dataset, '.pvar'))
signal_pos = gsub('^.*height.chr\\d*.', '', dataset)
pos = max(which(snps$POS[choose.idx[overlap.idx]] <= as.integer(signal_pos)))
X.idx = get_genotype(subset, ncol(X.sample), pos)
X.sample = X.sample[, X.idx]
maf.sample = maf.sample[X.idx]
if(all(!is.na(X.ref))){
  X.ref = X.ref[,X.idx]
  maf.ref = maf.ref[X.idx]
}

X.sample = center_scale(X.sample)
X.ref = center_scale(X.ref)
## Get PCA for X

library(rsvd)
out.pca <- rpca(X.sample,k = 5,center = TRUE,
                scale = TRUE,retx = TRUE)
pcs = out.pca$x
colnames(pcs) <- paste0("PC",1:5)

# Load covariates
# read phenotype files

# remove_Z = function(X, PC){
# 
#   Z = scale(PC, center=T, scale=T)
#   qrZ <- qr(Z) # Z = QR
#   Z = Z[, qrZ$pivot[1:qrZ$rank]] # remove rank deficient columns
#   
#   # Remove Z from X
#   qrZ.R = qr.R(qrZ)[1:qrZ$rank, 1:qrZ$rank]
#   qrZ.Q = qr.Q(qrZ)[, 1:qrZ$rank]
#   W = crossprod(qrZ.Q, X) # W = Q'X
#   SZX = backsolve(qrZ.R, W) # (Z'Z)^{-1} Z'X = R^{-1} Q'X
#   X.res = X - Z %*% SZX
# 
#   return(list(X=X.res, PCs=Z, xtxdiag = colSums(X^2), W=W))
# }

# X.sample.result = remove_Z(X.sample, pcs)
# X.sample.resid = X.sample.result$X

if(GWASsample == n){
  ld.matrix = as.matrix(fread(paste0(dataset,'.matrix')))
  ld.matrix = ld.matrix[choose.idx[overlap.idx[X.idx]], choose.idx[overlap.idx[X.idx]]]
  # XtX = sqrt(X.sample.result$xtxdiag) * t(ld.matrix*sqrt(X.sample.result$xtxdiag)) - 
  #   crossprod(X.sample.result$W) # W'W = X'Q Q'X = X'Z(Z'Z)^{-1}Z'X
  r.sample = ld.matrix
  # r.sample.Z = cov2cor(XtX)
}else{
  r.sample = cor(X.sample)
  # r.sample.Z = cor(X.sample.resid)
}

if (all(!is.na(X.ref))) {
  r.ref = cor(X.ref)
} else {
  r.ref = NA
}

if(all(!is.na(X.ref))){
  r.ref.2dist = Matrix::norm(r.sample - r.ref, type='2')
  r.ref.Mdist = max(abs(r.sample - r.ref))
}else{
  r.ref.2dist = r.ref.Mdist = NA
}

write.table(ind[in_sample,], in_sample_id_file, quote=F, col.names=F, row.names=F)
write.table(gsub('_[A-Z]$','',colnames(X.sample)), snps_id_file, quote=F, col.names=F, row.names=F)
write.table(r.sample, ld_sample_file, quote=F, col.names=F, row.names=F)
write.table(r.ref, ld_ref_file, quote=F, col.names=F, row.names=F)
