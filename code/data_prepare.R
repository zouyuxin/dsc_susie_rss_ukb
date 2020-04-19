library(data.table)
library(Matrix)
library(readr)
library(digest)
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

# Get subset of X (individuals)
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


# Load covariates
# read phenotype files
pheno.file <- "genotype_dir/height.csv.gz"
pheno        <- suppressMessages(read_csv(pheno.file))
class(pheno) <- "data.frame"
pheno$sex = factor(pheno$sex)
pheno$assessment_centre = factor(pheno$assessment_centre)
pheno$genotype_measurement_batch = factor(pheno$genotype_measurement_batch)
pheno$age2 = pheno$age^2
# match individual order with genotype file
ind = fread(paste0(dataset, '.psam'))
match.idx = match(ind$IID, pheno$id)
pheno = pheno[match.idx,]

remove_Z = function(X, pheno){
  Z = model.matrix(~ pc_genetic1 + pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 +
                     pc_genetic6 + pc_genetic7 + pc_genetic8 + pc_genetic9 + pc_genetic10 +
                     pc_genetic11 + pc_genetic12 + pc_genetic13 + pc_genetic14 + pc_genetic15 +
                     pc_genetic16 + pc_genetic17 + pc_genetic18 + pc_genetic19 + pc_genetic20, 
                   data = pheno)
  # Remove intercept
  Z = Z[,-1]
  Z = scale(Z)
  
  # Center X
  X.c = scale(X, center=T, scale=FALSE)
  
  # Remove Z from X
  qrZ <- qr(Z) # Z = QR
  Z = Z[, qrZ$pivot[1:qrZ$rank]] # remove rank deficient columns
  qrZ.R = qr.R(qrZ)[1:qrZ$rank, 1:qrZ$rank]
  qrZ.Q = qr.Q(qrZ)[, 1:qrZ$rank]
  W = crossprod(qrZ.Q, X.c) # W = Q'X
  SZX = backsolve(qrZ.R, W) # (Z'Z)^{-1} Z'X = R^{-1} Q'X
  X.res = X.c - Z %*% SZX
  
  return(list(X=X.c, X.res=X.res, Z=Z, xtxdiag = colSums(X.c^2), W=W))
}

pheno.sample = pheno[in_sample, ]
X.sample.result = remove_Z(X.sample, pheno.sample)
X.sample = X.sample.result$X
X.sample.resid = X.sample.result$X.res
Z.sample = X.sample.result$Z

if(GWASsample == n){
  ld.matrix = as.matrix(fread(paste0(dataset,'.matrix')))
  ld.matrix = ld.matrix[choose.idx[overlap.idx[X.idx]], choose.idx[overlap.idx[X.idx]]]
  XtX = sqrt(X.sample.result$xtxdiag) * t(ld.matrix*sqrt(X.sample.result$xtxdiag)) - 
    crossprod(X.sample.result$W) # W'W = X'Q Q'X = X'Z(Z'Z)^{-1}Z'X
  r.sample = ld.matrix
  r.sample.Z = cov2cor(XtX)
}else{
  r.sample = cor(X.sample)
  r.sample.Z = cor(X.sample.result$X.res)
}

if (all(!is.na(X.ref))) {
  pheno.ref = pheno[ref_sample, ]
  X.ref.result = remove_Z(X.ref, pheno.ref)
  X.ref = X.ref.result$X
  Z.ref = X.ref.result$Z
  r.ref = cor(X.ref)
  r.ref.Z = cor(X.ref.result$X.res)
} else {
  r.ref = NA
}

if(all(!is.na(X.ref))){
  r.Z.2dist = Matrix::norm(r.sample - r.sample.Z, type='2')
  r.ref.2dist = Matrix::norm(r.sample - r.ref, type='2')
  r.Z.Mdist = max(abs(r.sample - r.sample.Z))
  r.ref.Mdist = max(abs(r.sample - r.ref))
}else{
  r.Z.2dist = r.ref.2dist = r.Z.Mdist = r.ref.Mdist = NA
}

write.table(ind[in_sample,c(1,2)], in_sample_id_file, quote=F, col.names=F, row.names=F)
write.table(gsub('_[A-Z]$','',colnames(X.sample)), snps_id_file, quote=F, col.names=F, row.names=F)
write.table(r.sample, ld_sample_file, quote=F, col.names=F, row.names=F)
write.table(r.sample.Z, ld_sample_Z_file, quote=F, col.names=F, row.names=F)
write.table(r.ref, ld_ref_file, quote=F, col.names=F, row.names=F)
write.table(r.ref.Z, ld_ref_Z_file, quote=F, col.names=F, row.names=F)

