library(data.table)
library(Matrix)
library(readr)
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
if(all(!is.na(X.ref))){
  X.ref = X.ref[,X.idx]
}


# Remove covariates
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
  Z = model.matrix(~ sex + age + age2 + assessment_centre + genotype_measurement_batch +
                     pc_genetic1 + pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 +
                     pc_genetic6 + pc_genetic7 + pc_genetic8 + pc_genetic9 + pc_genetic10 +
                     pc_genetic11 + pc_genetic12 + pc_genetic13 + pc_genetic14 + pc_genetic15 +
                     pc_genetic16 + pc_genetic17 + pc_genetic18 + pc_genetic19 + pc_genetic20, 
                   data = pheno)
  # Remove intercept
  Z = Z[,-1]
  Z = scale(Z, center=T, scale=F)
  # standardize quantitative columns
  cols = which(colnames(Z) %in% c("age","pc_genetic1","pc_genetic2","pc_genetic3","pc_genetic4",
                                  "pc_genetic5","pc_genetic6","pc_genetic7","pc_genetic8","pc_genetic9", 
                                  "pc_genetic10","pc_genetic11","pc_genetic12","pc_genetic13","pc_genetic14",
                                  "pc_genetic15","pc_genetic16","pc_genetic17","pc_genetic18","pc_genetic19","pc_genetic20"))
  Z[,cols] = scale(Z[,cols])
  Z[,'age2'] = Z[,'age']^2
  
  # Center X
  X.c = scale(X, center=T, scale=FALSE)
  
  # Remove Z from X
  qrZ <- qr(Z) # Z = QR
  Z = Z[, qrZ$pivot[1:qrZ$rank]] # remove rank deficient columns
  qrZ.R = qr.R(qrZ)[1:qrZ$rank, 1:qrZ$rank]
  qrZ.Q = qr.Q(qrZ)[, 1:qrZ$rank]
  W = crossprod(qrZ.Q, X.c) # W = Q'X
  SZX = backsolve(qrZ.R, W) # (Z'Z)^{-1} Z'X = R^{-1} Q'X
  X = X.c - Z %*% SZX
  
  return(list(X=X, xtxdiag = colSums(X.c^2), W=W))
}

pheno.sample = pheno[in_sample, ]
X.sample.res = remove_Z(X.sample, pheno.sample)
X.sample = X.sample.res$X

if(GWASsample == n){
  ld.matrix = as.matrix(fread(paste0(dataset,'.matrix')))
  ld.matrix = ld.matrix[choose.idx[overlap.idx[X.idx]], choose.idx[overlap.idx[X.idx]]]
  XtX = sqrt(X.sample.res$xtxdiag) * t(ld.matrix*sqrt(X.sample.res$xtxdiag)) - 
    crossprod(X.sample.res$W) # W'W = X'Q Q'X = X'Z(Z'Z)^{-1}Z'X
  r.sample = cov2cor(XtX)
}else{
  r.sample = cor(X.sample)
}

if (all(!is.na(X.ref))) {
  pheno.ref = pheno[ref_sample, ]
  X.ref.res = remove_Z(X.ref, pheno.ref)
  X.ref = X.ref.res$X
  r.ref = cor(X.ref)
} else {
  r.ref = NA
}

write.table(r.sample, ld_sample_file, quote=F, col.names=F, row.names=F)
write.table(r.ref, ld_ref_file, quote=F, col.names=F, row.names=F)



