library(data.table)
library(Matrix)
library(readr)
# read genotypes
geno <- fread(paste0(dataset, '.raw.gz'),sep = "\t", header = TRUE, stringsAsFactors = FALSE)
class(geno) <- "data.frame"
# Extract the genotypes.
X <- as(as.matrix(geno[-(1:6)]),'dgCMatrix')

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

Z = model.matrix(~ sex + age + age2 + assessment_centre + genotype_measurement_batch +
                   pc_genetic1 + pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 +
                   pc_genetic6 + pc_genetic7 + pc_genetic8 + pc_genetic9 + pc_genetic10 +
                   pc_genetic11 + pc_genetic12 + pc_genetic13 + pc_genetic14 + pc_genetic15 +
                   pc_genetic16 + pc_genetic17 + pc_genetic18 + pc_genetic19 + pc_genetic20, data = pheno)

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
A   <- crossprod(Z)
SZX <- as.matrix(solve(A,t(Z) %*% X.c))
X   <- X.c - Z %*% SZX

xtxdiag = colSums(X.c^2)
ld.matrix = as.matrix(fread(paste0(dataset,'.matrix')))
# chol decomposition for (Z'Z)^(-1)
R = chol(solve(A)) # R'R = (Z'Z)^(-1)
W = R %*% crossprod(Z, X.c) # RZ'X
XtX = sqrt(xtxdiag) * t(ld.matrix*sqrt(xtxdiag)) - crossprod(W) # W'W = X'ZR'RZ'X = X'Z(Z'Z)^{-1}Z'X

write.table(XtX, XtX_full_file, quote=F, col.names=F, row.names=F)

snps = fread(paste0(dataset, '.pvar'))
signal_pos = gsub('^.*height.chr\\d*.', '', dataset)
pos = max(which(snps$POS <= as.integer(signal_pos)))

res = list(X = X, pos = pos)
