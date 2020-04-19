#' @title sim_gaussian simulates a normal y from given data matrix X
#' @param X an n by p matrix, centered and scaled
#' @param pve a scalar percentage variance explained
#' @param effect_num a scalar number of true nonzero effects
#' @param effects effects from neale's lab
#' @return train_n a scalar number of trainning samples
#' @return sim_y an n vector simulated gaussian y
#' @return beta a p vector of effects
#' @return mean_corX mean of correlations of X (lower triangular entries of correlation matrix of X)
sim_gaussian = function(X, Z, pve, Z_pve = 0, effect_num, effects){
  n = dim(X)[1]
  p = dim(X)[2]

  beta.idx = sample(p, effect_num)
  beta = rep(0,p)
  beta.values = sample(effects, effect_num)
  beta[beta.idx] = beta.values
  
  if (effect_num==1){
    mean_corX = 1
  } else {
    effectX = X[,beta.idx]
    corX = cor(effectX)
    mean_corX = mean(abs(corX[lower.tri(corX)]))
  }
  if(effect_num==0){
    sigma = 1
    sim.y = rnorm(n, 0, 1)
    y = (sim.y - mean(sim.y))/sd(sim.y)
  } else {
    y_genetic = X %*% beta
    pheno_var = var(y_genetic) / pve
    resid_var = pheno_var - var(y_genetic)
    if(Z_pve > 0){ 
      Z_sigma2 = pheno_var*Z_pve/sum(apply(Z, 2, var))
      Z_coef = rnorm(ncol(Z), 0, sqrt(Z_sigma2))
      resid_var = pheno_var*(1 - Z_pve - pve)
    }else{
      Z_coef = numeric(ncol(Z))
    }
    epsilon = rnorm(n, mean = 0, sd = sqrt(resid_var))
    y = y_genetic + Z%*%Z_coef + epsilon
  }

  return(list(Y = y, sigma2 = resid_var,
              beta = beta, mean_corX = mean_corX))
}

# A wrapper for simulating multiple Y's
sim_gaussian_multiple = function(X, Z, pve, Z_pve=0, effect_num, effects, n_traits=1, file_name, sample_file) {
  meta = list(residual_variance = vector())
  Y = NULL
  sample_id = data.table::fread(sample_file)
  for (r in 1:n_traits) {
    res = sim_gaussian(X, Z, pve, Z_pve, effect_num, effects)
    if (is.null(Y)) Y = as.matrix(res$Y)
    else Y = cbind(Y, as.matrix(res$Y))
    if (is.null(meta$true_coef)) meta$true_coef = as.matrix(res$beta)
    else meta$true_coef = cbind(meta$true_coef, as.matrix(res$beta))
    meta$residual_variance[r] = res$sigma2
    write.table(cbind(sample_id, Y[,r]), paste0(file_name,r), quote=F, col.names=F, row.names=F)
  }
  return(list(Y=Y, meta=meta))
}
