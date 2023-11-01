library(knockoff)
create_Yknockoff_para <- function(Y, method=c("asdp","equi","sdp"), shrink=F) {
  # Estimate the mean vectorand covariance matrix

  mu = colMeans(Y)
  
  # Estimate the covariance matrix
  if(!shrink) {
    Sigma = cov(Y)
    # Verify that the covariance matrix is positive-definite
    p = nrow(matrix(Sigma))
    tol=1e-9
    if (p<500) {
      lambda_min = min(eigen(Sigma)$values)
    }
    else {
      oldw <- getOption("warn")
      options(warn = -1)
      lambda_min = RSpectra::eigs(Sigma, 1, which="SM", opts=list(retvec = FALSE, maxitr=100, tol))$values
      options(warn = oldw)
      if( length(lambda_min)==0 ) {
        # RSpectra::eigs did not converge. Using eigen instead."
        lambda_min = min(eigen(Sigma)$values)
      }
    }
    if(lambda_min>tol*10){
      shrink=TRUE
      
    }
    
    #if(!is_posdef(Sigma)) {
    #  shrink=TRUE
    # }
  }
  if(shrink) {
    if (!requireNamespace('corpcor', quietly=T))
      stop('corpcor is not installed', call.=F)
    Sigma = tryCatch({suppressWarnings(matrix(as.numeric(corpcor::cov.shrink(Y,verbose=F)), nrow=ncol(Y)))},
                     warning = function(w){}, error = function(e) {
                       stop("SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0")
                     }, finally = {})
  }
  
  method = match.arg(method)
  if ((nrow(Sigma)<=500) && method=="asdp") {
    method="sdp"
  }
  
  
  
  # diag_s vector of length r, containing the pre-computed covariances between the original variables and the knockoffs
  diag_s = diag(switch(match.arg(method),
                       'equi' = create.solve_equi(Sigma),
                       'sdp'  = create.solve_sdp(Sigma),
                       'asdp' = create.solve_asdp(Sigma)))
  
  if (is.null(dim(diag_s))) {
    diag_s = diag(diag_s,length(diag_s))
  }
  
  # Sample the Gaussian knockoffs
  SigmaInv_s = solve(Sigma,diag_s)
  mu_k = Y - sweep(Y,2,mu,"-") %*% SigmaInv_s
  Sigma_k = 2*diag_s - diag_s %*% SigmaInv_s
  result=list(mu_k=mu_k,Sigma_k=Sigma_k)
  #Y_k = mu_k + matrix(rnorm(ncol(Y)*nrow(Y)),nrow(Y)) %*% chol(Sigma_k)
  return(result)
}

generateMultiKnockoff = function(Y, mu, Sigma, n=100, seed=NULL){
  if(is.null(seed)){
    seed = seq(1, n, 1)
  }else{
    seed = seed
  }
  cholSigma = chol(Sigma) 
  nrow = dim(Y)[1]
  ncol = dim(Y)[2]
  result = vector(mode = "list", length = n)
  for(i in 1:n){
    set.seed(seed[i])
    randMtx = rnorm(nrow * ncol, 0, 1)
    randMtx = matrix(randMtx, nrow = nrow, byrow = TRUE)
    Yk = mu + randMtx %*% cholSigma
    result[[i]] = Yk
  }
  return(result)
}
