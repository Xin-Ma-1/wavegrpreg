#######################################################################
## calculate the optimization problem with a single lam
## Args:
##      X.all:     all X covariates
##      Y.all:     all Y responses
##      obs:       number of observations in each dataset
##      lam:       penality of theta's
##      maxit:     maximal number of iterations
##      tol:       tolerance level
## Returns:
##      coe:       fitted coefficients
##      gamma:     fitted gamma's
##      iteration: number of iterations
##      converge:  indicator if convergence is achieved
##      lam:       penality used
##      diff:      last step difference

## adapted from Li, Quefeng, et al.
## "Meta-analysis based variable selection for gene expression data."
## Biometrics 70.4 (2014): 872-880.
########################################################################
library(glmnet)

metalasso <- function(X.all, Y.all, obs, lam1, maxit, tol){
  ## starting and ending index of each dataset
  start.idx   <- cumsum(obs) + 1 - obs    # starting index of each dataset
  end.idx     <- cumsum(obs)              # ending index of each dataset
  M           <- length(start.idx)        # number of datasets
  p           <- ncol(X.all)              # number of covariates
  N           <- sum(obs)                 # total number of obserations
  gamma       <- rep(1, p)                # vector of gamma
  theta       <- rep(NA, M * p)           # vector of theta
  X.tha       <- matrix(NA, N, p)         # colMultiply(X.all, theta)
  beta.hat    <- vector("list", M)        # list of beta.hat vectors
  itr         <- 1
  m.diff      <- NULL                     # marginal error
  for (m in 1:M){
    beta.hat[[m]] <- rep(NA, p)
  }

  while(!(itr > maxit)){
    ## Iterate as: theta --> gamma --> theta

    beta.hat.old = beta.hat
    for (m in 1:M){
      ## In each dataset, fit Y.all ~ colMultiply(X.all, gamma*rho)
      theta.fit <- glmnet(t(t(X.all[start.idx[m]:end.idx[m], ])*gamma),
                          Y.all[start.idx[m]:end.idx[m]],
                          lambda = lam1,
                          standardize = FALSE)

      theta[((m - 1) * p + 1):(m * p)] <- as.vector(theta.fit$beta)

      ## adjust X.all by colMultiply(X.all, theta) for further usage
      X.tha[start.idx[m]:end.idx[m], ] <- t(t(X.all[start.idx[m]:end.idx[m], ]) *
                                              as.vector(theta.fit$beta))

      beta.hat[[m]] <- theta[((m - 1) * p + 1):(m * p)] * gamma
      ## calculate iteration difference
      if (itr == 1) {
        m.diff[m] <- max(abs(beta.hat[[m]]))
      }
      else {
        m.diff[m] <- max(abs(beta.hat[[m]] - beta.hat.old[[m]]))
      }
    }

    if(max(m.diff) < tol) break         # break iterations if diff < tol

    if(sum(sapply(beta.hat, function(a) sum(a!=0)))==0) break  # break if all beta's are zero

    itr <- itr + 1

    gamma.fit <- glmnet(X.tha, Y.all,
                        lambda = 1,
                        weights = rep(1 / obs, obs),
                        lower.limits=0,
                        standardize = FALSE)
    gamma <- as.vector(gamma.fit$beta)

    if(sum(gamma!=0)==0){
      for(m in 1:M){beta.hat[[m]] <- gamma}
      break       # break if all gamma's are zero, then all beta's are zero
    }
  }

  ## determine if convergence is achieved
  if (itr == 1) {
    iteration <- itr
    converge  <- FALSE
  }
  else {
    if (itr > maxit) {
      iteration <- itr - 1
      converge  <- FALSE
    }
    else {
      iteration <- itr
      converge  <- TRUE
    }
  }

  # get intercept terms
  coe0 = rep(NA,M)
  for(m in 1:M){
    Ym.hat = apply(X.all[start.idx[m]:end.idx[m],],1,function(a) sum(a*beta.hat[[m]]))
    coe0[m] = mean(Y.all[start.idx[m]:end.idx[m]]-Ym.hat)
  }

  return(list(coe         = beta.hat,
              coe0        = coe0,
              gamma       = gamma,
              iteration   = iteration,
              converge    = converge,
              lam         = lam1,
              diff        = max(m.diff)
  ))
}


######################################
# cross validation for metalasso
######################################
cv.metalasso <- function(X,Y,sizes,lams,maxit,tol,k=5){
  if(missing(sizes)){cat("Error: missing group sizes.","\n")}
  if(missing(maxit)){maxit=200}
  if(missing(tol)){tol=1e-8}
  if(missing(lams)){cat("missing shrinkage lambdas")}
  n.lam = length(lams)

  ngroup = length(sizes)
  cumsize = c(0,cumsum(sizes))

  # split the data into k lists proportional to group sizes
  split = vector("list",k)
  batch = round(sizes/k)
  cumbatch = c(0,cumsum(batch))
  for(i in 1:ngroup){
    order = sample(1:sizes[i])
    for(j in 1:(k-1)){
      split[[j]] = c(split[[j]],(order[((j-1)*batch[i]+1):(j*batch[i])]+cumsize[i]))
    }
    split[[k]] = c(split[[k]],(order[-(1:((k-1)*batch[i]))]+cumsize[i]))
  }

  MSE = rep(0,n.lam)
  for(b in 1:n.lam){
    lam = lams[b]

    for(i in 1:k){
      train = sort(unlist(split[-i]))
      test = sort(unlist(split[[i]]))
      if(i<k){
        fit = metalasso(X[train,], Y[train], obs=(sizes-batch), lam1=lam, maxit=maxit, tol=tol)
        Ytest.hat=rep(NA,sum(batch))
        for(j in 1:ngroup){
          Ym.hat = apply(X[test,][(cumbatch[j]+1):cumbatch[(j+1)],],1,function(c) sum(c*fit$coe[[j]])) + fit$coe0[j]
          Ytest.hat[(cumbatch[j]+1):cumbatch[(j+1)]] = Ym.hat
        }
        mse = sum((Y[test]-Ytest.hat)^2)


      }else{
        fit = metalasso(X[train,], Y[train], obs=((k-1)*batch),
                          lam1=lam, maxit=maxit, tol=tol)
        Ytest.hat=rep(NA,sum(sizes-(k-1)*batch))
        for(j in 1:ngroup){
          Ym.hat = apply(X[test,][((cumsize[j]-(k-1)*cumbatch[j])+1):(cumsize[(j+1)]-(k-1)*cumbatch[(j+1)]),],1,function(c) sum(c*fit$coe[[j]])) + fit$coe0[j]
          Ytest.hat[((cumsize[j]-(k-1)*cumbatch[j])+1):(cumsize[(j+1)]-(k-1)*cumbatch[(j+1)])] = Ym.hat
        }
        mse = sum((Y[test]-Ytest.hat)^2)
      }


      MSE[b] = MSE[b] + mse
    }

  }

  names(MSE) = paste0("lambda=",lams)

  return(list(MSE = MSE,
              best.lam = lams[which(MSE==min(MSE))]
  ))
}
