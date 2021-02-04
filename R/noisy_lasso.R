###############################################################
# lasso for noisy input as in Loh & Wainwright
###############################################################
noisy_lasso <- function(X, Y, Sigma, z, lambda0, init=NA,
                        maxiter, tol0,
                        family, fil.num, min.level,
                        params=NA){
  par.las = list(gamma=1e-4, a_min=1e-30, a_max=1e30,
                 sig1=0.1, sig2=0.9, Msave=10,
                 max1=5e2, max2=2e3, tol=1e-2)
  if(sum(is.na(params))==0){par.las = params}

  dims = dim(X[[1]])
  p = prod(dims)
  center.val = mean(Y)
  Y.center = Y - center.val

  # wavelet transform on the design matrix W
  if(length(dims)==2){
    C = decomp.2D(X,family=family,fil.num=fil.num,min.level=min.level)
  }else if(length(dims)==3){
    C = decomp.3D(X,family=family,fil.num=fil.num,min.level=min.level)
  }

  # implement the group penalty algorithm
  init_val = NA
  if(sum(is.na(init))==0){
    if(length(dims)==2){
      init_val = c(t(decomp.2D(init,family=family,fil.num=fil.num,min.level=min.level)))
    }else if(length(dims)==3){
      init_val = c(t(decomp.3D(init,family=family,fil.num=fil.num,min.level=min.level)))
    }
  }
  fit = SPG_lasso(C, Y.center, Sigma,
                  z,lambda0,init_val,
                  par.las$gamma, par.las$a_min, par.las$a_max,
                  par.las$sig1, par.las$sig2, par.las$Msave,
                  par.las$max1, par.las$max2, par.las$tol, L1proj)

  # inverse wavelet transform to obtain the original coefficients beta
  if(length(dims)==2){
    beta = reconstr.2D(fit$eta,dims=dims,fil.num=fil.num,family=family,
                       min.level=min.level)
  }else if(length(dims)==3){
    beta = reconstr.3D(fit$eta,dims=dims,fil.num=fil.num,family=family,
                       min.level=min.level)
  }

  # return results
  return(list(beta       = beta,
              beta0      = center.val,
              z          = z,
              lambda0    = lambda0,
              family     = family,
              fil.num    = fil.num,
              min.level  = min.level,
              params     = params,
              iteration  = fit$iteration,
              diff       = fit$diff
  ))
}


###################################################
## cross validation based on prediction MSE
###################################################
cv.noisy_lasso <- function(X, Y, Sigma, zvec, lamvec, init=NA,
                           maxiter=50, tol0=1e-2,
                           family, fil.num, min.level,
                           params=NA){
  n.z = length(zvec)
  n.lam = length(lamvec)
  k = 5 # 5 fold CV

  # split the data into k lists proportional to group sizes
  split = vector("list",k)
  batch = round(length(Y)/k)
  order = sample(1:length(Y))
  for(j in 1:(k-1)){
    split[[j]] = order[((j-1)*batch+1):(j*batch)]
  }
  split[[k]] = order[-(1:((k-1)*batch))]

  MSE = matrix(0,n.z,n.lam)
  for(a in 1:n.z){
    for(b in 1:n.lam){
      z = zvec[a]
      lam = lamvec[b]

      message(cat("start iteration for z=",z,", lam=",lam))

      for(i in 1:k){
        train = sort(unlist(split[-i]))
        test = sort(unlist(split[[i]]))
        fit = noisy_lasso(X[train], Y[train], Sigma, z, lam, init,
                          maxiter, tol0,
                          family, fil.num, min.level,
                          params = params)
        Ytest.hat = sapply(X[test], function(c) sum(c*fit$beta)) + fit$beta0
        mse = sum((Y[test]-Ytest.hat)^2)

        message(cat("fold ",i,", mse=",mse))
        MSE[a,b] = MSE[a,b] + mse
      }
      message(cat("z=",z,", lam=",lam,", MSE=",MSE[a,b],"\n"))
    }
  }

  row.names(MSE) = paste0("z=",zvec)
  colnames(MSE) = paste0("lambda0=",lamvec)

  return(list(MSE = MSE,
              best.z = zvec[which(MSE==min(MSE),arr.ind=T)[1,1]],
              best.lam = lamvec[which(MSE==min(MSE),arr.ind=T)[1,2]],
              family = family,
              fil.num = fil.num,
              min.level = min.level
  ))

}
