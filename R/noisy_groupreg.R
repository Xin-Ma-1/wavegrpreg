###############################################################
# group lasso for noisy input
###############################################################
noisy_groupreg <- function(X, Y, Sigma, sizes, z, lambda0, init=NA,
                           method, maxiter, tol0,
                           family, fil.num, min.level,
                           params=NA){
  par.glas = list(gamma=1e-4, a_min=1e-30, a_max=1e30,
                  sig1=0.1, sig2=0.9, Msave=10,
                  max1=5e2, max2=2e3, tol=1e-2)
  par.gbrg = list(gamma=1e-4, a_min=1e-30, a_max=1e30,
                  sig1=0.1, sig2=0.9, Msave=10,
                  max1=5e1, max2=2e2, tol=1e-1)
  if(sum(is.na(params))==0){par.glas = par.gbrg = params}

  ngroup = length(sizes)
  cumsize = c(0,cumsum(sizes))
  dims = dim(X[[1]])
  p = prod(dims)
  Y.center = rep(NA,length(Y))
  center.val = rep(NA,ngroup)
  for(m in 1:ngroup){
    ym = Y[(cumsize[m]+1):cumsize[(m+1)]]
    center.val[m] = mean(ym)
    Y.center[(cumsize[m]+1):cumsize[(m+1)]] = ym - center.val[m]
  }

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
  if(method=="glasso"){
    fit = SPG_glasso(C,Y.center,Sigma,sizes,z,lambda0,init_val,
                     par.glas$gamma, par.glas$a_min, par.glas$a_max,
                     par.glas$sig1, par.glas$sig2, par.glas$Msave,
                     par.glas$max1, par.glas$max2, par.glas$tol, L12proj)
  }else if(method=="gbridge"){
    fit = SPG_gbridge(C,Y.center,Sigma,sizes,z,lambda0,init_val,
                      maxiter,tol0,
                      par.gbrg$gamma, par.gbrg$a_min, par.gbrg$a_max,
                      par.gbrg$sig1, par.gbrg$sig2, par.gbrg$Msave,
                      par.gbrg$max1, par.gbrg$max2, par.gbrg$tol, L1proj)
  }


  # inverse wavelet transform to obtain the original coefficients beta
  eta = list()
  for(m in 1:ngroup){
    eta[[m]] = fit$eta[((m-1)*p+1):(m*p)]
  }
  if(length(dims)==2){
    beta = lapply(eta,reconstr.2D,dims=dims,fil.num=fil.num,family=family,
                  min.level=min.level)
  }else if(length(dims)==3){
    beta = lapply(eta,reconstr.3D,dims=dims,fil.num=fil.num,family=family,
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
cv.noisy_groupreg <- function(X, Y, Sigma, sizes, zvec, lamvec, init=NA,
                            method, maxiter=50, tol0=1e-2,
                            family, fil.num, min.level,
                            params=NA){
  n.z = length(zvec)
  n.lam = length(lamvec)
  ngroup = length(sizes)
  cumsize = c(0,cumsum(sizes))
  k = 5 # 5 fold CV

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

  MSE = matrix(0,n.z,n.lam)
  for(a in 1:n.z){
    for(b in 1:n.lam){
      z = zvec[a]
      lam = lamvec[b]

      message(cat("start iteration for z=",z,", lam=",lam))

      for(i in 1:k){
        train = sort(unlist(split[-i]))
        test = sort(unlist(split[[i]]))
        if(i<k){
          fit = noisy_groupreg(X[train], Y[train], Sigma, sizes=(sizes-batch), z, lam, init,
                               method, maxiter, tol0,
                               family, fil.num, min.level,
                               params = params)
          Ytest.hat=rep(NA,sum(batch))
          for(j in 1:ngroup){
            Ytest.hat[(cumbatch[j]+1):cumbatch[(j+1)]] =
              sapply(X[test][(cumbatch[j]+1):cumbatch[(j+1)]],
                     function(c) sum(c*fit$beta[[j]])) + fit$beta0[j]
          }
          mse = sum((Y[test]-Ytest.hat)^2)


        }else{
          fit = noisy_groupreg(X[train], Y[train], Sigma, sizes=((k-1)*batch), z, lam, init,
                               method, maxiter, tol0,
                               family, fil.num, min.level,
                               params = params)
          Ytest.hat=rep(NA,sum(sizes-(k-1)*batch))
          for(j in 1:ngroup){
            Ytest.hat[((cumsize[j]-(k-1)*cumbatch[j])+1):(cumsize[(j+1)]-(k-1)*cumbatch[(j+1)])] =
              sapply(X[test][((cumsize[j]-(k-1)*cumbatch[j])+1):(cumsize[(j+1)]-(k-1)*cumbatch[(j+1)])],
                     function(c) sum(c*fit$beta[[j]])) + fit$beta0[j]
          }
          mse = sum((Y[test]-Ytest.hat)^2)
        }

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
